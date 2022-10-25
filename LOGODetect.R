suppressMessages(require(snowfall))
suppressMessages(require(data.table))
suppressMessages(require(optparse))
suppressMessages(require(BEDMatrix))
suppressMessages(require(XPASS))
options(stringsAsFactors = F)
options(datatable.fread.datatable = F)
warnings('off')



############################### input of LOGODetect
option_list = list(
  make_option('--sumstats', action = 'store', default = NA, type = 'character'),
  make_option('--n_gwas', action = 'store', default = NA, type = 'character'),
  make_option('--ref_dir', action = 'store', default = NA, type = 'character'),
  make_option('--pop', action = 'store', default = NA, type = 'character'),
  make_option('--ldsc_dir', action = 'store', default = NA, type = 'character'),
  make_option('--block_partition', action = 'store', default = NA, type = 'character'),
  make_option('--gc_snp', action = 'store', default = NA, type = 'character'),
  make_option('--out_dir', action = 'store', default = NA, type = 'character'),
  make_option('--n_cores', action = 'store', default = 20, type = 'numeric'),
  make_option('--target_pop', action = 'store', default = NA, type = 'character'),
  make_option('--n_topregion', action = 'store', default = NA, type = 'numeric'),
  make_option('--max_nsnps', action = 'store', default = 2000, type = 'integer'),
  make_option('--interval', action = 'store', default = 20, type = 'integer')
)
opt = parse_args(OptionParser(option_list=option_list))


sumstats = strsplit(opt$sumstats, ',')[[1]]
dat1 = data.frame(fread(sumstats[1]))
dat2 = data.frame(fread(sumstats[2]))
n_gwas = as.numeric(strsplit(opt$n_gwas, ',')[[1]])
n1 = n_gwas[1]
n2 = n_gwas[2]
ref_dir = opt$ref_dir
pop = strsplit(opt$pop, ',')[[1]]
if(length(pop) == 1){
  npop = 1
  pop1 = pop2 = pop
}else{
  npop = 2
  pop1 = pop[1]
  pop2 = pop[2]
}
if(npop == 1){
  ref = paste0(ref_dir, '/', pop, '/1000G_', pop, '_QC')
  ref_bim = data.frame(fread(paste0(ref, '.bim')))
  ref_bed = BEDMatrix(paste0(ref, '.bed'))
  colnames(ref_bed) = ref_bim$V2
  n_ref = nrow(ref_bed)
  ldscore = data.frame(fread(paste0(ref_dir, '/', pop, '/', pop, '.l2.ldscore')))
}
if(npop == 2){
  ref1 = paste0(ref_dir, '/', pop1, '/1000G_', pop1, '_QC')
  ref2 = paste0(ref_dir, '/', pop2, '/1000G_', pop2, '_QC')
  cov1 = paste0(ref_dir, '/', pop1, '/1000G_', pop1, '_QC_cov.txt')
  cov2 = paste0(ref_dir, '/', pop2, '/1000G_', pop2, '_QC_cov.txt')
  ref_bim = data.frame(fread(paste0(ref1, '.bim')))
  ref1_bed = BEDMatrix(paste0(ref1, '.bed'))
  ref2_bed = BEDMatrix(paste0(ref2, '.bed'))
  colnames(ref1_bed) = colnames(ref2_bed) = ref_bim$V2
  n_ref1 = nrow(ref1_bed)
  n_ref2 = nrow(ref2_bed)
  ldscore1 = data.frame(fread(paste0(ref_dir, '/', pop1, '/', pop1, '.l2.ldscore')))
  ldscore2 = data.frame(fread(paste0(ref_dir, '/', pop2, '/', pop2, '.l2.ldscore')))
}
block = read.table(opt$block_partition, header = T)
ldsc_dir = opt$ldsc_dir
if(!is.na(opt$gc_snp)){
  gc_snp = data.frame(fread(opt$gc_snp))[, 1]
}
out_dir = opt$out_dir
ncore = opt$n_cores
target_pop = opt$target_pop
n_topregion = opt$n_topregion
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}
Cn = opt$max_nsnps
inter = opt$interval

############################### QC
dat1 = dat1[complete.cases(dat1), ]
dat2 = dat2[complete.cases(dat2), ]
dat1$N = n1
dat2$N = n2
dat1$P[dat1$P < 1e-300] = 1e-300
dat2$P[dat2$P < 1e-300] = 1e-300
dat1$P[dat1$P > 1 - 1e-300] = 1 - 1e-300
dat2$P[dat2$P > 1 - 1e-300] = 1 - 1e-300
if('BETA' %in% colnames(dat1)){
  dat1$Z = -qnorm(dat1$P/2) * sign(dat1$BETA)
  dat2$Z = -qnorm(dat2$P/2) * sign(dat2$BETA)
}
if('OR' %in% colnames(dat1)){
  dat1$Z = -qnorm(dat1$P/2) * sign(log(dat1$OR))
  dat2$Z = -qnorm(dat2$P/2) * sign(log(dat2$OR))
}
dat1 = dat1[dat1$Z != Inf & dat1$Z != -Inf, ]
dat2 = dat2[dat2$Z != Inf & dat2$Z != -Inf, ]
dat1 = dat1[dat1$A1 %in% c('A', 'C', 'T', 'G') & dat1$A2 %in% c('A', 'C', 'T', 'G'), ]
dat2 = dat2[dat2$A1 %in% c('A', 'C', 'T', 'G') & dat2$A2 %in% c('A', 'C', 'T', 'G'), ]
dat1 = dat1[! ((dat1$A1=="G")&(dat1$A2=="C") | (dat1$A1=="C")&(dat1$A2=="G") | (dat1$A1=="A")&(dat1$A2=="T") | (dat1$A1=="T")&(dat1$A2=="A")), ]
dat2 = dat2[! ((dat2$A1=="G")&(dat2$A2=="C") | (dat2$A1=="C")&(dat2$A2=="G") | (dat2$A1=="A")&(dat2$A2=="T") | (dat2$A1=="T")&(dat2$A2=="A")), ]
dat1 = dat1[, c('CHR', 'SNP', 'BP', 'A1', 'A2', 'P', 'N', 'Z')]
dat2 = dat2[, c('CHR', 'SNP', 'BP', 'A1', 'A2', 'P', 'N', 'Z')]

### remove duplicated SNPs
dat1 = dat1[! dat1$SNP %in% dat1$SNP[duplicated(dat1$SNP)], ]
dat2 = dat2[! dat2$SNP %in% dat2$SNP[duplicated(dat2$SNP)], ]

### extract SNPs overlapped with reference SNPs
if(npop == 1){
  dat1_new = dat1[1, ][-1, ]
  dat2_new = dat2[1, ][-1, ]
  ref_bim_new = ref_bim[1, ][-1, ]
  ref_bed_new = NULL
  ldscore_new = ldscore[1, ][-1, ]
  for(chr in 1:22){
    f1 = dat1[dat1$CHR == chr, ]
    f2 = dat2[dat2$CHR == chr, ]
    ind = (ref_bim$V1 == chr)
    ref_bim_chr = ref_bim[ind, ]
    ref_bed_chr = ref_bed[, ind]
    ldscore_chr = ldscore[ind, ]
    
    f1 = f1[order(f1$BP), ]
    f1 = f1[! f1$BP %in% f1$BP[duplicated(f1$BP)], ]
    f1 = f1[f1$SNP %in% ref_bim_chr$V2 & f1$BP %in% ref_bim_chr$V4, ]
    temp = ref_bim_chr[ref_bim_chr$V2 %in% f1$SNP, ]
    pos = (f1$A1 == temp$V5)*(f1$A2 == temp$V6)
    neg = (f1$A1 == temp$V6)*(f1$A2 == temp$V5)
    flop = (neg == 1)
    f1$A1[flop] = temp$V5[flop]
    f1$A2[flop] = temp$V6[flop]
    f1$Z[flop] = -f1$Z[flop]
    ind = (pos + neg == 1)
    f1 = f1[ind, ]
    
    f2 = f2[order(f2$BP), ]
    f2 = f2[! f2$BP %in% f2$BP[duplicated(f2$BP)], ]
    f2 = f2[f2$SNP %in% ref_bim_chr$V2 & f2$BP %in% ref_bim_chr$V4, ]
    temp = ref_bim_chr[ref_bim_chr$V2 %in% f2$SNP, ]
    pos = (f2$A1 == temp$V5)*(f2$A2 == temp$V6)
    neg = (f2$A1 == temp$V6)*(f2$A2 == temp$V5)
    flop = (neg == 1)
    f2$A1[flop] = temp$V5[flop]
    f2$A2[flop] = temp$V6[flop]
    f2$Z[flop] = -f2$Z[flop]
    ind = (pos + neg == 1)
    f2 = f2[ind, ]
    
    overlap_snp = intersect(f1$SNP, f2$SNP)
    f1 = f1[f1$SNP %in% overlap_snp, ]
    f2 = f2[f2$SNP %in% overlap_snp, ]
    dat1_new = rbind(dat1_new, f1)
    dat2_new = rbind(dat2_new, f2)
    ind = ref_bim_chr$V2 %in% overlap_snp
    ref_bim_chr = ref_bim_chr[ind, ]
    ref_bed_chr = ref_bed_chr[, ind]
    ldscore_chr = ldscore_chr[ind, ]
    ref_bim_new = rbind(ref_bim_new, ref_bim_chr)
    ref_bed_new = cbind(ref_bed_new, ref_bed_chr)
    ldscore_new = rbind(ldscore_new, ldscore_chr)
  }
  dat1 = dat1_new
  dat2 = dat2_new
  ref_bim = ref_bim_new
  ref_bed = ref_bed_new
  ldscore = ldscore_new
  rm(dat1_new, dat2_new, ref_bim_new, ref_bed_new, ldscore_new)
}
if(npop == 2){
  dat1_new = dat1[1, ][-1, ]
  dat2_new = dat2[1, ][-1, ]
  ref_bim_new = ref_bim[1, ][-1, ]
  ref1_bed_new = NULL
  ref2_bed_new = NULL
  ldscore1_new = ldscore1[1, ][-1, ]
  ldscore2_new = ldscore2[1, ][-1, ]
  for(chr in 1:22){
    f1 = dat1[dat1$CHR == chr, ]
    f2 = dat2[dat2$CHR == chr, ]
    ind = (ref_bim$V1 == chr)
    ref_bim_chr = ref_bim[ind, ]
    ref1_bed_chr = ref1_bed[, ind]
    ref2_bed_chr = ref2_bed[, ind]
    ldscore1_chr = ldscore1[ind, ]
    ldscore2_chr = ldscore2[ind, ]
    
    f1 = f1[order(f1$BP), ]
    f1 = f1[! f1$BP %in% f1$BP[duplicated(f1$BP)], ]
    f1 = f1[f1$SNP %in% ref_bim_chr$V2 & f1$BP %in% ref_bim_chr$V4, ]
    temp = ref_bim_chr[ref_bim_chr$V2 %in% f1$SNP, ]
    pos = (f1$A1 == temp$V5)*(f1$A2 == temp$V6)
    neg = (f1$A1 == temp$V6)*(f1$A2 == temp$V5)
    flop = (neg == 1)
    f1$A1[flop] = temp$V5[flop]
    f1$A2[flop] = temp$V6[flop]
    f1$Z[flop] = -f1$Z[flop]
    ind = (pos + neg == 1)
    f1 = f1[ind, ]
    
    f2 = f2[order(f2$BP), ]
    f2 = f2[! f2$BP %in% f2$BP[duplicated(f2$BP)], ]
    f2 = f2[f2$SNP %in% ref_bim_chr$V2 & f2$BP %in% ref_bim_chr$V4, ]
    temp = ref_bim_chr[ref_bim_chr$V2 %in% f2$SNP, ]
    pos = (f2$A1 == temp$V5)*(f2$A2 == temp$V6)
    neg = (f2$A1 == temp$V6)*(f2$A2 == temp$V5)
    flop = (neg == 1)
    f2$A1[flop] = temp$V5[flop]
    f2$A2[flop] = temp$V6[flop]
    f2$Z[flop] = -f2$Z[flop]
    ind = (pos + neg == 1)
    f2 = f2[ind, ]
    
    overlap_snp = intersect(f1$SNP, f2$SNP)
    f1 = f1[f1$SNP %in% overlap_snp, ]
    f2 = f2[f2$SNP %in% overlap_snp, ]
    dat1_new = rbind(dat1_new, f1)
    dat2_new = rbind(dat2_new, f2)
    ind = ref_bim_chr$V2 %in% overlap_snp
    ref_bim_chr = ref_bim_chr[ind, ]
    ref1_bed_chr = ref1_bed_chr[, ind]
    ref2_bed_chr = ref2_bed_chr[, ind]
    ldscore1_chr = ldscore1_chr[ind, ]
    ldscore2_chr = ldscore2_chr[ind, ]
    ref_bim_new = rbind(ref_bim_new, ref_bim_chr)
    ref1_bed_new = cbind(ref1_bed_new, ref1_bed_chr)
    ref2_bed_new = cbind(ref2_bed_new, ref2_bed_chr)
    ldscore1_new = rbind(ldscore1_new, ldscore1_chr)
    ldscore2_new = rbind(ldscore2_new, ldscore2_chr)
  }
  dat1 = dat1_new
  dat2 = dat2_new
  ref_bim = ref_bim_new
  ref1_bed = ref1_bed_new
  ref2_bed = ref2_bed_new
  ldscore1 = ldscore1_new
  ldscore2 = ldscore2_new
  rm(dat1_new, dat2_new, ref_bim_new, ref1_bed_new, ref2_bed_new, ldscore1_new, ldscore2_new)
}



############################### caluculate heritability and genetic correlation
if(npop == 1){
  if(!dir.exists(paste0(out_dir, '/tmp_files'))){
    dir.create(paste0(out_dir, '/tmp_files'))
  }
  if(!dir.exists(paste0(out_dir, '/tmp_files/ldsc'))){
    dir.create(paste0(out_dir, '/tmp_files/ldsc'))
  }
  write.table(dat1, paste0(out_dir, '/tmp_files/dat1.txt'), row.names = F, col.names = T, quote = F)
  write.table(dat2, paste0(out_dir, '/tmp_files/dat2.txt'), row.names = F, col.names = T, quote = F)
  
  job = paste0('python ', ldsc_dir, '/munge_sumstats.py --sumstats ', out_dir, '/tmp_files/dat1.txt --chunksize 1000000 --out ', out_dir, '/tmp_files/ldsc/dat1_reformated --merge-alleles ', ref_dir, '/ldsc/w_hm3.snplist;
                python ', ldsc_dir, '/munge_sumstats.py --sumstats ', out_dir, '/tmp_files/dat2.txt --chunksize 1000000 --out ', out_dir, '/tmp_files/ldsc/dat2_reformated --merge-alleles ', ref_dir, '/ldsc/w_hm3.snplist;
                python ', ldsc_dir, '/ldsc.py --rg ', out_dir, '/tmp_files/ldsc/dat1_reformated.sumstats.gz,', out_dir, '/tmp_files/ldsc/dat2_reformated.sumstats.gz --ref-ld-chr ', ref_dir, '/ldsc/eur_w_ld_chr/ --w-ld-chr ', ref_dir, '/ldsc/eur_w_ld_chr/ --out ', out_dir, '/tmp_files/ldsc/ldsc_rg')
  system(job, ignore.stdout = T)
}
if(npop == 2){
  if(!dir.exists(paste0(out_dir, '/tmp_files'))){
    dir.create(paste0(out_dir, '/tmp_files'))
  }
  if(!dir.exists(paste0(out_dir, '/tmp_files/XPASS'))){
    dir.create(paste0(out_dir, '/tmp_files/XPASS'))
  }
  write.table(dat1, paste0(out_dir, '/tmp_files/dat1.txt'), row.names = F, col.names = T, quote = F)
  write.table(dat2, paste0(out_dir, '/tmp_files/dat2.txt'), row.names = F, col.names = T, quote = F)
  dat1_xpass = dat1[dat1$SNP %in% gc_snp, ]
  dat2_xpass = dat2[dat2$SNP %in% gc_snp, ]
  dat1_file = paste0(out_dir, '/tmp_files/XPASS/dat1_xpass.txt')
  dat2_file = paste0(out_dir, '/tmp_files/XPASS/dat2_xpass.txt')
  write.table(dat1_xpass, dat1_file, row.names = F, col.names = T, quote = F)
  write.table(dat2_xpass, dat2_file, row.names = F, col.names = T, quote = F)
  fit = XPASS(file_z1 = dat1_file, file_z2 = dat2_file,
              file_ref1 = ref1, file_ref2 = ref2,
              file_cov1 = cov1, file_cov2 = cov2,
              file_predGeno = NULL, compPRS = F,
              sd_method = "LD_block", compPosMean = F)
  xpass = data.frame(fit$H)
  write.table(xpass, paste0(out_dir, '/tmp_files/XPASS/XPASS_output.txt'), row.names = F, col.names = T, quote = F)
  rm(dat1_xpass, dat2_xpass, fit)
  gc()
}



############################### create tmp files in blocks
M = nrow(dat1)
if(!dir.exists(paste0(out_dir, '/tmp_files/block'))){
  dir.create(paste0(out_dir, '/tmp_files/block'))
}
write.table(M, paste0(out_dir, '/tmp_files/block/M.txt'), row.names = F, col.names = F, quote = F)
if(npop == 1){
  Read = file(paste0(out_dir, '/tmp_files/ldsc/ldsc_rg.log'), 'r')
  line = readLines(Read, n = 32)
  line = readLines(Read, n = 1)
  line = gsub("   ", " ", line)
  line = gsub("  ", " ", line)
  line = strsplit(line, ':')[[1]][2]
  h2_1 = as.numeric(strsplit(line, " ")[[1]][2])
  line = readLines(Read, n = 7)
  line = readLines(Read, n = 1)
  line = gsub("   ", " ", line)
  line = gsub("  ", " ", line)
  line = strsplit(line, ':')[[1]][2]
  h2_2 = as.numeric(strsplit(line, " ")[[1]][2])
  line = readLines(Read, n = 7)
  line = readLines(Read, n = 1)
  line = gsub("   ", " ", line)
  line = gsub("  ", " ", line)
  line = strsplit(line, ':')[[1]][2]
  gcov_total = as.numeric(strsplit(line, " ")[[1]][2])
  line = readLines(Read, n = 1)
  line = readLines(Read, n = 1)
  intercept = as.numeric(strsplit(line, " ")[[1]][2])
  close(Read)
  h2_snp_1 = h2_1/M
  h2_snp_2 = h2_2/M
  
  for(block.ind in 1:nrow(block)){
    chr = block$chr[block.ind]
    bp_start = block$start[block.ind]
    bp_stop = block$stop[block.ind]
    ind = (ref_bim$V1 == chr & ref_bim$V4 >= bp_start & ref_bim$V4 <= bp_stop)
    if(sum(ind) > 0){
      dat1_block = dat1[ind, ]
      dat2_block = dat2[ind, ]
      ref_bim_block = ref_bim[ind, ]
      ref_bed_block = ref_bed[, ind]
      ldsc_block = ldscore$L2[ind]
      ldsc_block = ldsc_block * (ldsc_block >= 1) + 1 * (ldsc_block < 1)
      write.table(dat1_block, paste0(out_dir, '/tmp_files/block/dat1_block_', block.ind, '.txt'), row.names = F, col.names = T, quote = F)
      write.table(dat2_block, paste0(out_dir, '/tmp_files/block/dat2_block_', block.ind, '.txt'), row.names = F, col.names = T, quote = F)
      write.table(ref_bim_block, paste0(out_dir, '/tmp_files/block/ref_bim_block_', block.ind, '.txt'), row.names = F, col.names = F, quote = F)
      write.table(ref_bed_block, paste0(out_dir, '/tmp_files/block/ref_bed_block_', block.ind, '.txt'), row.names = F, col.names = T, quote = F)
      write.table(ldsc_block, paste0(out_dir, '/tmp_files/block/ldsc_block_', block.ind, '.txt'), row.names = F, col.names = F, quote = F)
    }
  }
  rm(f1, f2, ind, dat1, dat2, ldscore, ldscore_chr, ref_bim, ref_bim_chr, ref_bed, ref_bed_chr)
}
if(npop == 2){
  h2_snp_1 = xpass$h1[1]/M
  h2_snp_2 = xpass$h2[1]/M
  for(block.ind in 1:nrow(block)){
    chr = block$chr[block.ind]
    bp_start = block$start[block.ind]
    bp_stop = block$stop[block.ind]
    ind = (ref_bim$V1 == chr & ref_bim$V4 >= bp_start & ref_bim$V4 <= bp_stop)
    if(sum(ind) > 0){
      dat1_block = dat1[ind, ]
      dat2_block = dat2[ind, ]
      ref_bim_block = ref_bim[ind, ]
      ref1_bed_block = ref1_bed[, ind]
      ref2_bed_block = ref2_bed[, ind]
      ldscore1_block = ldscore1[ind, ]
      ldscore2_block = ldscore2[ind, ]
      ldsc1 = n1 * h2_snp_1 * ldscore1_block$L2 + (1 - h2_snp_1 * M)
      ldsc2 = n2 * h2_snp_2 * ldscore2_block$L2 + (1 - h2_snp_2 * M)
      xldsc_block = abs(ldsc1 * ldsc2)
      write.table(dat1_block, paste0(out_dir, '/tmp_files/block/dat1_block_', block.ind, '.txt'), row.names = F, col.names = T, quote = F)
      write.table(dat2_block, paste0(out_dir, '/tmp_files/block/dat2_block_', block.ind, '.txt'), row.names = F, col.names = T, quote = F)
      write.table(ref_bim_block, paste0(out_dir, '/tmp_files/block/ref_bim_block_', block.ind, '.txt'), row.names = F, col.names = F, quote = F)
      write.table(ref1_bed_block, paste0(out_dir, '/tmp_files/block/ref1_bed_block_', block.ind, '.txt'), row.names = F, col.names = T, quote = F)
      write.table(ref2_bed_block, paste0(out_dir, '/tmp_files/block/ref2_bed_block_', block.ind, '.txt'), row.names = F, col.names = T, quote = F)
      write.table(xldsc_block, paste0(out_dir, '/tmp_files/block/xldsc_block_', block.ind, '.txt'), row.names = F, col.names = F, quote = F)
    }
  }
  rm(f1, f2, ind, dat1, dat2, ldscore1, ldscore2, ldscore1_chr, ldscore2_chr, ref_bim, ref_bim_chr, ref1_bed, ref1_bed_chr, ref2_bed, ref2_bed_chr)
}
gc()



############################### scan statistic
n_montecarlo = 5000
# Cn = 2000
# inter = 20
if(npop == 1){
  theta = c(0.5, 0.55, 0.6, 0.65, 0.7)
}
if(npop == 2){
  theta = c(0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75)
}
thre = 1000

scan = function(z1, z2, ldsc, Cn, inter, le, ri, theta){
  m = length(z1)
  z1_time_z2 = c(0,  cumsum(z1 * z2))
  ldscore_sum = c(0, cumsum(ldsc))
  if(m > inter){
    Qmax = -Inf
    window = seq(from = inter, to = floor(min(Cn, m)/inter)*inter, by = inter)
    for(I in window){
      j_count = ceiling((m-I+1)/inter)
      qq = rep(0, j_count)
      j = 1:j_count
      qq = (z1_time_z2[(j-1)*inter+1 +I] - z1_time_z2[(j-1)*inter+1])/(ldscore_sum[(j-1)*inter+1 +I] - ldscore_sum[(j-1)*inter+1])^theta
      qq = abs(qq)
      Q = max(qq)
      ind_Q = which.max(qq)
      if(Qmax < Q){
        Qmax = Q
        low = (ind_Q-1)*inter + 1
        high = (ind_Q-1)*inter + I
      }
    }
    re = list()
    re[[1]] = Qmax
    re[[2]] = c(low + le - 1, high + le - 1)
    return(re)
  }
  re = list()
  re[[1]] = z1_time_z2[m + 1]/(ldscore_sum[m + 1])^theta
  re[[2]] = c(le, ri)
  return(re)
}

svd.try = function(x, nu, nv){
  SVD = 'error'
  while(length(SVD) == 1){
    SVD = tryCatch(svd(x, nu, nv), error = function(e){ 'error' })
    x = x + diag(mean(abs(x))/10000, nrow = nrow(x))
  }
  return(SVD)
}

simulate_zscore_helper = function(segment.partition, geno.block, n_ref, thre){
  len = table(segment.partition$segment)
  segment.count = max(segment.partition$segment)
  
  if(segment.count == 1){
    ld_diag = ld2_diag = list()
    Sigma_p1 = Sigma_p2 = Sigma_V1 = Sigma_V2 = Sigma_D1 = Sigma_D2 = list()
    snp_1 = segment.partition$SNP[segment.partition$segment == 1]
    ld_diag[[1]] = cor(geno.block[, colnames(geno.block) %in% snp_1], use = 'pairwise.complete.obs')
    ld_diag[[1]][is.na(ld_diag[[1]])] = 0
    ld2_diag[[1]] = (n_ref-1)/(n_ref-2)*ld_diag[[1]]%*%ld_diag[[1]] - sum(len[1])/(n_ref-2)*ld_diag[[1]]
    
    SVD = svd.try(ld_diag[[1]], nu = 0, nv = nrow(ld_diag[[1]]))
    sv = SVD$d
    Sigma_p1[[1]] = sum(sv>sv[1]/thre)
    Sigma_D1[[1]] = sqrt(sv[1:Sigma_p1[[1]]])
    Sigma_V1[[1]] = SVD$v[,1:Sigma_p1[[1]]]
    SVD = svd.try(ld2_diag[[1]], nu = 0, nv = nrow(ld2_diag[[1]]))
    sv = SVD$d
    Sigma_p2[[1]] = sum(sv>sv[1]/thre)
    Sigma_D2[[1]] = sqrt(sv[1:Sigma_p2[[1]]])
    Sigma_V2[[1]] = SVD$v[,1:Sigma_p2[[1]]]
    
    return(list(p1 = Sigma_p1, p2 = Sigma_p2, V1 = Sigma_V1, V2 = Sigma_V2, D1 = Sigma_D1, D2 = Sigma_D2))
  }
  
  if(segment.count > 1){
    ld_diag = ld_offdiag = ld2_diag = ld2_offdiag = list()
    C1 = C2 = Sigma_p1 = Sigma_p2 = Sigma_V1 = Sigma_V2 = Sigma_D1 = Sigma_D2 = list()
    snp_1 = segment.partition$SNP[segment.partition$segment == 1]
    snp_2 = segment.partition$SNP[segment.partition$segment == 2]
    ld_diag[[1]] = cor(geno.block[, colnames(geno.block) %in% snp_1], use = 'pairwise.complete.obs')
    ld_diag[[2]] = cor(geno.block[, colnames(geno.block) %in% snp_2], use = 'pairwise.complete.obs')
    ld_offdiag[[1]] = cor(geno.block[, colnames(geno.block) %in% snp_1], geno.block[, colnames(geno.block) %in% snp_2], use = 'pairwise.complete.obs')
    ld_diag[[1]][is.na(ld_diag[[1]])] = ld_diag[[2]][is.na(ld_diag[[2]])] = ld_offdiag[[1]][is.na(ld_offdiag[[1]])] = 0
    ld2_diag[[1]] = (n_ref-1)/(n_ref-2)*( ld_diag[[1]]%*%ld_diag[[1]] + ld_offdiag[[1]]%*%t(ld_offdiag[[1]]) ) - sum(len[1:2])/(n_ref-2)*ld_diag[[1]]
    ld2_offdiag[[1]] = (n_ref-1)/(n_ref-2)*( ld_diag[[1]]%*%ld_offdiag[[1]] + ld_offdiag[[1]]%*%ld_diag[[2]] ) - sum(len[1:2])/(n_ref-2)*ld_offdiag[[1]]
    
    SVD = svd.try(ld_diag[[1]], nu = 0, nv = nrow(ld_diag[[1]]))
    sv = SVD$d
    Sigma_p1[[1]] = sum(sv>sv[1]/thre)
    Sigma_D1[[1]] = sqrt(sv[1:Sigma_p1[[1]]])
    Sigma_V1[[1]] = SVD$v[,1:Sigma_p1[[1]]]
    SVD = svd.try(ld2_diag[[1]], nu = 0, nv = nrow(ld2_diag[[1]]))
    sv = SVD$d
    Sigma_p2[[1]] = sum(sv>sv[1]/thre)
    Sigma_D2[[1]] = sqrt(sv[1:Sigma_p2[[1]]])
    Sigma_V2[[1]] = SVD$v[,1:Sigma_p2[[1]]]
    
    if(segment.count > 2){
      for(j in 2:(segment.count-1)){
        snp_1 = segment.partition$SNP[segment.partition$segment == j]
        snp_2 = segment.partition$SNP[segment.partition$segment == j + 1]
        ld_diag[[j+1]] = cor(geno.block[, colnames(geno.block) %in% snp_2], use = 'pairwise.complete.obs')
        ld_offdiag[[j]] = cor(geno.block[, colnames(geno.block) %in% snp_1], geno.block[, colnames(geno.block) %in% snp_2], use = 'pairwise.complete.obs')
        ld_diag[[j+1]][is.na(ld_diag[[j+1]])] = ld_offdiag[[j]][is.na(ld_offdiag[[j]])] = 0
        ld2_diag[[j]] = (n_ref-1)/(n_ref-2)*( ld_diag[[j]]%*%ld_diag[[j]] + t(ld_offdiag[[j-1]])%*%ld_offdiag[[j-1]] + ld_offdiag[[j]]%*%t(ld_offdiag[[j]]) ) - sum(len[(j-1):(j+1)])/(n_ref-2)*ld_diag[[j]]
        ld2_offdiag[[j]] = (n_ref-1)/(n_ref-2)*( ld_diag[[j]]%*%ld_offdiag[[j]] + ld_offdiag[[j]]%*%ld_diag[[j+1]] ) - sum(len[j:(j+1)])/(n_ref-2)*ld_offdiag[[j]]
        if(j > 2){
          ld_diag[[j-2]] = ld_offdiag[[j-2]] = ld2_diag[[j-2]] = ld2_offdiag[[j-2]] = 0 
        }
        
        svd1 = svd.try(ld_diag[[j-1]], nu = 0, nv = nrow(ld_diag[[j-1]]))
        sv1 = svd1$d
        p1 = sum(sv1>sv1[1]/thre)
        if(p1 == 1){
          D1 = 1/sv1[1:p1]
        }else{
          D1 = diag(1/sv1[1:p1])
        }
        V1 = as.matrix(svd1$v[,1:p1])
        C1[[j]] = t(ld_offdiag[[j-1]])%*%V1%*%D1%*%t(V1)
        Sigma1 = ld_diag[[j]] - C1[[j]]%*%ld_offdiag[[j-1]]
        svd2 = svd.try(ld2_diag[[j-1]], nu = 0, nv = nrow(ld2_diag[[j-1]]))
        sv2 = svd2$d
        p2 = sum(sv2>sv2[1]/thre)
        if(p2 == 1){
          D2 = 1/sv2[1:p2]
        }else{
          D2 = diag(1/sv2[1:p2])
        }
        V2 = as.matrix(svd2$v[,1:p2])
        C2[[j]] = t(ld2_offdiag[[j-1]])%*%V2%*%D2%*%t(V2)
        Sigma2 = ld2_diag[[j]] - C2[[j]]%*%ld2_offdiag[[j-1]]
        
        SVD = svd.try(Sigma1, nu = 0, nv = nrow(Sigma1))
        sv = SVD$d
        Sigma_p1[[j]] = sum(sv>sv[1]/thre)
        Sigma_D1[[j]] = sqrt(sv[1:Sigma_p1[[j]]])
        Sigma_V1[[j]] = SVD$v[,1:Sigma_p1[[j]]]
        SVD = svd.try(Sigma2, nu = 0, nv = nrow(Sigma2))
        sv = SVD$d
        Sigma_p2[[j]] = sum(sv>sv[1]/thre)
        Sigma_D2[[j]] = sqrt(sv[1:Sigma_p2[[j]]])
        Sigma_V2[[j]] = SVD$v[,1:Sigma_p2[[j]]]
      }
    }
    
    j = segment.count
    ld2_diag[[j]] = ld_diag[[j]]%*%ld_diag[[j]] + t(ld_offdiag[[j-1]])%*%ld_offdiag[[j-1]]
    svd1 = svd.try(ld_diag[[j-1]], nu = 0, nv = nrow(ld_diag[[j-1]]))
    sv1 = svd1$d
    p1 = sum(sv1>sv1[1]/thre)
    if(p1 == 1){
      D1 = 1/sv1[1:p1]
    }else{
      D1 = diag(1/sv1[1:p1])
    }
    V1 = as.matrix(svd1$v[,1:p1])
    C1[[j]] = t(ld_offdiag[[j-1]])%*%V1%*%D1%*%t(V1)
    Sigma1 = ld_diag[[j]] - C1[[j]]%*%ld_offdiag[[j-1]]
    svd2 = svd.try(ld2_diag[[j-1]], nu = 0, nv = nrow(ld2_diag[[j-1]]))
    sv2 = svd2$d
    p2 = sum(sv2>sv2[1]/thre)
    if(p2 == 1){
      D2 = 1/sv2[1:p2]
    }else{
      D2 = diag(1/sv2[1:p2])
    }
    V2 = as.matrix(svd2$v[,1:p2])
    C2[[j]] = t(ld2_offdiag[[j-1]])%*%V2%*%D2%*%t(V2)
    Sigma2 = ld2_diag[[j]] - C2[[j]]%*%ld2_offdiag[[j-1]]
    
    SVD = svd.try(Sigma1, nu = 0, nv = nrow(Sigma1))
    sv = SVD$d
    Sigma_p1[[j]] = sum(sv>sv[1]/thre)
    Sigma_D1[[j]] = sqrt(sv[1:Sigma_p1[[j]]])
    Sigma_V1[[j]] = SVD$v[,1:Sigma_p1[[j]]]
    SVD = svd.try(Sigma2, nu = 0, nv = nrow(Sigma2))
    sv = SVD$d
    Sigma_p2[[j]] = sum(sv>sv[1]/thre)
    Sigma_D2[[j]] = sqrt(sv[1:Sigma_p2[[j]]])
    Sigma_V2[[j]] = SVD$v[,1:Sigma_p2[[j]]]
    return(list(C1 = C1, C2 = C2, p1 = Sigma_p1, p2 = Sigma_p2, V1 = Sigma_V1, V2 = Sigma_V2, D1 = Sigma_D1, D2 = Sigma_D2))
  }
}

cal_qmax_1pop = function(bim, geno, n_ref, thre, n_montecarlo, n1, n2, h2_snp_1, h2_snp_2, intercept, M, ldsc, theta, Cn, inter){
  ### partition into segments
  bim$segment = as.numeric(as.factor(ceiling(bim$V4/1000000)))
  tab = which(table(bim$segment)>3000)
  while(length(tab)>0){
    ind = which(bim$segment == tab[1])
    K = ceiling(length(ind)/3000)
    quan = quantile(ind, seq(0, K)/K, type = 3)
    for(i in 2:K){
      bim$segment[(quan[i]+1):quan[i+1]] = bim$segment[(quan[i]+1):quan[i+1]] + 1/K*(i - 1)
    }
    bim$segment = as.numeric(as.factor(bim$segment))
    tab = which(table(bim$segment)>3000)
  }
  tab = which(table(bim$segment)<10)
  if(max(bim$segment) > 1){
    while(length(tab)>0){
      ind = which(bim$segment == tab[1])
      if(min(ind) == 1){
        bim$segment[ind] = bim$segment[max(ind) + 1]
      }else{
        bim$segment[ind] = bim$segment[min(ind) - 1]
      }
      bim$segment = as.numeric(as.factor(bim$segment))
      tab = which(table(bim$segment)<10)
    }
  }
  segment.partition = bim[, c(2, 7)]
  colnames(segment.partition) = c('SNP', 'segment')
  segment.count = max(segment.partition$segment)
  
  ### simulate zscore vectors under the null
  pop.helper = simulate_zscore_helper(segment.partition, geno, n_ref, thre)
  
  X = matrix(rnorm(n_montecarlo*3*pop.helper$p1[[1]]), nrow = pop.helper$p1[[1]], ncol = n_montecarlo*3)
  pop_ss = t(pop.helper$V1[[1]]%*%(X*pop.helper$D1[[1]]))
  X = matrix(rnorm(n_montecarlo*2*pop.helper$p2[[1]]), nrow = pop.helper$p2[[1]], ncol = n_montecarlo*2)
  pop_tt = t(pop.helper$V2[[1]]%*%(X*pop.helper$D2[[1]]))
  pop_SS = pop_ss
  pop_TT = pop_tt
  if(segment.count > 1){
    for(j in 2:segment.count){
      print(j)
      X = matrix(rnorm(n_montecarlo*3*pop.helper$p1[[j]]), nrow = pop.helper$p1[[j]], ncol = n_montecarlo*3)
      pop_ss_new = t(pop.helper$V1[[j]]%*%(X*pop.helper$D1[[j]])) + pop_ss%*%t(pop.helper$C1[[j]])
      X = matrix(rnorm(n_montecarlo*2*pop.helper$p2[[j]]), nrow = pop.helper$p2[[j]], ncol = n_montecarlo*2)
      pop_tt_new = t(pop.helper$V2[[j]]%*%(X*pop.helper$D2[[j]])) + pop_tt%*%t(pop.helper$C2[[j]])
      print(j)
      pop_SS = cbind(pop_SS, pop_ss_new)
      pop_TT = cbind(pop_TT, pop_tt_new)
      pop_ss = pop_ss_new
      pop_tt = pop_tt_new
    }
  }
  
  S1 = pop_SS[1:n_montecarlo, ]
  S2 = pop_SS[(n_montecarlo + 1):(2*n_montecarlo), ]
  S3 = pop_SS[(2*n_montecarlo + 1):(3*n_montecarlo), ]
  T1 = pop_TT[1:n_montecarlo, ]
  T2 = pop_TT[(n_montecarlo + 1):(2*n_montecarlo), ]
  Z1 = sqrt(n1*h2_snp_1)*T1 + sqrt(1 - h2_snp_1*M - abs(intercept))*S1 + sqrt(abs(intercept))*S3
  Z2 = sqrt(n2*h2_snp_2)*T2 + sqrt(1 - h2_snp_2*M - abs(intercept))*S2 + sqrt(abs(intercept))*sign(intercept)*S3
  
  sd1 = sd2 = rep(0, n_montecarlo)
  Qmax = matrix(0, nrow = n_montecarlo, ncol = length(theta))
  for(i in 1:n_montecarlo){
    z1 = Z1[i, ]
    z2 = Z2[i, ]
    m = length(z1)
    sd1[i] = sqrt(mean(z1^2))
    sd2[i] = sqrt(mean(z2^2))
    for(j in 1:length(theta)){
      re = scan(z1, z2, ldsc, Cn, inter, 1, m, theta[j])
      Qmax[i, j] = as.numeric(re[[1]])
    }
  }
  return(list(Qmax = Qmax, sd1 = sd1, sd2 = sd2))
}

cal_qmax_2pop = function(bim, geno1, geno2, n_ref1, n_ref2, thre, n_montecarlo, n1, n2, h2_snp_1, h2_snp_2, M, ldsc, theta, Cn, inter){
  ### partition into segments
  bim$segment = as.numeric(as.factor(ceiling(bim$V4/1000000)))
  tab = which(table(bim$segment)>3000)
  while(length(tab)>0){
    ind = which(bim$segment == tab[1])
    K = ceiling(length(ind)/3000)
    quan = quantile(ind, seq(0, K)/K, type = 3)
    for(i in 2:K){
      bim$segment[(quan[i]+1):quan[i+1]] = bim$segment[(quan[i]+1):quan[i+1]] + 1/K*(i - 1)
    }
    bim$segment = as.numeric(as.factor(bim$segment))
    tab = which(table(bim$segment)>3000)
  }
  tab = which(table(bim$segment)<10)
  if(max(bim$segment) > 1){
    while(length(tab)>0){
      ind = which(bim$segment == tab[1])
      if(min(ind) == 1){
        bim$segment[ind] = bim$segment[max(ind) + 1]
      }else{
        bim$segment[ind] = bim$segment[min(ind) - 1]
      }
      bim$segment = as.numeric(as.factor(bim$segment))
      tab = which(table(bim$segment)<10)
    }
  }
  segment.partition = bim[, c(2, 7)]
  colnames(segment.partition) = c('SNP', 'segment')
  segment.count = max(segment.partition$segment)
  
  ### simulate zscore vectors under the null
  pop1.helper = simulate_zscore_helper(segment.partition, geno1, n_ref1, thre)
  pop2.helper = simulate_zscore_helper(segment.partition, geno2, n_ref2, thre)
  
  X = matrix(rnorm(n_montecarlo*pop1.helper$p1[[1]]), nrow = pop1.helper$p1[[1]], ncol = n_montecarlo)
  pop1_ss = t(pop1.helper$V1[[1]]%*%(X*pop1.helper$D1[[1]]))
  X = matrix(rnorm(n_montecarlo*pop1.helper$p2[[1]]), nrow = pop1.helper$p2[[1]], ncol = n_montecarlo)
  pop1_tt = t(pop1.helper$V2[[1]]%*%(X*pop1.helper$D2[[1]]))
  pop1_SS = pop1_ss
  pop1_TT = pop1_tt
  X = matrix(rnorm(n_montecarlo*pop2.helper$p1[[1]]), nrow = pop2.helper$p1[[1]], ncol = n_montecarlo)
  pop2_ss = t(pop2.helper$V1[[1]]%*%(X*pop2.helper$D1[[1]]))
  X = matrix(rnorm(n_montecarlo*pop2.helper$p2[[1]]), nrow = pop2.helper$p2[[1]], ncol = n_montecarlo)
  pop2_tt = t(pop2.helper$V2[[1]]%*%(X*pop2.helper$D2[[1]]))
  pop2_SS = pop2_ss
  pop2_TT = pop2_tt
  
  if(segment.count > 1){
    for(j in 2:segment.count){
      X = matrix(rnorm(n_montecarlo*pop1.helper$p1[[j]]), nrow = pop1.helper$p1[[j]], ncol = n_montecarlo)
      pop1_ss_new = t(pop1.helper$V1[[j]]%*%(X*pop1.helper$D1[[j]])) + pop1_ss%*%t(pop1.helper$C1[[j]])
      X = matrix(rnorm(n_montecarlo*pop1.helper$p2[[j]]), nrow = pop1.helper$p2[[j]], ncol = n_montecarlo)
      pop1_tt_new = t(pop1.helper$V2[[j]]%*%(X*pop1.helper$D2[[j]])) + pop1_tt%*%t(pop1.helper$C2[[j]])
      pop1_SS = cbind(pop1_SS, pop1_ss_new)
      pop1_TT = cbind(pop1_TT, pop1_tt_new)
      pop1_ss = pop1_ss_new
      pop1_tt = pop1_tt_new
      
      X = matrix(rnorm(n_montecarlo*pop2.helper$p1[[j]]), nrow = pop2.helper$p1[[j]], ncol = n_montecarlo)
      pop2_ss_new = t(pop2.helper$V1[[j]]%*%(X*pop2.helper$D1[[j]])) + pop2_ss%*%t(pop2.helper$C1[[j]])
      X = matrix(rnorm(n_montecarlo*pop2.helper$p2[[j]]), nrow = pop2.helper$p2[[j]], ncol = n_montecarlo)
      pop2_tt_new = t(pop2.helper$V2[[j]]%*%(X*pop2.helper$D2[[j]])) + pop2_tt%*%t(pop2.helper$C2[[j]])
      pop2_SS = cbind(pop2_SS, pop2_ss_new)
      pop2_TT = cbind(pop2_TT, pop2_tt_new)
      pop2_ss = pop2_ss_new
      pop2_tt = pop2_tt_new
    }
  }
  
  Z1 = sqrt(n1*h2_snp_1)*pop1_TT + sqrt(1 - h2_snp_1*M)*pop1_SS
  Z2 = sqrt(n2*h2_snp_2)*pop2_TT + sqrt(1 - h2_snp_2*M)*pop2_SS
  sd1 = sd2 = rep(0, n_montecarlo)
  Qmax = matrix(0, nrow = n_montecarlo, ncol = length(theta))
  for(i in 1:n_montecarlo){
    z1 = Z1[i, ]
    z2 = Z2[i, ]
    m = length(z1)
    sd1[i] = sqrt(mean(z1^2))
    sd2[i] = sqrt(mean(z2^2))
    for(j in 1:length(theta)){
      re = scan(z1, z2, ldsc, Cn, inter, 1, m, theta[j])
      Qmax[i, j] = as.numeric(re[[1]])
    }
  }
  return(list(Qmax = Qmax, sd1 = sd1, sd2 = sd2))
}

apply.fun.1pop = function(block.ind){
  library(data.table)
  chr = block$chr[block.ind]
  bp_start = block$start[block.ind]
  bp_stop = block$stop[block.ind]
  if(file.exists(paste0(out_dir, '/tmp_files/block/ref_bim_block_', block.ind, '.txt'))){
    ref_bim_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/ref_bim_block_', block.ind, '.txt')))
    ref_bed_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/ref_bed_block_', block.ind, '.txt')))
    ldsc_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/ldsc_block_', block.ind, '.txt')))[, 1]
    
    result = cal_qmax_1pop(bim = ref_bim_block, geno = ref_bed_block, 
                           n_ref, thre, n_montecarlo, 
                           n1, n2, h2_snp_1, h2_snp_2, intercept, M, 
                           ldsc = ldsc_block, theta, Cn, inter)
    write.table(result$Qmax, paste0(out_dir, '/tmp_files/block/Qmax_block', block.ind, '.txt'), row.names = F, col.names = F, quote = F)
    write.table(result$sd1, paste0(out_dir, '/tmp_files/block/sd1_block', block.ind, '.txt'), row.names = F, col.names = F, quote = F)
    write.table(result$sd2, paste0(out_dir, '/tmp_files/block/sd2_block', block.ind, '.txt'), row.names = F, col.names = F, quote = F)
  }
}

apply.fun.2pop = function(block.ind){
  library(data.table)
  chr = block$chr[block.ind]
  bp_start = block$start[block.ind]
  bp_stop = block$stop[block.ind]
  if(file.exists(paste0(out_dir, '/tmp_files/block/ref_bim_block_', block.ind, '.txt'))){
    ref_bim_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/ref_bim_block_', block.ind, '.txt')))
    ref1_bed_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/ref1_bed_block_', block.ind, '.txt')))
    ref2_bed_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/ref2_bed_block_', block.ind, '.txt')))
    xldsc_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/xldsc_block_', block.ind, '.txt')))[, 1]
    
    result = cal_qmax_2pop(bim = ref_bim_block, geno1 = ref1_bed_block, geno2 = ref2_bed_block, 
                      n_ref1, n_ref2, thre, n_montecarlo, 
                      n1, n2, h2_snp_1, h2_snp_2, M, 
                      ldsc = xldsc_block, theta, Cn, inter)
    write.table(result$Qmax, paste0(out_dir, '/tmp_files/block/Qmax_block', block.ind, '.txt'), row.names = F, col.names = F, quote = F)
    write.table(result$sd1, paste0(out_dir, '/tmp_files/block/sd1_block', block.ind, '.txt'), row.names = F, col.names = F, quote = F)
    write.table(result$sd2, paste0(out_dir, '/tmp_files/block/sd2_block', block.ind, '.txt'), row.names = F, col.names = F, quote = F)
  }
}

sfInit(parallel = TRUE, cpus = ncore)
sfLibrary(snowfall)
if(npop == 1){
  sfExport('block', 'n_ref', 'thre', 'n_montecarlo', 'n1', 'n2', 'h2_snp_1', 'h2_snp_2', 'M', 'theta', 'Cn', 'inter', 'out_dir', 'intercept')
  sfExport('scan', 'svd.try', 'simulate_zscore_helper', 'cal_qmax_1pop')
  Result = sfLapply(1:nrow(block), apply.fun.1pop)
}
if(npop == 2){
  sfExport('block', 'n_ref1', 'n_ref2', 'thre', 'n_montecarlo', 'n1', 'n2', 'h2_snp_1', 'h2_snp_2', 'M', 'theta', 'Cn', 'inter', 'out_dir')
  sfExport('scan', 'svd.try', 'simulate_zscore_helper', 'cal_qmax_2pop')
  Result = sfLapply(1:nrow(block), apply.fun.2pop)
}
sfStop()



############################### scan statistic on real GWAS data
whole_scan = function(z1, z2, ldsc, Cn, inter, thre, le, ri, theta){
  re = scan(z1, z2, ldsc, Cn, inter, le, ri, theta)
  if(re[[1]] <= thre){
    return(list())
  }
  else{
    a = re[[2]][1]
    b = re[[2]][2]
    if(a == le & b == ri){return(re)}
    if(a == le & b < ri){
      return( c( re , whole_scan(z1[(b+2-le):(ri+1-le)], z2[(b+2-le):(ri+1-le)], ldsc[(b+2-le):(ri+1-le)], Cn, inter, thre, b+1, ri, theta) ) )
    }
    if(a > le & b == ri){
      return( c( whole_scan(z1[1:(a-le)], z2[1:(a-le)], ldsc[1:(a-le)], Cn, inter, thre, le, a-1, theta), re) )
    }
    return( c( whole_scan(z1[1:(a-le)], z2[1:(a-le)], ldsc[1:(a-le)], Cn, inter, thre, le, a-1, theta), re, whole_scan(z1[(b+2-le):(ri+1-le)], z2[(b+2-le):(ri+1-le)], ldsc[(b+2-le):(ri+1-le)], Cn, inter, thre, b+1, ri, theta) ) )
  }
}
dat1 = data.frame(fread(paste0(out_dir, '/tmp_files/dat1.txt')))
dat2 = data.frame(fread(paste0(out_dir, '/tmp_files/dat2.txt')))

Result = list()
for(j in 1:length(theta)){
  Result[[j]] = as.data.frame(matrix(0, nrow = 0, ncol = 7))
  colnames(Result[[j]]) = c('chr', 'begin_pos', 'stop_pos', 'count', 'stat', 'pval', 'block')
}
n.block = 0
if(npop == 1){
  for(block.ind in 1:nrow(block)){
    if(file.exists(paste0(out_dir, '/tmp_files/block/Qmax_block', block.ind, '.txt'))){
      n.block = n.block + 1
      dat1_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/dat1_block_', block.ind, '.txt')))
      dat2_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/dat2_block_', block.ind, '.txt')))
      ldsc_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/ldsc_block_', block.ind, '.txt')))[, 1]
      z1_block = dat1_block$Z
      z2_block = dat2_block$Z
      sd1_data = sqrt(mean(z1_block^2))
      sd2_data = sqrt(mean(z2_block^2))
      sd1_generated = fread(paste0(out_dir, '/tmp_files/block/sd1_block', block.ind, '.txt'))[, 1]
      sd2_generated = fread(paste0(out_dir, '/tmp_files/block/sd2_block', block.ind, '.txt'))[, 1]
      
      Qmax = as.matrix(read.table(paste0(out_dir, '/tmp_files/block/Qmax_block', block.ind, '.txt')))
      q = rep(0, length(theta))
      Qmax_scaled = matrix(0, nrow = n_montecarlo, ncol = length(theta))
      for(j in 1:length(theta)){
        Qmax_scaled[, j] = Qmax[, j]*(sd1_data * sd2_data)/(sd1_generated * sd2_generated)
        q[j] = quantile(Qmax_scaled[, j], 0.95)
      }
      
      for(j in 1:length(theta)){
        result = whole_scan(z1_block, z2_block, ldsc_block, Cn, inter, q[j], 1, length(z1_block), theta[j])
        detect = c()
        if(length(result)>0){
          for(i in 1:(length(result)/2)){
            ind.begin = result[[2*i]][1]
            ind.stop = result[[2*i]][2]
            count = ind.stop - ind.begin + 1
            stat = sum(z1_block[ind.begin:ind.stop]*z2_block[ind.begin:ind.stop])/(sum(ldsc_block[ind.begin:ind.stop]))^theta[j]
            pval = (sum(abs(stat) < Qmax_scaled[, j]) + 1)/(n_montecarlo + 1)
            chr = dat1_block$CHR[1]
            begin_pos = dat1_block$BP[ind.begin]
            stop_pos = dat1_block$BP[ind.stop]
            detect_new = c(chr, begin_pos, stop_pos, count, stat, pval, block.ind)
            detect = rbind(detect, detect_new)
          }
          detect = data.frame(detect, row.names = NULL)
          colnames(detect) = c('chr', 'begin_pos', 'stop_pos', 'count', 'stat', 'pval', 'block')
          Result[[j]] = rbind(Result[[j]], detect)
        }
      }
    }
  }
}
if(npop == 2){
  for(block.ind in 1:nrow(block)){
    if(file.exists(paste0(out_dir, '/tmp_files/block/Qmax_block', block.ind, '.txt'))){
      n.block = n.block + 1
      dat1_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/dat1_block_', block.ind, '.txt')))
      dat2_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/dat2_block_', block.ind, '.txt')))
      xldsc_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/xldsc_block_', block.ind, '.txt')))[, 1]
      z1_block = dat1_block$Z
      z2_block = dat2_block$Z
      sd1_data = sqrt(mean(z1_block^2))
      sd2_data = sqrt(mean(z2_block^2))
      sd1_generated = fread(paste0(out_dir, '/tmp_files/block/sd1_block', block.ind, '.txt'))[, 1]
      sd2_generated = fread(paste0(out_dir, '/tmp_files/block/sd2_block', block.ind, '.txt'))[, 1]
      
      Qmax = as.matrix(read.table(paste0(out_dir, '/tmp_files/block/Qmax_block', block.ind, '.txt')))
      q = rep(0, length(theta))
      Qmax_scaled = matrix(0, nrow = n_montecarlo, ncol = length(theta))
      for(j in 1:length(theta)){
        Qmax_scaled[, j] = Qmax[, j]*(sd1_data * sd2_data)/(sd1_generated * sd2_generated)
        q[j] = quantile(Qmax_scaled[, j], 0.95)
      }
      
      for(j in 1:length(theta)){
        result = whole_scan(z1_block, z2_block, xldsc_block, Cn, inter, q[j], 1, length(z1_block), theta[j])
        detect = c()
        if(length(result)>0){
          for(i in 1:(length(result)/2)){
            ind.begin = result[[2*i]][1]
            ind.stop = result[[2*i]][2]
            count = ind.stop - ind.begin + 1
            stat = sum(z1_block[ind.begin:ind.stop]*z2_block[ind.begin:ind.stop])/(sum(xldsc_block[ind.begin:ind.stop]))^theta[j]
            pval = (sum(abs(stat) < Qmax_scaled[, j]) + 1)/(n_montecarlo + 1)
            chr = dat1_block$CHR[1]
            begin_pos = dat1_block$BP[ind.begin]
            stop_pos = dat1_block$BP[ind.stop]
            detect_new = c(chr, begin_pos, stop_pos, count, stat, pval, block.ind)
            detect = rbind(detect, detect_new)
          }
          detect = data.frame(detect, row.names = NULL)
          colnames(detect) = c('chr', 'begin_pos', 'stop_pos', 'count', 'stat', 'pval', 'block')
          Result[[j]] = rbind(Result[[j]], detect)
        }
      }
    }
  }
}

for(j in 1:length(theta)){
  n0 = nrow(Result[[j]])
  if(n0 > 0){
    Result[[j]] = Result[[j]][order(abs(Result[[j]]$stat), decreasing = T),]
    Result[[j]] = Result[[j]][order(Result[[j]]$pval),]
    Result[[j]]$testing_count = n.block
    if(n0>1){
      for(i in 2:n0){
        temp = Result[[j]][1:(i-1), ]
        testing_add = (sum(temp$block == Result[[j]]$block[i]))>0
        Result[[j]]$testing_count[i:n0] = Result[[j]]$testing_count[i:n0] + testing_add
      }
    }
    p_val = Result[[j]]$pval
    BH = max(which(p_val<=0.05/Result[[j]]$testing_count*c(1:n0)))
    Result[[j]]$qval = Result[[j]]$pval * Result[[j]]$testing_count / c(1:n0)
    if(BH>0){
      Result[[j]] = Result[[j]][1:BH, c('chr', 'begin_pos', 'stop_pos', 'stat', 'pval', 'qval')]
      for(i in BH:1){
        Result[[j]]$qval[i] = min(Result[[j]]$qval[i:BH])
      }
    }else{
      Result[[j]] = Result[[j]][1, ][-1, ]
    }
  }
}



############################### calculate genetic covariane in aggregated regions
if(npop == 1){
  if(!dir.exists(paste0(out_dir, '/tmp_files/S-LDSC'))){
    dir.create(paste0(out_dir, '/tmp_files/S-LDSC'))
  }
  
  gcov = rep(0, length(theta))
  for(j in 1:length(theta)){
    if(nrow(Result[[j]])>0){
      re = Result[[j]][, c('chr', 'begin_pos', 'stop_pos')]
      re = re[order(re$begin_pos), ]
      re = re[order(re$chr), ]
      for(ch in 1:22){
        annot = data.frame(fread(paste0(ref_dir, '/ldsc/1000g_eur/1000G_EUR_QC_chr', ch, '.bim')))
        annot = annot[, 1:4]
        colnames(annot) = c('CHR', 'SNP', 'CM', 'BP')
        annot$C0 = 1
        annot$C1 = 0
        re_new = re[re$chr == ch, ]
        if(nrow(re_new) > 0){
          for(w in 1:nrow(re_new)){
            annot$C0[re_new$begin_pos[w]<=annot$BP & annot$BP<=re_new$stop_pos[w]] = 0
            annot$C1[re_new$begin_pos[w]<=annot$BP & annot$BP<=re_new$stop_pos[w]] = 1
          }
        }
        write.table(annot, paste0(out_dir, '/tmp_files/S-LDSC/LOGODetect_theta_', theta[j], '_chr', ch, '.annot'), col.names = T, row.names = F, quote = F)
      }
      
      ### S-LDSC
      for(ch in 1:22){
        job = paste0('python ', ldsc_dir, '/ldsc.py --l2 --bfile ', ref_dir, '/ldsc/1000g_eur/1000G_EUR_QC_chr', ch, ' --ld-wind-cm 1 --annot ', out_dir, '/tmp_files/S-LDSC/LOGODetect_theta_', theta[j], '_chr', ch, '.annot --out ', out_dir, '/tmp_files/S-LDSC/LOGODetect_theta_', theta[j], '_chr', ch, ' --print-snps ', ref_dir, '/ldsc/hapmap3_snps/hm.', ch, '.snp;')
        system(job, ignore.stdout = T)
      }
      job = paste0('python ', ldsc_dir, '/ldsc.py --rg ', out_dir, '/tmp_files/ldsc/dat1_reformated.sumstats.gz,', out_dir, '/tmp_files/ldsc/dat2_reformated.sumstats.gz --ref-ld-chr ', out_dir, '/tmp_files/S-LDSC/LOGODetect_theta_', theta[j], '_chr --w-ld-chr ', ref_dir, '/ldsc/eur_w_ld_chr/ --out ', out_dir, '/tmp_files/S-LDSC/LOGODetect_theta_', theta[j])
      system(job, ignore.stdout = T)
    }
  }
  
  for(j in 1:length(theta)){
    Read = file(paste0(out_dir, '/tmp_files/S-LDSC/LOGODetect_theta_', theta[j], '.log'), 'r')
    line = readLines(Read, n = 66)
    line = readLines(Read, n = 1)
    line = gsub("   ", " ", line)
    line = gsub("  ", " ", line)
    if(length(line)>0){
      line = strsplit(line, ':')[[1]][2]
      gcov[j] = abs(as.numeric(strsplit(line, ' ')[[1]][3]))
    }
    close(Read)
  }
}

if(npop == 2){
  xpass_theta = list()
  for(j in 1:length(theta)){
    if(nrow(Result[[j]])>0){
      re = Result[[j]][, c('chr', 'begin_pos', 'stop_pos')]
      re = re[order(re$begin_pos), ]
      re = re[order(re$chr), ]
      dat1_LOGO = dat1[1, ][-1, ]
      dat2_LOGO = dat2[1, ][-1, ]
      for(ch in 1:22){
        temp = re[re$chr == ch, ]
        dat1_chr = dat1[dat1$CHR == ch, ]
        dat2_chr = dat2[dat2$CHR == ch, ]
        dat1_LOGO_new = dat1_chr[1, ][-1, ]
        dat2_LOGO_new = dat2_chr[1, ][-1, ]
        if(nrow(temp)>0){
          for(k in 1:nrow(temp)){
            begin_pos = temp$begin_pos[k]
            stop_pos = temp$stop_pos[k]
            ind = dat1_chr$BP >= begin_pos & dat1_chr$BP <= stop_pos
            dat1_LOGO_new = rbind(dat1_LOGO_new, dat1_chr[ind, ])
            dat2_LOGO_new = rbind(dat2_LOGO_new, dat2_chr[ind, ])
          }
        }
        dat1_LOGO = rbind(dat1_LOGO, dat1_LOGO_new)
        dat2_LOGO = rbind(dat2_LOGO, dat2_LOGO_new)
      }
      write.table(dat1_LOGO, paste0(out_dir, '/tmp_files/XPASS/dat1_LOGODetect_theta_', theta[j], '.txt'), row.names = F, col.names = T, quote = F)
      write.table(dat2_LOGO, paste0(out_dir, '/tmp_files/XPASS/dat2_LOGODetect_theta_', theta[j], '.txt'), row.names = F, col.names = T, quote = F)
      
      dat1_LOGO = paste0(out_dir, '/tmp_files/XPASS/dat1_LOGODetect_theta_', theta[j], '.txt')
      dat2_LOGO = paste0(out_dir, '/tmp_files/XPASS/dat2_LOGODetect_theta_', theta[j], '.txt')
      fit = XPASS(file_z1 = dat1_LOGO, file_z2 = dat2_LOGO,
                  file_ref1 = ref1, file_ref2 = ref2,
                  file_cov1 = cov1, file_cov2 = cov2,
                  file_predGeno = NULL, compPRS = F,
                  sd_method = "LD_block", compPosMean = F)
      xpass_theta[[j]] = fit$H
      xpass_theta[[j]] = data.frame(xpass_theta[[j]])
    }else{
      xpass_theta[[j]] = NA
    }
  }
  
  gcov_total = xpass$h12[1]
  gcov = rep(0, length(theta))
  for(j in 1:length(theta)){
    if(class(xpass_theta[[j]]) == 'data.frame'){
      gcov[j] = xpass_theta[[j]]$h12[1]
    }
  }
}



############################### select theta and obtain output regions
if(sum(gcov != 0) > 0){
  j = which.max(gcov/gcov_total)
  theta_selected = theta[j]
  
  ### merge regions that are <100KB away
  re = Result[[j]]
  re = re[order(re$begin_pos), ]
  re = re[order(re$chr), ]
  re_merged = re[1, ][-1, ]
  for(ch in 1:22){
    temp = re[re$chr == ch, ]
    if(nrow(temp)>0){
      re_merged_chr = temp[1, ]
      k = 2
      while(k <= nrow(temp)){
        if(temp$begin_pos[k] - re_merged_chr$stop_pos[nrow(re_merged_chr)] < 100000){
          re_merged_chr$stop_pos[nrow(re_merged_chr)] = temp$stop_pos[k]
          re_merged_chr$stat[nrow(re_merged_chr)] = max(re_merged_chr$stat[nrow(re_merged_chr)], temp$stat[k])
          re_merged_chr$pval[nrow(re_merged_chr)] = min(re_merged_chr$pval[nrow(re_merged_chr)], temp$pval[k])
          re_merged_chr$qval[nrow(re_merged_chr)] = min(re_merged_chr$qval[nrow(re_merged_chr)], temp$qval[k])
          k = k + 1
        }else{
          re_merged_chr = rbind(re_merged_chr, temp[k, ])
          k = k + 1
        }
      }
      re_merged = rbind(re_merged, re_merged_chr)
    }
  }
  write.table(re_merged, paste0(out_dir, '/LOGODetect_regions.txt'), col.names = T, row.names = F, quote = F)
}else{
  write.table('No significant region is identified.', paste0(out_dir, '/LOGODetect_regions.txt'), col.names = T, row.names = F, quote = F)
}



### create annotations
if(!is.na(n_topregion)){
  topregion = as.data.frame(matrix(0, nrow = 0, ncol = 6))
  colnames(topregion) = c('chr', 'begin_pos', 'stop_pos', 'stat', 'pval', 'block')
  for(block.ind in 1:nrow(block)){
    if(file.exists(paste0(out_dir, '/tmp_files/block/Qmax_block', block.ind, '.txt'))){
      dat1_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/dat1_block_', block.ind, '.txt')))
      dat2_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/dat2_block_', block.ind, '.txt')))
      xldsc_block = data.frame(fread(paste0(out_dir, '/tmp_files/block/xldsc_block_', block.ind, '.txt')))[, 1]
      z1_block = dat1_block$Z
      z2_block = dat2_block$Z
      sd1_data = sqrt(mean(z1_block^2))
      sd2_data = sqrt(mean(z2_block^2))
      sd1_generated = fread(paste0(out_dir, '/tmp_files/block/sd1_block', block.ind, '.txt'))[, 1]
      sd2_generated = fread(paste0(out_dir, '/tmp_files/block/sd2_block', block.ind, '.txt'))[, 1]
      Qmax = as.matrix(read.table(paste0(out_dir, '/tmp_files/block/Qmax_block', block.ind, '.txt')))
      Qmax_scaled = Qmax[, j]*(sd1_data * sd2_data)/(sd1_generated * sd2_generated)
      q = min(Qmax_scaled)
      result = whole_scan(z1_block, z2_block, xldsc_block, Cn, inter, q, 1, length(z1_block), theta[j])
      detect = c()
      if(length(result)>0){
        for(i in 1:(length(result)/2)){
          ind.begin = result[[2*i]][1]
          ind.stop = result[[2*i]][2]
          stat = sum(z1_block[ind.begin:ind.stop]*z2_block[ind.begin:ind.stop])/(sum(xldsc_block[ind.begin:ind.stop]))^theta[j]
          pval = (sum(abs(stat) < Qmax_scaled) + 1)/(n_montecarlo + 1)
          chr = dat1_block$CHR[1]
          begin_pos = dat1_block$BP[ind.begin]
          stop_pos = dat1_block$BP[ind.stop]
          detect_new = c(chr, begin_pos, stop_pos, stat, pval, block.ind)
          detect = rbind(detect, detect_new)
        }
        detect = data.frame(detect, row.names = NULL)
        colnames(detect) = c('chr', 'begin_pos', 'stop_pos', 'stat', 'pval', 'block')
        topregion = rbind(topregion, detect)
      }
    }
  }
  write.table(topregion, paste0(out_dir, '/tmp_files/LOGODetect_topregions.txt'), col.names = T, row.names = F, quote = F)
  topregion = topregion[order(abs(topregion$stat), decreasing = T), ]
  topregion = topregion[order(topregion$pval), ]
  topregion = topregion[1:min(n_topregion, nrow(topregion)), ]
  topregion = topregion[topregion$stat * gcov_total>0, ]
  
  if(pop1 == target_pop){
    annot = data.frame(fread(sumstats[1]))
    annot$Anno = 1
    write.table(annot[, c('CHR', 'SNP', 'A1', 'A2', 'Anno')], paste0(out_dir, '/annot_', pop1, '.txt'), col.names = T, row.names = F, quote = F)
    
    annot = data.frame(fread(sumstats[2]))
    annot$Anno = 0
    annot_new = annot[1, ][-1, ]
    for(ch in 1:22){
      temp = topregion[topregion$chr == ch, ]
      if(nrow(temp)>0){
        annot_new_chr = annot[annot$CHR == ch, ]
        for(i in 1:nrow(temp)){
          annot_new_chr$Anno[annot_new_chr$BP >= temp$begin_pos[i] & annot_new_chr$BP <= temp$stop_pos[i]] = 1
        }
        annot_new = rbind(annot_new, annot_new_chr)
      }
    }
    annot = annot_new[, c('CHR', 'SNP', 'A1', 'A2', 'Anno')]
    write.table(annot, paste0(out_dir, '/annot_', pop2, '.txt'), col.names = T, row.names = F, quote = F)
  }
  if(pop2 == target_pop){
    annot = data.frame(fread(sumstats[1]))
    annot$Anno = 0
    annot_new = annot[1, ][-1, ]
    for(ch in 1:22){
      temp = topregion[topregion$chr == ch, ]
      if(nrow(temp)>0){
        annot_new_chr = annot[annot$CHR == ch, ]
        for(i in 1:nrow(temp)){
          annot_new_chr$Anno[annot_new_chr$BP >= temp$begin_pos[i] & annot_new_chr$BP <= temp$stop_pos[i]] = 1
        }
        annot_new = rbind(annot_new, annot_new_chr)
      }
    }
    annot = annot_new[, c('CHR', 'SNP', 'A1', 'A2', 'Anno')]
    write.table(annot, paste0(out_dir, '/annot_', pop1, '.txt'), col.names = T, row.names = F, quote = F)
    
    annot = data.frame(fread(sumstats[2]))
    annot$Anno = 1
    write.table(annot[, c('CHR', 'SNP', 'A1', 'A2', 'Anno')], paste0(out_dir, '/annot_', pop2, '.txt'), col.names = T, row.names = F, quote = F)
  }
  
}

unlink(paste0(out_dir, '/tmp_files'), recursive = T)
