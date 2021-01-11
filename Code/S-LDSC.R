options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)
dir_out = as.character(args[1])
if(!dir.exists(paste0(dir_out, '/Temp/S-LDSC'))){
  dir.create(paste0(dir_out, '/Temp/S-LDSC'))
}

theta = c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7)
for(j in 1:length(theta)){
  if(file.exists(paste0(dir_out, '/Temp/Result_different_theta/Result_theta_', theta[j], '.txt'))){
    ### Create annotation file
    re = read.table(paste0(dir_out, '/Temp/Result_different_theta/Result_theta_', theta[j], '.txt'), header = T)
    re = re[, c('chr', 'begin_pos', 'stop_pos')]
    re$begin_pos = re$begin_pos*1000000
    re$stop_pos = re$stop_pos*1000000
    for(ch in 1:22){
      annot = read.table(paste0('./Data/1000G_EUR_Phase3_plink/1000G.EUR.QC.', ch, '.bim'))
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
      write.table(annot, paste0(dir_out, '/Temp/S-LDSC/LOGODetect_theta_', theta[j], '_chr', ch, '.annot'), col.names = T, row.names = F, quote = F)
    }
    
    ### S-LDSC
    job = paste0('source activate ldsc; 
    for ch in $(seq 22)
    do
        ./Data/ldsc/ldsc.py --l2 --bfile ./Data/1000G_EUR_Phase3_plink/1000G.EUR.QC.${ch} --ld-wind-cm 1 --annot ', dir_out, '/Temp/S-LDSC/LOGODetect_theta_', theta[j], '_chr${ch}.annot --out ', dir_out, '/Temp/S-LDSC/LOGODetect_theta_', theta[j], '_chr${ch} --print-snps ./Data/hapmap3_snps/hm.${ch}.snp;
    done
    ./Data/ldsc/ldsc.py --rg ', dir_out, '/ldsc/dat1_reformated.sumstats.gz,', dir_out, '/ldsc/dat2_reformated.sumstats.gz --ref-ld-chr ', dir_out, '/Temp/S-LDSC/LOGODetect_theta_', theta[j], '_chr --w-ld-chr ./Data/ldsc/weights_hm3_no_hla/weights. --out ', dir_out, '/Temp/S-LDSC/LOGODetect_theta_', theta[j])
    system(job)
  }
}


gencov_prop = rep(0, length(theta))
for(j in 1:length(theta)){
  if(file.exists(paste0(dir_out, '/Temp/S-LDSC/LOGODetect_theta_', theta[j], '.log'))){
    Read = file(paste0(dir_out, '/Temp/S-LDSC/LOGODetect_theta_', theta[j], '.log'), 'r')
    line = readLines(Read, n = 69)
    line = readLines(Read, n = 1)
    line = gsub("   ", " ", line)
    line = gsub("  ", " ", line)
    if(length(line)>0){
      line = strsplit(line, ':')[[1]][2]
      gencov_prop[j] = abs(as.numeric(strsplit(line, ' ')[[1]][3]))
    }else{
      gencov_prop[j] = 0
    }
    close(Read)
  }
}
if(max(abs(gencov_prop))>0){
  j = which.max(abs(gencov_prop))
  dat = read.table(paste0(dir_out, '/Temp/Result_different_theta/Result_theta_', theta[j], '.txt'), header = T)
  write.table(dat, paste0(dir_out, '/LOGODetect_result.txt'), col.names = T, row.names = F, quote = F)
}else{
  dat = 'No region with significant local genetic correlation is identified.'
  write.table(dat, paste0(dir_out, '/LOGODetect_result.txt'), col.names = T, row.names = F, quote = F)
}




