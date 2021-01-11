library(data.table)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)
dir1 = as.character(args[1])
dir2 = as.character(args[2])
dir_out = as.character(args[3])

dat1 = as.data.frame(fread(dir1, header = T))
dat2 = as.data.frame(fread(dir2, header = T))

dat1_chr = list()
dat1_new = NULL
for(i in 1:22){
  dat1_chr[[i]] = dat1[dat1$chr == i, ]
  dat1_chr[[i]] = dat1_chr[[i]][order(dat1_chr[[i]]$pos), ]
  dat1_new = rbind(dat1_new, dat1_chr[[i]])
}
rm(dat1_chr)
dat1 = dat1_new
dat2_chr = list()
dat2_new = NULL
for(i in 1:22){
  dat2_chr[[i]] = dat2[dat2$chr == i, ]
  dat2_chr[[i]] = dat2_chr[[i]][order(dat2_chr[[i]]$pos), ]
  dat2_new = rbind(dat2_new, dat2_chr[[i]])
}
rm(dat2_chr)
dat2 = dat2_new


## quality control
dat1 = dat1[complete.cases(dat1),]
dat2 = dat2[complete.cases(dat2),]
dat1 = dat1[! ((dat1$A1=="G")&(dat1$A2=="C") | (dat1$A1=="C")&(dat1$A2=="G") | (dat1$A1=="A")&(dat1$A2=="T") | (dat1$A1=="T")&(dat1$A2=="A")),]
dat2 = dat2[! ((dat2$A1=="G")&(dat2$A2=="C") | (dat2$A1=="C")&(dat2$A2=="G") | (dat2$A1=="A")&(dat2$A2=="T") | (dat2$A1=="T")&(dat2$A2=="A")),]

# remove duplicated SNPs
dat1 = dat1[! dat1$SNP %in% dat1$SNP[duplicated(dat1$SNP)],]
dat2 = dat2[! dat2$SNP %in% dat2$SNP[duplicated(dat2$SNP)],]

## take overlap SNPs
dat1_new = dat1[1,][-1,]
dat2_new = dat2[1,][-1,]
for(i in 1:22){
  f1 = dat1[dat1$chr == i, ]
  f2 = dat2[dat2$chr == i, ]
  f1 = f1[! f1$pos %in% f1$pos[duplicated(f1$pos)],]
  f2 = f2[! f2$pos %in% f2$pos[duplicated(f2$pos)],]
  overlap_chr = intersect(f1$pos, f2$pos)
  f1 = f1[f1$pos %in% overlap_chr, ]
  f2 = f2[f2$pos %in% overlap_chr, ]
  dat1_new = rbind(dat1_new, f1)
  dat2_new = rbind(dat2_new, f2)
}
dat1 = dat1_new
dat2 = dat2_new
rm(dat1_new, dat2_new)


#quality control
dat1$A1 = as.character(dat1$A1)
dat1$A2 = as.character(dat1$A2)
dat2$A1 = as.character(dat2$A1)
dat2$A2 = as.character(dat2$A2)
pos = (dat1$A1 == dat2$A1)*(dat1$A2 == dat2$A2)
neg = (dat1$A1 == dat2$A2)*(dat1$A2 == dat2$A1)
flop = (neg == 1)
dat1$A1[flop] = dat2$A1[flop]
dat1$A2[flop] = dat2$A2[flop]
dat1$Z[flop] = -dat1$Z[flop]
ind = (pos + neg == 1)
dat1 = dat1[ind, ]
dat2 = dat2[ind, ]

dat1_final = NULL
dat2_final = NULL
dat1_chr = list()
dat2_chr = list()

if(!dir.exists(paste0(dir_out, '/Data_QC'))){
  dir.create(paste0(dir_out, '/Data_QC'))
}
for(ch in 1:22){
  dat1_chr[[ch]] = dat1[dat1$chr==ch, ]
  dat2_chr[[ch]] = dat2[dat2$chr==ch, ]
  ref_snp = read.table(paste0('Data/ldsc/l2/chr', ch, '.l2.ldscore'), header = T)$SNP
  snp = dat1_chr[[ch]]$SNP
  overlap_snp = snp[snp%in%ref_snp]
  dat1_chr[[ch]] = dat1_chr[[ch]][snp%in%ref_snp, ]
  dat2_chr[[ch]] = dat2_chr[[ch]][snp%in%ref_snp, ]
  write.table(dat1_chr[[ch]], paste0(dir_out, '/Data_QC/dat1_chr', ch, '.txt'), row.names = F, col.names = T, quote = F)
  write.table(dat2_chr[[ch]], paste0(dir_out, '/Data_QC/dat2_chr', ch, '.txt'), row.names = F, col.names = T, quote = F)
  dat1_final = rbind(dat1_final, dat1_chr[[ch]])
  dat2_final = rbind(dat2_final, dat2_chr[[ch]])
}
write.table(dat1_final, paste0(dir_out, '/Data_QC/dat1.txt'), row.names = F, col.names = T, quote = F)
write.table(dat2_final, paste0(dir_out, '/Data_QC/dat2.txt'), row.names = F, col.names = T, quote = F)

