
library(data.table)
warnings('off')
options(stringsAsFactors=F)
options(datatable.fread.datatable=F)
args = commandArgs(trailingOnly=TRUE)
ch = as.numeric(args[1])
n1 = as.numeric(args[2])  ##sample size for trait1
n2 = as.numeric(args[3])  ##sample size for trait2
dir_out = as.character(args[4])


################# Read ldsc results including heritability for two traits and intercept(explaining sample overlaps)
Read = file(paste0(dir_out, '/ldsc/ldsc_trait1_trait2.log'), 'r')
line = readLines(Read, n = 34)
line = readLines(Read, n = 1)
line = gsub("   ", " ", line)
line = gsub("  ", " ", line)
line = strsplit(line, ':')[[1]][2]
h2_1 = as.numeric(strsplit(line, " ")[[1]][2:23])
line = readLines(Read, n = 15)
line = readLines(Read, n = 1)
line = gsub("   ", " ", line)
line = gsub("  ", " ", line)
line = strsplit(line, ':')[[1]][2]
h2_2 = as.numeric(strsplit(line, " ")[[1]][2:23])
line = readLines(Read, n = 21)
line = readLines(Read, n = 1)
intercept = as.numeric(strsplit(line, " ")[[1]][2])
close(Read)
var.beta = max(h2_1[ch], 0)
var.gamma = max(h2_2[ch], 0)


##################################### BiScan null distribution
LDSC = fread(paste0('Data/ldsc/l2/chr', ch, '.l2.ldscore'), header = T)
ref_snp = LDSC$SNP
dat = fread(paste0(dir_out, '/Data_QC/dat1_chr', ch, '.txt'), header = T)
dat_snp = dat$SNP ##snplist in the data of interest
ldsc = LDSC$L2
len = fread(paste0('Data/LD_matrix/chr', ch, '/snplist/block_len.txt'))[, 1]  ##size for each block
count = length(len) ##blocks count
N = 5000   ##size of scan stat empirical distribution
M = length(dat_snp)
ind = ref_snp%in%dat_snp
ldsc = ldsc[ind]
for(i in 1:M) ldsc[i] = max(ldsc[i], 1)

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

## change parameter
Cn = 2000
inter = 20
theta = c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7)

dat_merge = fread(paste0('Data/LD_block/ldblock_merged_chr', ch, '.txt'), header = T)
frag_count = nrow(dat_merge)
Qmax = list()
for(j in 1:length(theta)){
  Qmax[[j]] = matrix(0, nrow = N, ncol = frag_count)
}
sd1 = rep(0, N)
sd2 = rep(0, N)
frag = matrix(0, nrow = frag_count, ncol = 2)
ind_noninf = rep(0, frag_count)
for(k in 1:frag_count){
  frag[k, 1] = min(which( dat$pos >= dat_merge[k, 1] ))
  frag[k, 2] = max(which( dat$pos <= dat_merge[k, 2] ))
  ind_noninf[k] = is.finite(frag[k, 1]) & is.finite(frag[k, 2])
}

pb = txtProgressBar(0, N, style = 3)
for(i in 1:N){
  s1 = fread(paste0('Data/random_ld/chr', ch, '/s_', 2*i-1, '.txt'))[, 1]
  s2 = fread(paste0('Data/random_ld/chr', ch, '/s_', 2*i, '.txt'))[, 1]
  t1 = fread(paste0('Data/random_ld2/chr', ch, '/t_', 2*i-1, '.txt'))[, 1]
  t2 = fread(paste0('Data/random_ld2/chr', ch, '/t_', 2*i, '.txt'))[, 1]
  r = fread(paste0('Data/random_ld/chr', ch, '/s_', 2*N + i, '.txt'))[, 1]
  
  u1 = sqrt(n1/M*var.beta)*t1 + sqrt(1 - var.beta - abs(intercept))*s1
  u2 = sqrt(n2/M*var.gamma)*t2 + sqrt(1 - var.gamma - abs(intercept))*s2
  v1 = r * sqrt(abs(intercept))
  v2 = r * sqrt(abs(intercept)) * sign(intercept)
  z1 = u1 + v1
  z2 = u2 + v2
  z1 = z1[ind]
  z2 = z2[ind]
  sd1[i] = sqrt(mean(z1^2))
  sd2[i] = sqrt(mean(z2^2))
  for(j in 1:length(theta)){
    for(k in c(1:frag_count)[as.logical(ind_noninf)]){
      re = scan(z1[frag[k,1]:frag[k,2]], z2[frag[k,1]:frag[k,2]], ldsc[frag[k,1]:frag[k,2]], Cn, inter, frag[k,1], frag[k,2], theta[j])
      Qmax[[j]][i, k] = as.numeric(re[[1]])
    }
  }
  setTxtProgressBar(pb, i)
}

if(!dir.exists(paste0(dir_out, '/Temp'))){
  dir.create(paste0(dir_out, '/Temp'))
}
if(!dir.exists(paste0(dir_out, '/Temp/sd'))){
  dir.create(paste0(dir_out, '/Temp/sd'))
}
if(!dir.exists(paste0(dir_out, '/Temp/Qmax'))){
  dir.create(paste0(dir_out, '/Temp/Qmax'))
}
write.table(sd1, paste0(dir_out, '/Temp/sd/sd1_chr', ch, '.txt'), quote = F, col.names = F, row.names = F)
write.table(sd2, paste0(dir_out, '/Temp/sd/sd2_chr', ch, '.txt'), quote = F, col.names = F, row.names = F)
for(j in 1:length(theta)){
  write.table(Qmax[[j]], paste0(dir_out, '/Temp/Qmax/Qmax_chr', ch, '_theta_', theta[j], '.txt'), quote = F, col.names = F, row.names = F)
}





