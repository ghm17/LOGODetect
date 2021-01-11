library(data.table)
library(MASS)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)
dir_out = as.character(args[1])

scan = function(z1, z2, ldsc, Cn, inter, le, ri, theta){
  m = length(z1)
  z1_time_z2 = rep(0, m + 1)
  ldscore_sum = rep(0, m + 1)
  for(i in 1:m){
    z1_time_z2[i+1] = z1_time_z2[i]+ z1[i]*z2[i]
    ldscore_sum[i+1] = ldscore_sum[i] + ldsc[i]
  }
  if(m > Cn){
    Qmax = -Inf
    window = seq(from = inter, to = floor(Cn/inter)*inter, by = inter)
    for(I in window){
      j_count = ceiling((m-I+1)/inter)
      qq = rep(0, j_count)
      for(j in 1:j_count){
        qq[j] = (z1_time_z2[(j-1)*inter+1 +I] - z1_time_z2[(j-1)*inter+1])/(ldscore_sum[(j-1)*inter+1 +I] - ldscore_sum[(j-1)*inter+1])^theta
        qq[j] = abs(qq[j])
      }
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
  Qmax = -Inf
  for(I in 1:m){
    for(j in 1:(m-I+1)){
      Q = (z1_time_z2[j+I] - z1_time_z2[j])/(ldscore_sum[j+I] - ldscore_sum[j])^theta
      Q = abs(Q)
      if(Qmax < Q){
        Qmax = Q
        low = j
        high = j + I - 1
      }
    }
  }
  re = list()
  re[[1]] = Qmax
  re[[2]] = c(low + le - 1, high + le - 1)
  return(re)
}
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

Cn = 2000
inter = 20
theta = c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7)
N = 5000   ##size of scan stat empirical distribution

Result = list()
for(j in 1:length(theta)){
  Result[[j]] = as.data.frame(matrix(0, nrow = 0, ncol = 9))
  colnames(Result[[j]]) = c('chr', 'begin', 'stop', 'count', 'stat', 'pval', 'begin_pos', 'stop_pos', 'screening_frag')
}

dat_merge = list()
Frag_count = rep(0, 22)
for(ch in 1:22){
  dat_merge[[ch]] = read.table(paste0('Data/LD_block/ldblock_merged_chr', ch, '.txt'), header = T)
  Frag_count[ch] = nrow(dat_merge[[ch]])
}
testing_count = sum(Frag_count)
dat1 = list()
dat2 = list()
for(ch in 1:22){
  frag_count = Frag_count[ch]
  LDSC = read.table(paste0('Data/ldsc/l2/chr', ch, '.l2.ldscore'), header = T)
  ref_snp = LDSC$SNP
  dat_snp = read.table(paste0(dir_out, '/Data_QC/dat1_chr', ch, '.txt'), header = T)$SNP
  ldsc = LDSC$L2
  ind = ref_snp%in%dat_snp
  ldsc = ldsc[ind]
  M = length(dat_snp)
  for(i in 1:M) ldsc[i] = max(ldsc[i], 1)
  dat1[[ch]] = read.table(paste0(dir_out, '/Data_QC/dat1_chr', ch, '.txt'), header = T)
  dat2[[ch]] = read.table(paste0(dir_out, '/Data_QC/dat2_chr', ch, '.txt'), header = T)
  z1 = dat1[[ch]]$Z
  z2 = dat2[[ch]]$Z
  Qmax = list()
  for(j in 1:length(theta)){
    Qmax[[j]] = as.matrix(read.table(paste0(dir_out, '/Temp/Qmax/Qmax_chr', ch, '_theta_', theta[j], '.txt')))
  }
  sd1_data = sqrt(mean(z1^2))
  sd2_data = sqrt(mean(z2^2))
  sd1_generated = read.table(paste0(dir_out, '/Temp/sd/sd1_chr', ch, '.txt'))[, 1]
  sd2_generated = read.table(paste0(dir_out, '/Temp/sd/sd2_chr', ch, '.txt'))[, 1]
  
  q = list()
  Qmax_scaled = list()
  for(j in 1:length(theta)){
    Qmax_scaled[[j]] = matrix(0, nrow = N, ncol = frag_count)
    q[[j]] = rep(0, frag_count)
    for(k in 1:frag_count){
      Qmax_scaled[[j]][, k] = Qmax[[j]][, k]*(sd1_data * sd2_data)/(sd1_generated * sd2_generated)
      q[[j]][k] = quantile(Qmax_scaled[[j]][,k], 0.95)
    }
  }
  
  frag = matrix(0, nrow = frag_count, ncol = 2)
  ind_noninf = rep(0, frag_count)
  for(k in 1:frag_count){
    frag[k, 1] = min(which( dat1[[ch]]$pos >= dat_merge[[ch]][k, 1] ))
    frag[k, 2] = max(which( dat1[[ch]]$pos <= dat_merge[[ch]][k, 2] ))
    ind_noninf[k] = is.finite(frag[k, 1]) & is.finite(frag[k, 2])
  }
  for(j in 1:length(theta)){
    for(k in c(1:frag_count)[as.logical(ind_noninf)]){
      result = whole_scan(z1[frag[k,1]:frag[k,2]], z2[frag[k,1]:frag[k,2]], ldsc[frag[k,1]:frag[k,2]], Cn, inter, q[[j]][k], frag[k,1], frag[k,2], theta[j])
      detect = c()
      if(length(result)>0){
        for(i in 1:(length(result)/2)){
          begin = result[[2*i]][1]
          stop = result[[2*i]][2]
          count = stop - begin + 1
          stat = sum(z1[begin:stop]*z2[begin:stop])/(sum(ldsc[begin:stop]))^theta[j]
          pval = (sum(abs(stat) < Qmax_scaled[[j]][,k]) + 1)/(N + 1)
          begin_pos = round(dat1[[ch]]$pos[begin]/1000000, 3)
          stop_pos = round(dat1[[ch]]$pos[stop]/1000000, 3)
          detect_new = c(begin, stop, count, stat, pval, begin_pos, stop_pos)
          detect = rbind(detect, detect_new)
        }
        detect = data.frame(rep(ch, length(result)/2), detect, rep(k, length(result)/2), row.names = NULL)
        colnames(detect) = c('chr', 'begin', 'stop', 'count', 'stat', 'pval', 'begin_pos', 'stop_pos', 'screening_frag')
        Result[[j]] = rbind(Result[[j]], detect)
      }
    }
  }
}

if(!dir.exists(paste0(dir_out, '/Temp/Result_different_theta'))){
  dir.create(paste0(dir_out, '/Temp/Result_different_theta'))
}

for(j in 1:length(theta)){
  Result[[j]] = Result[[j]][order(Result[[j]]$stat, decreasing = T),]
  Result[[j]] = Result[[j]][order(Result[[j]]$pval),]
  Result[[j]]$testing_count = testing_count
  n0 = nrow(Result[[j]])
  if(n0>1){
    for(i in 2:n0){
      temp = Result[[j]][1:(i-1), ]
      testing_add = (sum(temp$chr == Result[[j]]$chr[i] & temp$screening_frag == Result[[j]]$screening_frag[i]))>0
      Result[[j]]$testing_count[i:n0] = Result[[j]]$testing_count[i:n0] + testing_add
    }
  }
  
  p_val = Result[[j]]$pval
  BH = max(which(p_val<=0.05/Result[[j]]$testing_count*c(1:n0)))
  Result[[j]]$qval = Result[[j]]$pval * Result[[j]]$testing_count / c(1:n0)
  if(BH>0){
    write.table(Result[[j]][1:BH, c('chr', 'begin_pos', 'stop_pos', 'stat', 'pval', 'qval')], paste0(dir_out, '/Temp/Result_different_theta/Result_theta_', theta[j], '.txt'), col.names = T, row.names = F, quote = F)
  }
}

