library(data.table)
library(MASS)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

scan = function(z1, z2, ldsc, Cn, inter, le, ri, alpha){
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
        qq[j] = (z1_time_z2[(j-1)*inter+1 +I] - z1_time_z2[(j-1)*inter+1])/(ldscore_sum[(j-1)*inter+1 +I] - ldscore_sum[(j-1)*inter+1])^alpha
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
      Q = (z1_time_z2[j+I] - z1_time_z2[j])/(ldscore_sum[j+I] - ldscore_sum[j])^alpha
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
whole_scan = function(z1, z2, ldsc, Cn, inter, thre, le, ri, alpha){
  re = scan(z1, z2, ldsc, Cn, inter, le, ri, alpha)
  if(re[[1]] <= thre){
    return(list())
  }
  else{
    a = re[[2]][1]
    b = re[[2]][2]
    if(a == le & b == ri){return(re)}
    if(a == le & b < ri){
      return( c( re , whole_scan(z1[(b+2-le):(ri+1-le)], z2[(b+2-le):(ri+1-le)], ldsc[(b+2-le):(ri+1-le)], Cn, inter, thre, b+1, ri, alpha) ) )
    }
    if(a > le & b == ri){
      return( c( whole_scan(z1[1:(a-le)], z2[1:(a-le)], ldsc[1:(a-le)], Cn, inter, thre, le, a-1, alpha), re) )
    }
    return( c( whole_scan(z1[1:(a-le)], z2[1:(a-le)], ldsc[1:(a-le)], Cn, inter, thre, le, a-1, alpha), re, whole_scan(z1[(b+2-le):(ri+1-le)], z2[(b+2-le):(ri+1-le)], ldsc[(b+2-le):(ri+1-le)], Cn, inter, thre, b+1, ri, alpha) ) )
  }
}

Cn = 2000
inter = 20
alpha = 1/2
N = 5000   ##size of scan stat empirical distribution

Result = as.data.frame(matrix(0, nrow = 1, ncol = 8))
colnames(Result) = c('chr', 'begin', 'stop', 'count', 'stat', 'pval', 'begin_pos', 'stop_pos')
Result = Result[-1, ]
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
  dat_snp = read.table(paste0('Temp/Data_QC/dat1_chr', ch, '.txt'))[, 1]
  ldsc = LDSC$L2
  ind = ref_snp%in%dat_snp
  ldsc = ldsc[ind]
  M = length(dat_snp)
  for(i in 1:M) ldsc[i] = max(ldsc[i], 1)
  dat1[[ch]] = read.table(paste0('Temp/Data_QC/dat1_chr', ch, '.txt'), header = T)
  dat2[[ch]] = read.table(paste0('Temp/Data_QC/dat2_chr', ch, '.txt'), header = T)
  z1 = dat1[[ch]]$Z
  z2 = dat2[[ch]]$Z
  Qmax = as.matrix(read.table(paste0('Temp/Qmax/Qmax_chr', ch, '.txt')))
  sd1_data = sqrt(mean(z1^2))
  sd2_data = sqrt(mean(z2^2))
  sd1_generated = read.table(paste0('Temp/sd/sd1_chr', ch, '.txt'))[, 1]
  sd2_generated = read.table(paste0('Temp/sd/sd2_chr', ch, '.txt'))[, 1]
  
  q = rep(0, frag_count)
  for(k in 1:frag_count){
    Qmax[,k] = Qmax[,k]*(sd1_data * sd2_data)/(sd1_generated * sd2_generated)
    q[k] = quantile(Qmax[,k], 0.95)
  }
  
  frag = matrix(0, nrow = frag_count, ncol = 2)
  ind_noninf = rep(0, frag_count)
  for(k in 1:frag_count){
    frag[k, 1] = min(which( dat1[[ch]]$pos >= dat_merge[[ch]][k, 1] ))
    frag[k, 2] = max(which( dat1[[ch]]$pos <= dat_merge[[ch]][k, 2] ))
    ind_noninf[k] = is.finite(frag[k, 1]) & is.finite(frag[k, 2])
  }
  for(k in c(1:frag_count)[as.logical(ind_noninf)]){
    result = whole_scan(z1[frag[k,1]:frag[k,2]], z2[frag[k,1]:frag[k,2]], ldsc[frag[k,1]:frag[k,2]], Cn, inter, q[k], frag[k,1], frag[k,2], alpha)
    detect = c()
    if(length(result)>0){
      for(i in 1:(length(result)/2)){
        begin = result[[2*i]][1]
        stop = result[[2*i]][2]
        count = stop - begin + 1
        stat = sum(z1[begin:stop]*z2[begin:stop])/(sum(ldsc[begin:stop]))^alpha
        pval = (sum(abs(stat) < Qmax[,k]) + 1)/(N + 1)
        begin_pos = round(dat1[[ch]]$pos[begin]/1000000, 3)
        stop_pos = round(dat1[[ch]]$pos[stop]/1000000, 3)
        detect_new = c(begin, stop, count, stat, pval, begin_pos, stop_pos)
        detect = rbind(detect, detect_new)
      }
      detect = data.frame(rep(ch, length(result)/2), detect, row.names = NULL)
      colnames(detect) = c('chr', 'begin', 'stop', 'count', 'stat', 'pval', 'begin_pos', 'stop_pos')
      Result = rbind(Result, detect)
    }
  }
}


Result = Result[order(Result$stat, decreasing = T),]
Result = Result[order(Result$pval),]
p_val = Result$pval
BH = max(which(p_val<=0.05/testing_count*c(1:nrow(Result))))
Result$testing_count = testing_count
Result$qval = Result$pval * Result$testing_count / c(1:nrow(Result))
if(BH>0){
  write.table(Result[1:BH, c(1, 5, 6, 7, 8, 10)], 'Result.txt', row.names = F, col.names = T, quote = F)
}
if(BH<=0){
  write.table('detect nothing', 'Result.txt', row.names = F, col.names = T, quote = F)
}

