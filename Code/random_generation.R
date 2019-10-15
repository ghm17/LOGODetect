library(data.table)
warnings('off')
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)
ch = as.numeric(args[1])
len = read.table(paste0('Data/LD_matrix/chr', ch, '/snplist/block_len.txt'))[, 1]
count = length(len)
N = 5000
N0 = 503

############## create directory
if(!dir.exists('Temp/random_ld')){
  dir.create('Temp/random_ld')
}
if(!dir.exists('Temp/random_ld2')){
  dir.create('Temp/random_ld2')
}
if(!dir.exists(paste0('Temp/random_ld/chr', ch))){
  dir.create(paste0('Temp/random_ld/chr', ch))
}
if(!dir.exists(paste0('Temp/random_ld2/chr', ch))){
  dir.create(paste0('Temp/random_ld2/chr', ch))
}
for(i in 1:(3*N)){
  if(!dir.exists(paste0('Temp/random_ld/chr', ch, '/', i))){
    dir.create(paste0('Temp/random_ld/chr', ch, '/', i))
  }
}
for(i in 1:(2*N)){
  if(!dir.exists(paste0('Temp/random_ld2/chr', ch, '/', i))){
    dir.create(paste0('Temp/random_ld2/chr', ch, '/', i))
  }
}



thre = 1000
rmv = function(n, mu, Sigma){
  SVD = svd(Sigma, nu = 0, nv = nrow(Sigma))
  sv = SVD$d
  p = sum(sv>0)
  D = sqrt(sv[1:p])
  V = SVD$v[,1:p]
  X = matrix(rnorm(n*p), nrow = p, ncol = n)
  out = t(V%*%(X*D))
  return(out)
}

LD = list() 
ld_diag = list()
ld_offdiag = list()
ld2_diag = list()
ld2_offdiag = list()
pb = txtProgressBar(0, count, style = 3)
LD[[1]] = as.matrix(fread(paste0('Data/LD_matrix/chr', ch, '/LD_', 1, '_', 2, '.ld')))
LD[[1]][is.na(LD[[1]])] = 0
ld_diag[[1]] = LD[[1]][1:len[1], 1:len[1]]
ld_offdiag[[1]] = LD[[1]][1:len[1], (len[1]+1):(len[1]+len[2])]
ld_diag[[2]] = LD[[1]][(len[1]+1):(len[1]+len[2]), (len[1]+1):(len[1]+len[2])]
ld2_diag[[1]] = (N0-1)/(N0-2)*( ld_diag[[1]]%*%ld_diag[[1]] + ld_offdiag[[1]]%*%t(ld_offdiag[[1]]) ) - sum(len[1:2])/(N0-2)*ld_diag[[1]]
ld2_offdiag[[1]] = (N0-1)/(N0-2)*( ld_diag[[1]]%*%ld_offdiag[[1]] + ld_offdiag[[1]]%*%ld_diag[[2]] ) - sum(len[1:2])/(N0-2)*ld_offdiag[[1]]

ss = rmv(3*N, rep(0, len[1]), ld_diag[[1]])
tt = rmv(2*N, rep(0, len[1]), ld2_diag[[1]])
ss = cbind(1:(3*N), ss)
tt = cbind(1:(2*N), tt)
wri_ss = function(dat){
  i = dat[1]
  dat = as.data.frame(dat[-1])
  fwrite(dat, file = paste0('Temp/random_ld/chr', ch, '/', i, '/ss_1.txt'), sep = ' ', row.names = F, col.names = F)
}
wri_tt = function(dat){
  i = dat[1]
  dat = as.data.frame(dat[-1])
  fwrite(dat, file = paste0('Temp/random_ld2/chr', ch, '/', i, '/tt_1.txt'), sep = ' ', row.names = F, col.names = F)
}
apply(ss, 1, wri_ss)
apply(tt, 1, wri_tt)
ss = ss[, -1]
tt = tt[, -1]
setTxtProgressBar(pb, 1)

for(j in 2:(count-1)){
  LD[[j]] = as.matrix(fread(paste0('Data/LD_matrix/chr', ch, '/LD_', j, '_', j + 1, '.ld')))
  LD[[j]][is.na(LD[[j]])] = 0
  ld_diag[[j+1]] = LD[[j]][(len[j]+1):(len[j]+len[j+1]), (len[j]+1):(len[j]+len[j+1])]
  ld_offdiag[[j]] = matrix(LD[[j]][1:len[j], (len[j]+1):(len[j]+len[j+1])], nrow = len[j], ncol = len[j+1])
  ld2_diag[[j]] = (N0-1)/(N0-2)*( ld_diag[[j]]%*%ld_diag[[j]] + t(ld_offdiag[[j-1]])%*%ld_offdiag[[j-1]] + ld_offdiag[[j]]%*%t(ld_offdiag[[j]]) ) - sum(len[(j-1):(j+1)])/(N0-2)*ld_diag[[j]]
  ld2_offdiag[[j]] = matrix( (N0-1)/(N0-2)*( ld_diag[[j]]%*%ld_offdiag[[j]] + ld_offdiag[[j]]%*%ld_diag[[j+1]] ) - sum(len[j:(j+1)])/(N0-2)*ld_offdiag[[j]], nrow = len[j], ncol = len[j+1])
  
  svd1 = svd(ld_diag[[j-1]], nu = 0, nv = nrow(ld_diag[[j-1]]))
  sv1 = svd1$d
  p1 = sum(sv1>sv1[1]/thre)
  D1 = diag(1/sv1[1:p1])
  V1 = svd1$v[,1:p1]
  S11 = ld_diag[[j-1]]
  S12 = ld_offdiag[[j-1]]
  S22 = ld_diag[[j]]
  C1 = t(S12)%*%V1%*%D1%*%t(V1)
  Sigma1 = S22 - C1%*%S12
  
  svd2 = svd(ld2_diag[[j-1]], nu = 0, nv = nrow(ld2_diag[[j-1]]))
  sv2 = svd2$d
  p2 = sum(sv2>sv2[1]/thre)
  D2 = diag(1/sv2[1:p2])
  V2 = svd2$v[,1:p2]
  T11 = ld2_diag[[j-1]]
  T12 = ld2_offdiag[[j-1]]
  T22 = ld2_diag[[j]]
  C2 = t(T12)%*%V2%*%D2%*%t(V2)
  Sigma2 = T22 - C2%*%T12
  
  ss_new = rmv(3*N, rep(0, len[j]), Sigma1) + ss%*%t(C1)
  tt_new = rmv(2*N, rep(0, len[j]), Sigma2) + tt%*%t(C2)
  ss = ss_new
  tt = tt_new
  ss = cbind(1:(3*N), ss)
  tt = cbind(1:(2*N), tt)
  wri_ss = function(dat){
    i = dat[1]
    dat = as.data.frame(dat[-1])
    fwrite(dat, file = paste0('Temp/random_ld/chr', ch, '/', i, '/ss_', j, '.txt'), sep = ' ', row.names = F, col.names = F)
  }
  wri_tt = function(dat){
    i = dat[1]
    dat = as.data.frame(dat[-1])
    fwrite(dat, file = paste0('Temp/random_ld2/chr', ch, '/', i, '/tt_', j, '.txt'), sep = ' ', row.names = F, col.names = F)
  }
  apply(ss, 1, wri_ss)
  apply(tt, 1, wri_tt)
  ss = ss[, -1]
  tt = tt[, -1]
  
  LD[[j-1]] = 0
  ld_diag[[j-1]] = 0
  ld_offdiag[[j-1]] = 0
  ld2_diag[[j-1]] = 0
  ld2_offdiag[[j-1]] = 0
  setTxtProgressBar(pb, j)
}

j = count
ld2_diag[[j]] = ld_diag[[j]]%*%ld_diag[[j]] + t(ld_offdiag[[j-1]])%*%ld_offdiag[[j-1]]
svd1 = svd(ld_diag[[j-1]], nu = 0, nv = nrow(ld_diag[[j-1]]))
sv1 = svd1$d
p1 = sum(sv1>sv1[1]/thre)
D1 = diag(1/sv1[1:p1])
V1 = svd1$v[,1:p1]
S11 = ld_diag[[j-1]]
S12 = ld_offdiag[[j-1]]
S22 = ld_diag[[j]]
C1 = t(S12)%*%V1%*%D1%*%t(V1)
Sigma1 = S22 - C1%*%S12

svd2 = svd(ld2_diag[[j-1]], nu = 0, nv = nrow(ld2_diag[[j-1]]))
sv2 = svd2$d
p2 = sum(sv2>sv2[1]/thre)
D2 = diag(1/sv2[1:p2])
V2 = svd2$v[,1:p2]
T11 = ld2_diag[[j-1]]
T12 = ld2_offdiag[[j-1]]
T22 = ld2_diag[[j]]
C2 = t(T12)%*%V2%*%D2%*%t(V2)
Sigma2 = T22 - C2%*%T12
ss_new = rmv(3*N, rep(0, len[j]), Sigma1) + ss%*%t(C1)
tt_new = rmv(2*N, rep(0, len[j]), Sigma2) + tt%*%t(C2)
ss = ss_new
tt = tt_new
ss = cbind(1:(3*N), ss)
tt = cbind(1:(2*N), tt)
wri_ss = function(dat){
  i = dat[1]
  dat = as.data.frame(dat[-1])
  fwrite(dat, file = paste0('Temp/random_ld/chr', ch, '/', i, '/ss_', j, '.txt'), sep = ' ', row.names = F, col.names = F)
}
wri_tt = function(dat){
  i = dat[1]
  dat = as.data.frame(dat[-1])
  fwrite(dat, file = paste0('Temp/random_ld2/chr', ch, '/', i, '/tt_', j, '.txt'), sep = ' ', row.names = F, col.names = F)
}
apply(ss, 1, wri_ss)
apply(tt, 1, wri_tt)
ss = ss[, -1]
tt = tt[, -1]
setTxtProgressBar(pb, count)
rm(LD, ld_diag, ld_offdiag, ld2_diag, ld2_offdiag)

