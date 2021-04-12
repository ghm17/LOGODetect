library(data.table)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)
ch = as.numeric(args[1])
len = read.table(paste0('Data/LD_matrix/chr', ch, '/snplist/block_len.txt'))[, 1]  ##size for each block
count = length(len)
N = 5000
ref_snp = read.table(paste0('Data/ldsc/l2/chr', ch, '.l2.ldscore'), header = T)$SNP

############## create directory
if(!dir.exists('Data/random_ld')){
  dir.create('Data/random_ld')
}
if(!dir.exists('Data/random_ld2')){
  dir.create('Data/random_ld2')
}
if(!dir.exists(paste0('Data/random_ld/chr', ch))){
  dir.create(paste0('Data/random_ld/chr', ch))
}
if(!dir.exists(paste0('Data/random_ld2/chr', ch))){
  dir.create(paste0('Data/random_ld2/chr', ch))
}





pb = txtProgressBar(0, (3*N), style = 3)
for(i in 1:(3*N)){
  s = rep(0, length(ref_snp))
  for(j in 1:count){
    ss = fread(paste0('Temp/random_ld/chr', ch, '/', i, '/ss_', j, '.txt'))[[1]]
    s[( sum(len[1:j]) - len[j] + 1 ) : sum(len[1:j]) ] = ss
  }
  unlink(paste0('Temp/random_ld/chr', ch, '/', i), recursive = T)
  s = data.frame(s)
  fwrite(s, file = paste0('Data/random_ld/chr', ch, '/s_', i, '.txt'), sep = ' ', row.names = F, col.names = F)
  setTxtProgressBar(pb, i)
}

pb = txtProgressBar(0, (2*N), style = 3)
for(i in 1:(2*N)){
  t = rep(0, length(ref_snp))
  for(j in 1:count){
    tt = fread(paste0('Temp/random_ld2/chr', ch, '/', i, '/tt_', j, '.txt'))[[1]]
    t[( sum(len[1:j]) - len[j] + 1 ) : sum(len[1:j]) ] = tt
  }
  unlink(paste0('Temp/random_ld2/chr', ch, '/', i), recursive = T)
  t = data.frame(t)
  fwrite(t, file = paste0('Data/random_ld2/chr', ch, '/t_', i, '.txt'), sep = ' ', row.names = F, col.names = F)
  setTxtProgressBar(pb, i)
}
