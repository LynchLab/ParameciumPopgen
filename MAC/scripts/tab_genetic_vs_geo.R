fgeo <- 'data/table_S1.csv'
indr1 <- 'data/geodist'
indr2 <- 'data/mean_pdist_d4'
outdr <- 'data/genetic_vs_geo'

fls1 <- list.files(indr1)
fls2 <- list.files(indr2)
fls <- intersect( fls1, fls2)
sps <- gsub('.csv', '', fls)

nsp <- length(sps)

for( i in 1:nsp){
  f1 <- file.path( indr1, paste0(sps[i],'.csv'))
  f2 <- file.path( indr2, paste0(sps[i],'.csv'))
  fout <- file.path( outdr, paste0(sps[i],'.csv'))
  
  dgeo <- read.csv( file = f1, check.names = F, row.names = 1)
  dgen <- read.csv( file = f2, check.names = F, row.names = 1)
  dmat <- as.matrix( as.dist(t(dgen)))
  samples <- intersect( rownames(dgen), rownames(dgeo))
  nsp <- length(samples)
  out <- NULL
  for( j in 1:(nsp-1)){
    s1 <- samples[j]
    for( k in (j+1):nsp){
      s2 <- samples[k]
      entry <- c( s1, s2, dgeo[s1,s2], dmat[s1,s2])
      out <- rbind( out, entry)
    }
  }
  colnames(out) <- c('Sample1','Sample2','Geographic_Distance','Genetic_Distance')
  write.csv( x = out, file = fout, quote = F, row.names = F)
}