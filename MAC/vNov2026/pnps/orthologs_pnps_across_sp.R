fog <- 'orthologs_transcripts.tsv' # transcript id table of orthologs
fsm <- 'sp_ref_map.txt' # species reference map
indr <- 'data' # dir of species tables of genes' pnps
fout <- 'orthologs_species_pnps.csv'

dog <- read.csv(file = fog, sep = '\t', row.names = 1)
dsm <- read.csv(file = fsm, header = F, sep = ' ')
nsp <- nrow(dsm)

out <- matrix( NA, nrow(dog), nsp, dimnames = list( rownames(dog), dsm[,1]))

for(i in 1:nsp){
  fin <- file.path( indr, paste0(dsm[i,1],".csv"))
  df <- read.csv(file = fin, row.names = 1, check.names = F)
  out[,i] <- df[ dog[ , dsm[i,2]], 'pN/pS']
}

write.csv( x = out, file = fout, quote = F)