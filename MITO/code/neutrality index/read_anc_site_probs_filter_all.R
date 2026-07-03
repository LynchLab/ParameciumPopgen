## remove sites with low posterior prob of ancestral state in any species

args <- commandArgs(trailingOnly = T)
fin <- args[1]
fout <- args[2]

## threshold
a <- 0.85

df <- read.csv( file = fin, header = T, check.names = F, row.names = 1)

probs <- t( apply( df, 1, function(x) as.numeric( sapply( strsplit(x, split = ' ', fixed = T), '[', 2))))

colnames(probs) <- sapply( strsplit( x = colnames(df), split = '-', fixed = T), '[', 1)
colnames(probs) <- sapply( strsplit( x = colnames(probs), split = '_', fixed = T), '[', 2)

probs <- as.data.frame(probs)

## list sites with probability less than threshold
fs <- apply( probs, 2, function(x) unname(which( x < a)))

## sites filtered in any species
comfs <- sort( Reduce( union, fs))

out <- lapply( fs, function(X) comfs)

out <- as.data.frame( sapply( out, paste, collapse = '-'))
colnames(out) <- NA

write.table( x = out, file = fout, quote = F, row.names = T, col.names = F, sep = " ")