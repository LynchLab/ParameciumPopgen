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

out <- as.data.frame( sapply( fs, paste, collapse = '-'))
colnames(out) <- NA

write.table( x = out, file = fout, quote = F, row.names = T, col.names = F, sep = " ")