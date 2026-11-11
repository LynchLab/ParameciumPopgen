#### Are pNpS of genes correlated across species? ####
fin <- 'orthologs_species_pnps.csv'
fout <- 'figure_S8.png'

df <- read.csv(file = fin, row.names = 1)

sps <- colnames(df)
nsp <- length(sps)

spshort <- sapply( strsplit(x = sps, split = '', fixed = T), function(x) paste(x[1:4],collapse=''))

png(fout, width = 5, height = 5, units = 'in', pointsize = 8, res = 150)
par( mfrow = c(nsp-1,nsp-1), mar = c(0,0,0,0), oma = c(5,4,2,2))

k = 1

cobs <- NULL

# for each unique pair of species
for(i in 1:(nsp-1)){
  for(j in 2:nsp){
    if( j <= i){
      # do not draw unless j > 1
      plot(NA, xlim = c(-3,1), ylim = c(-3,1), xlab = '', ylab = '', axes = F)
    }
    else{
      # take the subset from the full dataframe
      sub <- df[ , c(i,j)]
      # take genes with data for both species, excluding genes with undefined pNpS
      sub <- sub[ !apply( sub, 1, function(x) any(is.infinite(x) | is.na(x))), ]
      lsub <- sub[ !apply( sub, 1, function(x) any(x==0)), ]
      lsub$log1 <- log10(lsub[,1])
      lsub$log2 <- log10(lsub[,2])
      plot( log1 ~ log2, data = lsub, xlim = c(-4,2), ylim = c(-4,2), xlab = '', ylab = '', axes = F)
      box(col='grey')
      # linear regression
      lmob <- lm( log2 ~ log1, data = lsub)
      abline(lmob,lty=2)
      # spearman correlation
      scorob <- cor.test(x = sub[,1], y = sub[,2], method = 's')
      r <- round( scorob$estimate, 2)
      p <- ifelse( scorob$p.value < 0.05, '*', '')
      cobs <- rbind( cobs, c( sps[ c(i,j)], scorob$estimate, scorob$p.value))
      mtext( text = bquote( paste( 
        .(spshort[j]), ', ', .(r), .(p))), 
             side = 3, adj = 0.1, cex = 0.7, line = -1.2)
      
      k = k+1
      if(j == i+1)
        mtext(text = spshort[i], side = 2, line = 1, cex = 0.8)
    }
  }
}

mtext(text = bquote( paste( "Species 1 ", pi[N]/pi[S])), side = 1, outer = T, cex = 1.3)
mtext(text = bquote( paste( "Species 2 ", pi[N]/pi[S])), side = 2, outer = T, cex = 1.3)

dump <- dev.off()

## no. of pairwise comparisions that are significantly positive
cobs <- as.data.frame(cobs)
colnames(cobs) <- c('sp1','sp2', 'cor','p')
cobs$cor <- as.numeric( cobs$cor)
cobs$p <- as.numeric( cobs$p)

cobs$p.adj <- p.adjust(p = cobs$p, method = 'holm')

npairs <- nrow(cobs)
nsigpos <- sum(cobs$cor > 0 & cobs$p.adj < 0.05)

cat(sprintf("%i out of %i pairwise comparisons of pN/pS are significantly positive at P < 0.05 with multiple-testing correction\n",nsigpos,npairs))
