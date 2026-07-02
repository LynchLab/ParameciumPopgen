#### Are pNpS of genes correlated across species? ####
source('scripts/sci_notation.R')

fin <- 'data/aurelia_core_orthologs_pnps.csv'
fout <- 'images/figure_4B.png'

df <- read.csv(file = fin, row.names = 1)

sps <- colnames(df)
nsp <- length(sps)

spshort <- sapply( strsplit(x = sps, split = '', fixed = T), function(x) paste(x[1:4],collapse=''))

cobs <- NULL

# for each unique pair of species
for(i in 1:(nsp-1)){
  for(j in 2:nsp){
    if( j <= i) next
    # take the subset from the full dataframe
    sub <- df[ , c(i,j)]
    # take genes with data for both species, excluding genes with undefined pNpS
    sub <- sub[ !apply( sub, 1, function(x) any(is.infinite(x) | is.na(x))), ]
    # spearman correlation
    scorob <- suppressWarnings( cor.test(x = sub[,1], y = sub[,2], method = 's'))
    cobs <- rbind( cobs, c( sps[ c(i,j)], scorob$estimate, scorob$p.value))
  }
}

## no. of pairwise comparisons that are significantly positive
cobs <- as.data.frame(cobs)
colnames(cobs) <- c('sp1','sp2', 'cor','p')
cobs$cor <- as.numeric( cobs$cor)
cobs$p <- as.numeric( cobs$p)

cobs$p.adj <- p.adjust(p = cobs$p, method = 'holm')

npairs <- nrow(cobs)
nsigpos <- sum(cobs$cor > 0 & cobs$p.adj < 0.05)

cat(sprintf("%i out of %i pairwise comparisons of pN/pS are significantly positive at P < 0.05 with multiple-testing correction\n",nsigpos,npairs))

### aurelia v/s others
# one species must be aurelia
acobs <- cobs[ grepl("aure", cobs$sp1) | grepl("aure", cobs$sp2), ]
# within aurelia
wnxs <- which( grepl("aure", acobs$sp1) & grepl("aure", acobs$sp2))
wncobs <- acobs[ wnxs, ]
# aurelia with others
otcobs <- acobs[ !(1:nrow(acobs)) %in% wnxs, ]

wncors <- wncobs$cor
otcors <- otcobs$cor

n1 <- length(wncors)
n2 <- length(otcors)

## wilcox test
wob <- wilcox.test( x = wncors, y = otcors, alternative = 't')

## plotting
png( fout, width = 2.5, height = 3.4, units = 'in', res = 150, pointsize = 8)
par( mar = c(4,4.5,2,2))

cxl = 1.2
cxx = 1.2
cxp = 1.1

p <- scinot(wob$p.value)

boxplot( wncors, otcors, notch = F, boxwex = 0.5, 
         names = NA, ylim = c(0.1,0.65), yaxt = 'n',
         ylab = 
           bquote( paste( 
             "Pairwise Spearman correlation in ", pi[N]/pi[S])),
         cex.lab = cxl, cex.axis = cxx)

# box labels
mtext(text = "Among", side = 1, at = 1, line = 1, cex = cxx)
mtext(text = expression(italic("P. aurelia")), side = 1, at = 1, line = 2, cex = cxx)
mtext(text = expression(italic("P. aurelia")), side = 1, at = 2, line = 1, cex = cxx)
mtext(text = "with others", side = 1, at = 2, line = 2, cex = cxx)
# y-axis
axis(side = 2, at = seq(0.2,0.6,0.2), cex.axis = cxx)
axis(side = 2, at = c(0.1,0.3,0.5), labels = F, lwd.ticks = 0.5, cex.axis = cxx)
# sample size
#mtext( text = bquote( paste( italic(n), ' = ')), side = 3, line = -2, at = 0.5, 
#       cex = cxx, adj = 0)
mtext( text = c(n1,n2), side = 1, line = -1, at = 1:2, cex = cxx)
# p-value
mtext(text = bquote( paste( italic(P), ' = ', .(p[1])^.(p[2]))), adj = 1, 
      cex = cxx, line = 0.1)
dump <- dev.off()