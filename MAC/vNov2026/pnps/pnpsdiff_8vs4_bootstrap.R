source("../../scripts/sci_notation.R")

fin <- 'pnps_reldif_aure_burs_ncopies.csv'
fout <- 'figure_S7.png'

set.seed(11211)

reps <- 1000000

df <- read.csv( file = fin)

df4 <- df$difrelburs[ df$ngenes == 4] 
df8 <- df$difrelburs[ df$ngenes == 8] 

n4 <- length(df4)
n8 <- length(df8)

meds4 <- replicate( n = reps, expr = median( 
  sample( x = df4, size = n8, replace = T)))

med8 <- median(df8)

P <- signif( (sum(meds4 <= med8)+1)/(reps+1), digits = 2)

ps <- scinot(P)

wob <- wilcox.test( x = meds4, mu = med8)
wob

tob <- t.test( x = meds4, mu = med8)
tob

png( fout, width = 4, height = 4, units = 'in', res = 150, 
     pointsize = 8)
par( mar = c(5,5,2,2))

cxl = 1.3
cxx = 1.2
cxt = 1.1

ptype = '.'
acol = 'slategrey'
tcol <- apply( col2rgb(col = acol), 2, function(x) rgb( x[1], x[2], x[3], 100, maxColorValue = 255))

hist( meds4, border = NA, xlim = c(-1,5), probability = T,
      main = '', cex.axis = cxx,
      xlab = bquote( paste( (bar(pi[N]/pi[S]))[aure]/(pi[N]/pi[S])[burs]-1)),
      cex.lab = cxl)

abline( v = med8, lty = 2)

mtext(text = bquote( paste( italic(P), ' = ', .(ps[1])^.(ps[2]))), 
      side = 3, adj = 1, cex = cxt)

box()

dev.off()
