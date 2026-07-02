fin <- 'data/pnps_reldif_aure_burs_ncopies.csv'
fout <- 'images/figure_4A.png'

df <- read.csv( file = fin)
df <- df[ df$ngenes <= 8, ]


lsdf <- split( df$difrelburs, df$ngenes)

class_size <- sapply( lsdf, length)


png( fout, width = 3.4, height = 3.4, units = 'in', res = 150, pointsize = 8)
par( mar = c(4.5,5,2,2))

cxl = 1.4
cxx = 1.3
cxp = 0.1
ptype = '.'
acol = 'slategrey'
tcol <- apply( col2rgb(col = acol), 2, function(x) rgb( x[1], x[2], x[3], 100, maxColorValue = 255))

boxplot( lsdf, ylim = c(-5,10), notch = T, boxwex = 0.6,
      xlab = bquote( paste( "Number of gene copies, ",italic(n))),
      ylab = bquote( paste( (bar(pi[N]/pi[S]))[aure]/(pi[N]/pi[S])[burs]-1)), 
      cex.lab = cxl, col = tcol, pch = ptype, cex = cxp, cex.axis = cxx)
abline( h = 0, lty = 2)

mtext( text = class_size[1:8], side = 1, line = -1.2, at = 1:8, cex = 1)

ytks <- seq(-5,10,1)
axis(side = 2, at = ytks, labels = F, lwd.ticks = 0.5)

dump <- dev.off()