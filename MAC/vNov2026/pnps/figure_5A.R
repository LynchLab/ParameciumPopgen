fin <- 'pnps_reldif_aure_burs_ncopies.csv'
fout <- 'figure_5A.svg'

df <- read.csv( file = fin)
df <- df[ df$ngenes <= 8, ]


lsdf <- split( df$difrelburs, df$ngenes)

class_size <- sapply( lsdf, length)


svg( fout, width = 3.5, height = 3)

par( mar = c(3.3,4,1.5,1))

cxl = 1.0
cxx = 1.0
cxt = 0.6

ptype = '.'

acol = 'darkslategrey'

tcol <- apply( col2rgb(col = acol), 2, function(x) rgb( x[1], x[2], x[3], 100, maxColorValue = 255))

## boxes
boxplot( lsdf, ylim = c(-2,7), notch = T, boxwex = 0.6,
      xlab = "",
      ylab = "", 
      axes = F, outline = T, border = NA,
      col = tcol, cex = cxp)

## medians
meds <- sapply( lsdf, median)
points( x = 1:8, y = meds, pch = '_', cex = 1.2)

## ranges
qts <- sapply( lsdf, function(X) quantile(x = X, probs=c(0,0.25,0.75,0.9)))
arrows( x0 = 1:8, y0 = qts[3,], x1 = 1:8, y1 = qts[4,], col = acol, 
        angle = 90, length = 0.05)
arrows( x0 = 1:8, y0 = qts[2,], x1 = 1:8, y1 = qts[1,], col = acol,
        angle = 90, length = 0.05)

## zero
abline( h = 0, lty = 3, lwd = 0.6)

mtext( text = class_size[1:8], side = 1, line = -1.2, at = 1:8, cex = cxt)

xtks <- 1:8
axis(side = 1, at = xtks)


ytks <- seq(0,6,2)
axis(side = 2, at = ytks, las = 2)
ytks <- seq(-1,6,1)
axis(side = 2, at = ytks, labels = F, tcl=0.2)

box()

mtext( text = "Number of gene copies in aurelia",
       side = 1, line = 2, cex = cxl)

mtext( text = bquote( paste( 
  pi[N]/pi[S], ' in aurelia')), 
       side = 2, line = 2.6, cex = cxl)
mtext( text = bquote( paste( 
  'relative to bursaria')), 
  side = 2, line = 1.7, cex = 0.9)




dump <- dev.off()