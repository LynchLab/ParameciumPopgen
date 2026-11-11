fin <- 'orthologs_species_pnps.csv'
fom <- 'ref_omega_M8.csv'
fout <- 'figure_5C.svg'

df <- read.csv(file = fin, row.names = 1)
sps <- colnames(df)
nsp <- ncol(df)

# keep only aurelias ( because these genes are always present in aurelias)
auredf <- df[ , grep("aure",colnames(df))]
mpnps <- apply( auredf, 1, function(x) mean(x[!is.na(x) & !is.infinite(x)]))

# load omega values (dN/dS)
dom <- read.csv(file = fom, row.names = 1)

genes <- intersect( names(mpnps), rownames(dom))

rpd <- data.frame( dnds = dom[ genes, 'omega'], pnps = mpnps[ genes])

# remove genes with pnps = 0
rpd <- rpd[ rpd$pnps > 0, ]

rpd$lpnps <- log10(rpd$pnps)
rpd$ldnds <- log10(rpd$dnds)

## linear regression on log-log scale
lob <- lm( lpnps ~ ldnds, data = rpd)
lob.sm <- summary(lob)

# goodness-of-fit
r2 <- round( lob.sm$adj.r.squared, 2)

# coefficients on log-10 scale
b0 <- round( lob.sm$coefficients[1,1], 2)
d0 <- round( lob.sm$coefficients[1,2], 2)
b1 <- round( lob.sm$coefficients[2,1], 2)
d1 <- round( lob.sm$coefficients[2,2], 2)

# coefficients on original scale
K <- round( 10^b0, 2)
Kmin <- round( 10^(b0-d0), 2)
Kmax <- round( 10^(b0+d0), 2)


## plotting
svg( fout, width = 3, height = 3)

#par( mar = c(5,4.5,2,2))

par( mar = c(3.3,4,1.5,0.8))

cxl = 1.0
cxx = 1.0
cxt = 0.6

ptype = 1

ltks <- -3:1
lims <- c(-2.5,0.5)

acol = 'darkslategrey'

tcol <- apply( col2rgb(col = acol), 2, function(x) 
  rgb( x[1], x[2], x[3], 100, maxColorValue = 255))

plot( lpnps ~ ldnds, data = rpd, col = tcol, cex = cxp, cex.lab = cxl,
      xlim = lims, ylim = lims, axes = F, cex.axis = cxx,
      xlab = '', ylab = '')

# regression line
abline( lob, lty = 2)
# bisector
abline(0,1,lty=3)

mtext( text = bquote( paste( beta, ' = ', .(b1), ' \u00B1', .(d1))), side = 3, adj = 1, cex = cxl)
mtext( text = bquote( paste( italic(R)^2, ' = ', .(r2))), side = 3, adj = 0, cex = cxl)

## axes
axis( side = 1, at = ltks, cex.axis = cxx,
      labels = as.expression( sapply( ltks, function(x) 
        bquote( paste( 10^.(x))))))

axis( side = 2, at = ltks, las = 2,  cex.axis = cxx,
      labels = as.expression( sapply( ltks, function(x) 
        bquote( paste( 10^.(x))))))
box()

op <- par( family = 'mono-space')

## axis labels
mtext( text = bquote( paste( "Mean ", pi[N]/pi[S], " in aurelia")), 
       side = 2, line = 2.4, cex = cxl)

mtext( text = bquote( paste( "Genus-wide ", italic(d)[N]/italic(d)[S])), 
       side = 1, line = 2.4, cex = cxl)

dump <- dev.off()