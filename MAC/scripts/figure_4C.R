fin <- 'data/aurelia_core_orthologs_pnps.csv'
fom <- 'data/ref_omega_M8.csv'
fout <- 'images/figure_4C.png'

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

rpd$lpnps <- log10(rpd$pnps)
rpd$ldnds <- log10(rpd$dnds)

lob <- lm( lpnps ~ ldnds, data = rpd)
lob.sm <- summary(lob)

png( fout, width = 3.4, height = 3.4, units = 'in', res = 150, pointsize = 8)
par( mar = c(5,4.5,2,2))

cxl = 1.4
cxx = 1.3
cxp = 0.8
ptype = 1

ltks <- -3:1
lims <- c(-2.5,0.5)

plot( lpnps ~ ldnds, data = rpd, col = 'grey50', cex = cxp, cex.lab = cxl,
      xlim = lims, ylim = lims, axes = F, cex.axis = cxx,
      xlab = '', ylab = '')

abline( lob, lty = 2)
abline(0,1,lty=3)

# coefficients on log-10 scale
b0 <- round( lob.sm$coefficients[1,1], 2)
d0 <- round( lob.sm$coefficients[1,2], 2)
b1 <- round( lob.sm$coefficients[2,1], 2)
d1 <- round( lob.sm$coefficients[2,2], 2)

# coefficients on original scale
K <- round( 10^b0, 2)
Kmin <- round( 10^(b0-d0), 2)
Kmax <- round( 10^(b0+d0), 2)


mtext( text = bquote( paste( pi[N]/pi[S], ' = ', .(K), (italic(d)[N]/italic(d)[S])^.(b1))), side = 3, adj = 1, cex = cxl)


axis( side = 1, at = ltks, 
      labels = as.expression( sapply( ltks, function(x) 
        bquote( paste( 10^.(x))))))

axis( side = 2, at = ltks, 
      labels = as.expression( sapply( ltks, function(x) 
        bquote( paste( 10^.(x))))))
box()

op <- par( family = 'mono-space')

mtext( text = bquote( paste( pi[N]/pi[S], ", averaged over ", italic("P. aurelia"))), 
       side = 2, line = 3, cex = cxx)

mtext( text = bquote( paste( italic(d)[N]/italic(d)[S], ", estimated over all species")), 
       side = 1, line = 3, cex = cxx)

dump <- dev.off()