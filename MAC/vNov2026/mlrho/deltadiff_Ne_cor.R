f1 <- 'delta_pars_md250_pi4.csv'
f2 <- '../sitehets/species_pop_mu_Ne_maxpi4.csv'

minpi <- 0.004

fout <- 'delta_diff_Ne.svg'

d1 <- read.csv(file = f1, header = T, row.names = 3)

d1 <- d1[ d1$pi4 >= minpi, ]

d1 <- d1[ !is.na(d1$population), ]
d1$spnames <- apply( d1[ , c("species","population")], 1, paste, collapse="__")

d2 <- read.csv(file = f2)
d2 <- d2[ d2$population != 'all', ]
rownames(d2) <- apply( d2[ , c("species","population")], 1, paste, collapse="__")

d1$dd <- d1$delta3/d1$delta4

spdds <- split( d1$dd, d1$spnames)

dd <- as.data.frame( t( sapply( spdds, function(X) 
  c( length(X), min(X), mean(X), max(X))
  )))

colnames(dd) <- c("n","min","mean","max")

dd$Ne <- d2[ rownames(dd), 'Ne']
dd <- dd[ !is.na(dd$Ne), ]

dd$lNe <- log10(dd$Ne*1e07)
dd$lmin <- log10(dd$min) 
dd$lmean <- log10(dd$mean) 
dd$lmax <- log10(dd$max) 

lob <- lm( lmean ~ lNe, data = dd)
lob.sm <- summary(lob)
b <- round( lob.sm$coefficients[2,1], 2)
d <- round( lob.sm$coefficients[2,2], 2)
r2 <- round( lob.sm$adj.r.squared, 2)
regstats <- c(b,d)

cob <- cor.test( dd$mean, dd$Ne, method = 's')
r <- round( cob$estimate, 2)
p <- round( cob$p.value, 4)
corstats <- c(cob$estimate, cob$p.value)

### plotting ##############################
svg( fout, width = 2, height = 2)
par( mar = c(3.3,3.1,1.2,0.6))

cxl = 1.0
cxx = 1.0
cxp = 1.2
cxt = 0.9

# Set color palette
fcol <- '../../species_colors_pts.csv'
cols <- read.csv( file = fcol, row.names = 1)

dd$sp <- gsub("__pop.*", "", rownames(dd))
dd$cols <- cols[ dd$sp, 'color']

# points
dd$pch <- cols[ dd$sp, 'pch']

# line types
ltypes <- c(2,3)

# axis ticks
lxtks <- 6:9
ytks <- c(1,2,3,5)
lytks <- log10(ytks)

plot( NA, NA, 
      xlim = c(6.5,9), ylim = range(lytks), 
      xlab = '',
      ylab = '', 
      axes = F)

abline( 0, 1, lty=3, col='grey40')
abline( lob, lty = 2, col='grey40')

## segments for n > 1
subdd <- dd[ dd$n > 1, ]
segments(x0 = subdd$lNe, y0 = subdd$lmin, x1 = subdd$lNe, y1 = subdd$lmax, lty = 1, 
         lwd = 1.5, col = subdd$cols)

# scatter of points
points( lmean ~ lNe, data = dd, col = dd$cols, cex.lab = cxl, 
          pch = dd$pch, cex = cxp)
  
## axes
axis( side = 1, at = lxtks, labels = as.expression( 
  sapply( lxtks, function(x) bquote( paste( 10^.(x))))), cex.axis = cxx)
axis( side = 2, at = lytks, labels = ytks, cex.axis = cxx, las = 2)
box()

## axis labels
mtext( text = expression(italic(N)[e]), side = 1, line = 2.4, cex = cxl)
mtext( text = bquote( paste( Delta[3], ' / ', Delta[4])),
       side = 2, line = 1.8, cex = cxl)


## stats
mtext( text = bquote( paste( beta, ' = ', .(regstats[1]), 
                             ' \u00B1 ', .(regstats[2]))),
       lty = ltypes, bty = 'n', cex = cxt, adj = 1)

dev.off()