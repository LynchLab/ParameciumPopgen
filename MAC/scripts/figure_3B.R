f1 <- 'data/species_pop_class_het_gfe.csv'
indr <- 'data/avg_delta_rho'
fout <- 'images/figure_3B.png'

## extract delta at d = 1,10,100,1000
dists <- 10^c(0,3)
fls <- list.files( indr, pattern = 'pop')
pops <- gsub( ".csv", "", fls)
nf <- length(fls)
dd <- matrix( NA, nf, length(dists), dimnames = list( pops, paste0("d",dists)))

for(i in 1:nf){
  dl <- read.csv( file = file.path( indr, fls[i]), row.names = 1)
  dd[i,] <- dl[ dists, 'delta']
}

dd <- as.data.frame( log10(dd))

## extract p0/p4
df <- read.csv(file = f1)
df$spnames <- apply( df[ , c('species','population')], 1, paste, collapse="_")
df$pnps <- log10(df$p1p4)

## merge
dd$pnps <- df$pnps[ match( rownames(dd), df$spnames)]
dd <- dd[ !apply( dd, 1, anyNA), ]

png( fout, width = 3, height = 3, units = 'in', 
     res = 150, pointsize = 8)
par( mar = c(4.5,4.5,2.7,1))

cxl = 1.3
cxx = 1.2
cxp = 1.3

# Set color palette
fcol <- 'species_colors_pts.csv'
cols <- read.csv( file = fcol, row.names = 1)

dd$sp <- gsub("_pop.*", "", rownames(dd))
dd$cols <- cols[ dd$sp, 'color']

# points
dd$pch <- cols[ dd$sp, 'pch']

# line types
ltypes <- c(2,3)

# positions for slopes to appear
pos <- rbind( c(-1.5,log10(0.04)), 
              c(-2.7,log10(0.3)))

# hold regression stats
regstats <- matrix( NA, 2, 2)
# hold correlation stats
corstats <- matrix( NA, 2, 2)

lxtks <- seq( -3, 0, 1)
ytks <- c(0.02,0.05,0.2,0.5)
lytks <- log10(ytks)

plot( NA, NA, 
      xlim = range(lxtks), ylim = range(lytks), 
      xlab = expression(Delta[italic(d)]),
      ylab = expression(pi[0]/pi[4]), cex.lab = cxl,
      axes = F)

for( i in 1:length(dists)){
  # distance column name
  dname <- colnames(dd)[i]
  # scatter of points
  points( dd[,dname], dd[,'pnps'], col = dd$cols, 
          cex.lab = cxl, pch = dd$pch, cex = cxp, lwd = 1.2)
  # linear regression
  lob <- lm( dd[,'pnps'] ~ dd[,dname] )
  abline( lob, lty = ltypes[i], col = 'grey40')
  lob.sm <- summary(lob)
  b <- round( lob.sm$coefficients[2,1], 2)
  d <- round( lob.sm$coefficients[2,2], 2)
  r2 <- round( lob.sm$adj.r.squared, 2)
  regstats[i,] <- c(b,d)

  # spearman correlation
  cob <- cor.test( dd[,dname], dd$pnps, method = 's')
  r <- round( cob$estimate, 2)
  p <- round( cob$p.value, 4)
  corstats[i,] <- c(cob$estimate, cob$p.value)
}

## axes
axis( side = 1, at = lxtks, labels = as.expression( 
  sapply( lxtks, function(x) bquote( paste( 10^.(x))))), cex.axis = cxx)
axis( side = 2, at = lytks, labels = ytks, cex.axis = cxx)
## legend
  legend('topright', inset = c(-0.03,-0.2), xpd = T, legend = as.expression( sapply( 
  1:length(dists), function(x) bquote( 
    paste(italic(d), ' = ', .(dists[x]), ', ', 
          beta, ' = ', .(regstats[x,1]), ' \u00B1 ', .(regstats[x,2]))))), 
  lty = ltypes, bty = 'n', cex = cxx)

box()
dump <- dev.off()