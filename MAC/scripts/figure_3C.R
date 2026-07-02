frh <- 'data/summary_mlrho_least_squares.csv'
fpi <- 'data/species_pop_class_het_gfe.csv'
fout <- 'images/figure_3C.png'

## mlrho delta v/s dist profile least-squares fit summary stats of parameters
drh <- read.csv( file = frh)
drh <- drh[ drh$stats == '50%', ]
drh <- drh[ drh$params %in% c('rho1', 'cbyu'), ]

drh <- drh[,c("species","params","value")]
drh <- reshape(data = drh, idvar = 'species', timevar = 'params', 
              v.names = 'value', direction = 'wide')
rownames(drh) <- drh$species
drh <- drh[ , -1]
colnames(drh) <- c('cu','rho1')

## p0p4 data for pop
dpi <- read.csv( file = fpi)
dpi <- dpi[ dpi$population != 'all', ]
rownames(dpi) <- apply( dpi[ c('species','population')], 1, paste, collapse = '_')

## merge data
species <- intersect( rownames(drh), rownames(dpi))
df <- data.frame( cu = drh[species,'cu'], rho1 = drh[species,'rho1'], p0p4 = dpi[species,'p1p4'])
rownames(df) <- species
df <- df[ order(df$cu), ]

## take logs
df$lcu <- log10(df$cu)
df$lrho1 <- log10(df$rho1)
df$lp0p4 <- log10(df$p0p4)

## spearman correlation
cob <- cor.test( df$cu, df$p0p4, method = 's')

## linear regression
lob <- lm( lp0p4 ~ lcu, data = df)
lob.sm <- summary(lob)

## plotting
png( fout, width = 3, height = 3, units = 'in', res = 150, pointsize = 8)
par( mar = c(4.5,4.5,2,1))

## graphical parameters

# Set color palette
fcol <- 'species_colors_pts.csv'
cols <- read.csv( file = fcol, row.names = 1)

df$sp <- gsub("_pop.*", "", rownames(df))
df$cols <- cols[ df$sp, 'color']

# points
df$pch <- cols[ df$sp, 'pch']

# font sizes
cxl = 1.3
cxx = 1.2
cxp = 1.3


# axis limits
ytks <- c(0.02,0.05,0.2,0.5)
lytks <- log10(ytks)

lxtks <- seq( -2, 0, 1)

# scatter
plot( lp0p4 ~ lcu, data = df, col = df$cols, 
      xlim = range(lxtks), ylim = range(lytks), 
      xlab = expression(italic(c)/mu),
      ylab = expression(pi[0]/pi[4]),
      axes = F, cex.lab = cxl, pch = df$pch, cex = cxp, lwd = 1.2)

# regression line
abline( lob, lty = 2,col = 'grey40')

# regression stats
b <- round( lob.sm$coefficients[2,1], 2)
d <- round( lob.sm$coefficients[2,2], 2)
r2 <- round( lob.sm$adj.r.squared, 2)
mtext( text = bquote( paste( beta, ' = ', .(b), ' (', '\u00B1', .(d), ')')), side = 3, adj = 1, cex = cxx)

# correlation stats
r <- round( cob$estimate, 2)
p <- round( cob$p.value, 4)
mtext( text = bquote( paste( italic(r), ' = ', .(r))), side = 3, adj = 0, cex = cxx)

axis( side = 1, at = lxtks, labels = as.expression( 
  sapply( lxtks, function(x) bquote( paste( 10^.(x))))), cex.axis = cxx)
axis( side = 2, at = lytks, labels = ytks, cex.axis = cxx)

box()
dump <- dev.off()
