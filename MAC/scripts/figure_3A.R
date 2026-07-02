indr1 <- "data/avg_delta_rho/"
fpm <- 'data/least_squares_avg_delta.csv'
fout <- "images/figure_3A.png"
maxd <- 5000

### Lynch et al. 2014 -- Expected Delta given distance d and parameters

expected_delta <- function(d, theta, rho1, x, L) {
  # Eq. 13a: distance-dependent recombination rate (population-scaled)
  rho_d <- rho1 * (x * d + (1 - x) * L * (1 - exp(-d / L)))
  # Eq. 8: expected correlation of zygosity
  theta * (18 + rho_d) / (18 + 13 * rho_d + rho_d^2)
}

## summary stats of estimated parameters from the least-squares fit
dpm <- read.csv(file = fpm)
dpm <- dpm[ grep("pop",dpm$X,invert = T), ]
rownames(dpm) <- dpm$X
dpm <- dpm[ , -1]

### average delta v/s distance profiles
fls <- list.files(indr1)
fls <- grep("pop",fls,invert = T,value = T)
sps <- gsub(pattern = ".csv", replacement = "", x = fls, fixed = T)
nsp <- length(sps)
dfs <- sapply( sps, function(x) NULL)

for(i in 1:nsp){
  df <- read.csv( file.path( indr1, fls[i]))
  df <- df[ df$dist <= maxd, ]
  df$ldist <- log10(df$dist)
  df$ldelta <- log10(df$delta)
  dfs[[i]] <- df[ order(df$dist), ]
}

## set up groups
grps <- list( bursaria = c( 'bursaria','tribursaria'), 
              multimicronucleatum = c( 'multimicronucleatum','fokini'),
              caudatum = c( 'caudatum_cladeA','caudatum_cladeB'),
              primaurelia = c( 'primaurelia_cladeA','primaurelia_cladeB'),
              biaurelia = 'biaurelia',
              quindecaurelia = 'quindecaurelia',
              triaurelia = 'triaurelia',
              tetraurelia = 'tetraurelia',
              sexaurelia = 'sexaurelia')

ng <- length(grps)

grp_dfs <- sapply( grps, function(x) dfs[x])

######################

### plotting ############
png(fout, width = 6, height = 5, units = 'in', res = 150)
par( mfrow = c(3,3), mar = c(0,0,0,0), oma = c(5.1,5.1,1,1))

## graphical parameters
cxl = 1.2
cxx = 1.3
cxp = 1.3

# Set color palette
fcol <- 'species_colors_pts.csv'
dcol <- read.csv(file = fcol,row.names = 1)
cols <- dcol[ sps, ]

# axis ticks
xtks <- 0:4
ytks <- -3:0

# for each group of species
for( i in 1:ng){
  X <- grp_dfs[[i]]
  # number of species in the group
  nx <- length(X)
  plot( NA, NA, xlim = range(xtks), ylim = range(ytks), axes = F, 
        xlab = '', ylab = '', type = 'n')
  # grids
  abline( v = xtks, lty = 3, col = 'grey50', lwd = 0.5)
  abline( h = ytks, lty = 3, col = 'grey50', lwd = 0.5)
  box()
  # for each species
  for(j in 1:nx){
    # get species data
    sX <- X[[j]]
    # species name
    sp <- names(X)[j]
    # plot scatter of species data points
    points( ldelta ~ ldist, data = sX, col = cols[ sp, 'color'], pch = cols[ sp, 'pch'])
    # predicted delta based on fit
    prdelta <- expected_delta(d = sX$dist, theta = dpm[ sp, "theta_fit"], 
                              rho1 = dpm[ sp, "rho1"], x = dpm[ sp, "x"], L = dpm[ sp, "avg_L"])
    lines( sX$ldist, log10(prdelta), 
           lty = 1, col = 'darkslategrey', lwd = 1.3)
  }
  #abline( v = log10(400), lty = 2)
  # legend
  spnames <- paste("P. ", names(X))
  spnames <- gsub("_clade", " clade ", spnames)
  legend('topright', legend = spnames, text.font = 3, bty = 'n', 
         cex = cxl, col = cols[ names(X), 'color'], pch = cols[ names(X), 'pch'])
  # x-axis
  if(i %in% 7:9){
    axis(side = 1, at = xtks[-1], outer = T, cex.axis = cxx,  
         labels = as.expression( sapply( xtks[-1], function(x) bquote( paste( 10^.(x))))))
  }
  # y-axis
  if(i %in% c(1,4,7)){
    axis( side = 2, at = ytks[-1], cex.axis = cxx,
          labels = as.expression( sapply( ytks[-1], function(x) bquote( paste( 10^.(x))))))
  }
}

## axis labels
mtext(text = bquote( paste( 'Distance between sites (bp), ',italic(d))), side = 1, cex = cxl, line = 4, outer = T)
mtext(text = bquote( paste( "Zygosity correlation, ", Delta)), 
      side = 2, cex = cxl, line = 3.2, outer = T)
box()

#

## legend
if(F){
  op <- par(family='mono-space')
  legend("topright", inset = c(-0.37,-1.1), xpd = NA, text.font = 3, 
         legend = cols$sp[1:7], ncol = 1, bty = 'n', 
         col = cols$hex[1:7], pch = cols$pch[1:7], cex = 1.1, y.intersp = 1.2)
  
  legend("topright", inset = c(-0.48,0.01), xpd = NA, text.font = 3, 
         legend = cols$sp[8:13], ncol = 1, bty = 'n', 
         col = cols$hex[8:13], pch = cols$pch[8:13], cex = 1.1, y.intersp = 1.2)
}
dump <- dev.off()