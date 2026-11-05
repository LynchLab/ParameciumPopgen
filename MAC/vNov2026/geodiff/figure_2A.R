# ---- setup & inputs ----
options(scipen = 999)

suppressPackageStartupMessages({
  library(dplyr)
  library(tools)
})


## inputs
indr = 'genetic_vs_geo'
floc <- 'isolate_location_coords.csv'
fcor1 <- 'mantel_indiv_spearman.csv'
fcor2 <- 'mantel_one_spearman.csv'

# output
fout <- 'figure_2A.png'

alpha <- 0.05

dloc <- read.csv(file = floc)
dloc <- dloc[ !is.na(dloc$latitude), ]

# ---- read & combine all species files; add species name from filename ----
files <- list.files(indr)
sps <- gsub('.csv', '', files)
nsp <- length(sps)

df_list <- sapply( sps, function(X) NULL)
for( i in 1:nsp){
  f <- file.path( indr, paste0(sps[i],'.csv'))
  d <- read.csv(f)
  # replace geo dist 0 with 1 ( for log transformation)
  d$Geographic_Distance[d$Geographic_Distance == 0] = 1
  
  d$Geographic_Distance <- log10(d$Geographic_Distance)
  d$Genetic_Distance <- log10(d$Genetic_Distance)
  # remove samples that are missing from location table
  df_list[[i]] <- d[ (d$Sample1 %in% dloc$isolate) & 
                       (d$Sample2 %in% dloc$isolate), ]
}

# ---- Load Spearman correlations per species (rho & p) ----
cors1 <- read.csv(fcor1, row.names = 1)
cors1$statistic <- round( cors1$statistic, 2)

cors2 <- read.csv( fcor2, row.names = 1)


## set up groups
grps <- list( bursaria = c('protobursaria', 'tritobursaria'),
              multimicronucleatum = c('multimicronucleatum', 'fokini'),
              caudatum = c('caudatum_cladeA', 'caudatum_cladeB'),
              primaaurelia = c( 'primaurelia_cladeA', 'primaurelia_cladeB', 
                                'triaurelia'),
              biaurelia = c( 'biaurelia', 'quindecaurelia'),
              tetraurelia = c('tetraurelia','sexaurelia'))

ng <- length(grps)

grp_dfs <- sapply( grps, function(x) df_list[x])


######################
## plotting
png( fout, width = 4, height = 2.6, units = 'in', res = 300, pointsize = 7)
par( mfrow = c(2,3), mar = c(0,0,0,0), oma = c(5.5,5.5,2,2))
cxl = 1.1
cxx = 1.5
cxp = 1
cxt = 1.2

clrs <- read.csv('../species_colors_pts.csv', row.names = 1)

nisos <- NULL

# for each group
for( i in 1:ng){
  X <- grp_dfs[[i]]
  sps <- names(X)
  nx <- length(X)
  ylims <- range( unlist( lapply( X, function(x) x$Genetic_Distance)))
  # start empty plot
  plot( NA, NA, xlim = c(0,5), ylim = c(-5,-1), axes = F,
        xlab = '', ylab = '')
  # grids
  abline(v = 0:5, lty = 3, col = 'grey')
  abline(h = -5:0, lty = 3, col = 'grey')
  # for each species in the group
  for(j in 1:nx){
    # no. of isolates
    nisos <- c( nisos, sum( dloc$species == sps[j]))
    sX <- X[[j]]
    points( Genetic_Distance ~ Geographic_Distance, data = sX, 
              col = clrs[ sps[j], 'color'], pch = clrs[ sps[j], 'pch'])
  }
  # x-axis
  if(i %in% 4:6){
    axis( side = 1, at = 0, outer = T, cex.axis = cxx)
    axis(side = 1, at = c(2,4), outer = T, cex.axis = cxx,  
         labels = as.expression( sapply( c(2,4), function(x) bquote( paste( 10^.(x))))))
  }
  # y-axis
  if(i %in% c(1,4)){
    ytks <- c(-4:-1)
    axis( side = 2, at = ytks, cex.axis = cxx, las = 2,
          labels = as.expression( sapply( ytks, function(x) bquote( paste( 10^.(x))))))
  }
  box(col='grey40',lwd=0.8)
}

snames <- unname(unlist( sapply( grp_dfs, names)))
names(nisos) <- snames

## legend
nsp <- length(snames)
# species labels
slabs <- gsub("_clade", " ", snames)


if(F){
legend( 'topright', inset = c(-1.1,-1.15), xpd = NA, 
        bty = 'n', title = "", text.font = 3, 
        cex = cxt, y.intersp = 1.5, ncol = 1,
        legend = slabs, 
        col = c( clrs[ snames, 'color'], rep(NA,nsp*2)), 
        pch = clrs[ snames, 'pch']
        )

legend( 'topright', inset = c(-1.8,-1.15), xpd = NA, 
        bty = 'n', title = "", text.font = 1, text.width = 0.5,
        cex = cxt, y.intersp = 1.5, ncol = 3,
        legend = c( sprintf( "%02d", nisos[ snames]), format( cors_tbl[ snames, "statistic"], nsmall = 2), 
                    cors_tbl[ snames, "pvalue"]), 
        col = NA, 
        pch = NA
)
mtext(text = expression( italic(n)), font = 3, outer = T, adj = 1.42, line = -0.3)
mtext(text = expression( italic("\U03C1")), font = 3, outer = T, adj = 1.52, line = -0.3)
mtext(text = expression( italic(P)), font = 3, outer = T, adj = 1.64, line = -0.3)
}


## axis labels
mtext(text = 'Geographic distance (km)', side = 1, outer = T, line = 3.5, cex = cxl)
mtext(text = 'Genetic distance (/bp)', side = 2, outer = T, line = 3.5, cex = cxl)
dump <- dev.off()