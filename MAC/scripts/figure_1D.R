# ---- setup & inputs ----
suppressPackageStartupMessages({
  library(dplyr)
  library(tools)
})

indr = 'data/genetic_vs_geo'
floc <- 'data/table_S1.csv'
fout <- 'images/figure_1D.png'

dloc <- read.csv(file = floc)
dloc <- dloc[ dloc$source != 'USA', ]
dloc <- dloc[ !((dloc$source == '') & (dloc$collection %in% c('Germany','Japan'))), ]
dloc <- dloc[ !is.na(dloc$latitude) & !is.na(dloc$longitude), ]

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
  df_list[[i]] <- d[ (d$Sample1 %in% dloc$sample) & (d$Sample2 %in% dloc$sample), ]
}

# ---- Spearman correlation per species (rho & p) ----
cors_tbl <- t( sapply( df_list, function(X) {
  ob <- suppressWarnings( cor.test( X$Geographic_Distance, X$Genetic_Distance, method = 's'));
  c( ob$estimate, ob$p.value);
}))
rownames(cors_tbl) <- names(df_list)
sigpos <- apply( cbind( round(cors_tbl[,1],2), ifelse(cors_tbl[,2] < 0.05, '**', '')), 1,
                 paste, collapse = ' ')

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

grp_dfs <- sapply( grps, function(x) df_list[x])


######################
## plotting
png( fout, width = 3.5, height = 3.5, units = 'in', res = 150, pointsize = 9)
par( mfrow = c(3,3), mar = c(0,0,0,0), oma = c(5.5,5,1,1))

cxl = 1.3
cxx = 1.3
cxp = 1.3

clrs <- read.csv('species_colors_pts.csv', row.names = 1)

# for each group
for( i in 1:ng){
  X <- grp_dfs[[i]]
  sps <- names(X)
  nx <- length(X)
  ylims <- range( unlist( lapply( X, function(x) x$Genetic_Distance)))
  # for each species in the group
  for(j in 1:nx){
    sX <- X[[j]]
    if(j == 1){
      plot( Genetic_Distance ~ Geographic_Distance, data = sX, 
            xlim = c(0,5), ylim = c(-4.5,-0.5), axes = F,
            col = clrs[ sps[j], 'color'], pch = clrs[ sps[j], 'pch'], xlab = '', ylab = '', xaxt = 'n')
    }
    else
      points( Genetic_Distance ~ Geographic_Distance, data = sX, 
              col = clrs[ sps[j], 'color'], pch = clrs[ sps[j], 'pch'])
  }
  if(i %in% 7:9){
    axis( side = 1, at = 0, outer = T, cex.axis = cxx)
    axis(side = 1, at = c(2,4), outer = T, cex.axis = cxx,  
         labels = as.expression( sapply( c(2,4), function(x) bquote( paste( 10^.(x))))))
  }
  if(i %in% c(1,4,7)){
    ytks <- c(-4:-1)
    axis( side = 2, at = ytks, cex.axis = cxx,
          labels = as.expression( sapply( ytks, function(x) bquote( paste( 10^.(x))))))
  }
  ## legend 
  # species labels
  splabs <- sapply( strsplit(x = sps, split = '', fixed = T), function(x)
    paste( c('P',x[1:4]),collapse = ''))
  clades <- sapply( strsplit(sps,split = "_clade", fixed = T), '[', 2)
  clades <- clades[!is.na(clades)]
  legs <- apply( cbind( splabs, clades, sigpos[sps]), 1, paste, collapse = ' ')
  legend( 'bottomright', legend = legs, 
          box.lty = 1, box.lwd = 0.5, box.col = 'grey40',
         col = clrs[ sps, 'color'], pch = clrs[ sps, 'pch'], cex = 1.15)
  box(col='grey40',lwd=0.8)
}

mtext(text = 'Geographic distance (km)', side = 1, outer = T, line = 3.5, cex = cxl)
mtext(text = 'Genetic distance (/bp)', side = 2, outer = T, line = 3, cex = cxl)
dump <- dev.off()