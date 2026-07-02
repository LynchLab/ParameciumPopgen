suppressMessages(require(RColorBrewer))

indr <- 'data/pca'
species <- c('biaurelia','tetraurelia','sexaurelia')

fin <- 'genolike_deg4.cov'
fpop <- 'sample_loclab.txt'
fplot <- 'images/figure_1E.png'

nsp <- length(species)

png( fplot, width = 2.966, height = 3.934, units = 'in', res = 150, pointsize = 7)
par( mfrow = c(3,1), mar = c(4.5,4.5,1,1), oma = c(0.5,0.5,0.5,9))

cxl = 1.4
cxx = 1.2
cxp = 2
ptype = 1

cols <- brewer.pal(n = 9, name = 'Set1')
insets <- c( -0.27, -0.39, -0.39)
slabs <- gsub("^","P. ",species)
los <- 0.05
  
for(i in 1:nsp){
  # load covariance matrix
  C <- as.matrix( read.table( file.path( indr, species[i], fin)))
  # do PCA
  eg <- eigen(C)

  # prop var explained
  evals <- eg$values
  pvars <- round( evals/sum(evals)*100, 1)

  pop <- read.csv( file.path( indr, species[i], fpop), sep=' ', header = F)

  pop$pc1 <- eg$vectors[,1]
  pop$pc2 <- eg$vectors[,2]

  pop <- pop[ pop$V2 != '', ]
  plabs <- levels( factor(pop$V2))
  pop$col <- cols[factor(pop$V2)]

  xlims <- range(pop$pc1)
  ylims <- range(pop$pc2)
  
  xlims[1] <- xlims[1] - los
  xlims[2] <- xlims[2] + los
  ylims[1] <- ylims[1] - los
  ylims[2] <- ylims[2] + los
  
  plot( pc2 ~ pc1, data = pop, col = pop$col, cex = cxp,
      xlim = xlims, ylim = ylims,
      lwd = 1, pch = c(1,2,5,6), 
      xlab = paste0( 'PC1 (',pvars[1],'% )'),
      ylab = paste0( 'PC2 (',pvars[2],'% )'),
      main = '', cex.lab = cxl, cex.axis = cxx)

  legend( "topright", inset = c( insets[i], -0.07), 
        xpd = NA, fill = cols, plabs, bty = 'n', cex = cxl)
  
  mtext(text = slabs[i], side = 3, line = -1.5, cex = 1.1, font = 3)
}

dump <- dev.off()