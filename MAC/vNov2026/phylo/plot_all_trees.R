require(ape)
require(phytools)

indr <- 'njtrees'
fout <- 'trees.pdf'

fls <- list.files(indr)
nf <- length(fls)


pdf(fout)
par( mfrow = c(2,2), mar = c(0,0,1,0), oma = c(0,0,0,0), xpd = NA)

cxl = 0.6

ofs <- rep( 0.5, nf)
ofs[c(1,2,3,4,6)] <- c(1,1,1.5,1,1.5)

rts <- rep( 0, nf)
rts[c(1,2,4,6,7,8,10,11)] <- c(-45,-35,-15,45,30,60,30,20)

lims <- c(-4,4)

for(i in 1:nf){
  sp <- gsub(pattern = ".nwk", replacement = "", x = basename(fls[i]))
  tr <- read.tree(file = file.path( indr, fls[i]))
  if( i == 8){
    tr$edge.length[ tr$edge[,2] == 9] <- 1e-06
  }
  tr$edge.length <- tr$edge.length / max(tr$edge.length)  
  plot( tr, edge.width = 0.5, type = 'u', no.margin = T, show.tip.label = T, lab4ut = "axial", 
        label.offset = ofs[i], cex = cxl, rotate.tree = rts[i],
        x.lim = lims, y.lim = lims)
  mtext(text = paste("P.", sp), side = 3, line = -1, font = 3)
}


dump <- dev.off()