require(ape)

indr <- 'njtrees'

fls <- list.files(indr)

sps <- gsub(pattern = ".nwk", replacement = "", x = fls, fixed = T)

nsp <- length(sps)

trees <- lapply( fls, function(f) read.tree( file.path( indr, f)))
names(trees) <- sps

# sort trees by number of tips
trees <- trees[ order(sapply( trees, Ntip))]

# get lengths of internal branches
ibls <- lapply( trees, function(tr)
  tr$edge.length[ !apply( tr$edge, 1, function(x) any(x %in% 1:Ntip(tr)))])

libls <- sapply( ibls, log10)



boxplot( ibls)
abline( h = 0.008)

boxplot(libls)
 

## min-max scaling
mmibls <- sapply( ibls, function(x) ( x - min(x))/ ( max(x) - min(x)))
boxplot( mmibls)

## function to locate internal branches
find_ib <- function(tr){
  # plot to confirm labels
  plot( tr, use.edge.length = T)
  nodelabels(); tiplabels()
  # first N nodes are tips 
  # exclude edges contaning those
  ixs <- which( !apply( tr$edge, 1, function(x) any(x %in% 1:Ntip(tr))))
  # internal edge matrix
  emat <- tr$edge[ ixs, ]
  # internal edge lengths
  elen <- tr$edge.length[ ixs]
  # scale lengths by mean and SD
  sl <- (elen - mean(elen)) / sd(elen)
  # order internal branches in decreasing order of length
  ord <- order( sl, decreasing = T)
  cbind( emat[ ord, ], elen[ ord], sl) 
}

find_ib(trees$fokini)


### z-score ###############
png("figure_S2.png", width = 7, height = 5, units = 'in', res = 150)
par( mar = c(4.5,4.5,3,2))
cxt = 1.1

zibls <- sapply( ibls, function(x) (x - mean(x))/sd(x))
boxplot( zibls, names = F, boxwex = 0.7)
abline( h = 3, lty = 2)
slabs <- sapply( strsplit(x = names(trees), split = '', fixed = T), function(x)
  paste(c('P',x[1:4]), collapse = ''))
#mtext(text = slabs, side = 1, at = 1:nsp, line = 1, cex = 1)
text( 1:nsp, par("usr")[3] - 0.8, slabs, srt = 30, xpd = T, cex = cxt)
mtext( text = bquote( paste ( 'Scaled branch length, (',italic(b) - bar( italic(b)),') / ',sigma[italic(b)])), 
       side = 2, line = 2.5, cex = cxt)
dev.off()