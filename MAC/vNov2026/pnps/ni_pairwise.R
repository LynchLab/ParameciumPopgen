fin <- 'orthologs_species_pnps.csv'
fom <- '~/github/ParameciumPopgen/MAC/data/constraints_across_species/gene_pairwise_dN_dS.csv'
fout <- 'gene_sp_ni_pairwise.csv'

maxds = 0.8

# load pnps for single-copy orthologs
df <- read.csv(file = fin)
colnames(df)[grep("protobursaria",colnames(df))] <- "bursaria"
sps <- colnames(df)[-1]
nsp <- length(sps)
colnames(df)[1] <- "gene"

# 4-letter version of species names
spshort <- sapply( strsplit(x = sps, split = '', fixed = T), function(x) paste(x[1:4],collapse=''))
colnames(df)[-1] <- spshort

# convert to a long table of gene, species, pnps
longdf <- reshape( df, varying = spshort, v.names = "pnps", timevar = "species", idvar = "gene", times = spshort, direction = "long")
rownames(longdf) <- NULL
# only keep gene-species pairs with pnps > 0 & not infinite
longdf <- longdf[ ( longdf$pnps > 0 & !is.infinite(longdf$pnps)), ]

### load dnds values
dom <- read.csv( file = fom)
cat( sprintf( "dNdS data has %i entries\n",nrow(dom)))
## remove comparisons with dS > maxds
dom <- dom[ dom$dS <= maxds, ]
cat( sprintf( "dS <= %.1f for %i entries\n",maxds,nrow(dom)))
## both dN and dS must be > 0
dom <- dom[ (dom$dS>0) & (dom$dN>0), ]
cat( sprintf( "dN and dS > 0 for %i entries\n",nrow(dom)))

## limit data to species present in both
comsps <- intersect( unique(longdf$species), unique( c( dom$species1, dom$species2)))
cat( sprintf("Species with valid pN/pS and dN/dS data = %s\n", 
             paste(comsps,collapse = ', ')))

longdf <- longdf[ longdf$species %in% comsps, ]
dom <- dom[ (dom$species1 %in% comsps) & (dom$species2 %in% comsps), ]

## limit data to genes present in both
genes <- intersect( longdf$gene, dom$gene)
ng <- length(genes)
cat( sprintf( "# genes with both pNpS and dNdS data = %i\n", ng))
longdf <- longdf[ longdf$gene %in% genes, ]
dom <- dom[ dom$gene %in% genes, ]

cat( sprintf( "%i entries have their gene species pairs in pNpS data\n",nrow(dom)))

### extract pnps values for each row in the pairwise dnds table
afunc <- function(row){
  gene <- row[1]
  sp1 <- row[2]
  sp2 <- row[3]
  sub <- longdf[ (longdf$gene == gene), ]
  out <- sub[ match( c(sp1,sp2), sub$species), 'pnps']
  return(out)
}

pnps <- t( apply( dom, 1, afunc))
dom$pnps1 <- pnps[,1]
dom$pnps2 <- pnps[,2]

dom <- dom[ !apply( dom, 1, anyNA), ]
cat( sprintf( "%i entries have pNpS values\n",nrow(dom)))

dom$ni1 <- dom$pnps1/dom$dN.dS
dom$ni2 <- dom$pnps2/dom$dN.dS

ls_gene_dom <- split( dom, dom$gene)
genes <- names(ls_gene_dom)
ng <- length(genes)
cat( sprintf( "These come from %i genes\n",ng))
out <- NULL

for(i in 1:ng){
  gene_dom <- ls_gene_dom[[i]]
  sps <- unique(c(gene_dom$species1,gene_dom$species2))
  nsp <- length(sps)
  for(j in 1:nsp){
    sp_dom <- gene_dom[ (gene_dom$species1 == sps[j]) | (gene_dom$species2 == sps[j]), ]
    sp_dom <- sp_dom[ order(sp_dom$dS), ]
    minds_row <- t( sp_dom[ 1, ])
    sp_ix <- which(minds_row == sps[j])-1
    ot_ix <- ifelse( sp_ix == 1, 2, 1)
    ot_sp <- minds_row[ paste0("species",ot_ix), 1]
    sp_pnps_ni <- minds_row[ paste0( c( "pnps", "ni"), sp_ix), 1]
    entry <- c( genes[i], sps[j], ot_sp, minds_row[4:6,1], sp_pnps_ni)
    out <- rbind( out, entry)
  }
}

colnames(out) <- c('gene','focal_species','closest_species','dS','dN','dN/dS','pN/pS','NI')
rownames(out) <- NULL
out <- as.data.frame(out)
out$NI <- as.numeric(out$NI)
nni <- nrow(out)
cat(sprintf("# gene-species combinations with NI calculated = %i\n",nni))
write.csv( x = out, file = fout, quote = F, row.names = F)

## instances of neutral or purifying selection
nneg <- sum(out$NI > 1)
pneg <- nneg/nni * 100

cat(sprintf("# gene-species combinations with NI > 1 = %i / %i (%.1f %%)\n",
            nneg,nni,pneg))


### get quantiles of NI for each gene
gqts <- as.data.frame( t( sapply( split( out$NI, out$gene), quantile, 
                                  probs = c(0,0.05,0.5,0.95,1), na.rm=T)))


## genes with all NI < 1
allpos <- sapply( split( out$NI, out$gene), function(x) all(x<1))
ng <- length(allpos)
allpos <- names(allpos)[allpos]
npos <- length(allpos)
P <- npos/ng
cat(sprintf("%i out of %i genes have their 95 percentile NI < 1 i.e P = %f\n",npos,ng,P))
gfun <- read.csv("~/github/ParameciumPopgen/MAC/data/constraints_across_species/gene_cat_class.tsv",sep='\t')
gfun <- gfun[ gfun$class != "uncategorized", ]
posfun <- allpos[ allpos %in% gfun$gene]
cat(sprintf("%i of these are categorized by function\n",length(posfun)))
pos_fun <- gfun[ match( posfun, gfun$gene), 1:2]
write.table(x = pos_fun, file = "", quote = F, sep='\t', row.names = F, col.names = F)