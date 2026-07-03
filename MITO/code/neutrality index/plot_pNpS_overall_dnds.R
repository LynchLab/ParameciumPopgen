setwd('~/ASU Dropbox/Farhan ALI/Paramecium/intersp/final_315_samples/mito')

f1 <- 'dnds/overall/overall_dnds_gene.csv'
f2 <- 'pnps/pnps_mean_all_species.csv'

fout <- 'gene_species_pnps_overall_dnds.pdf'

df1 <- read.csv( file = f1, row.names = 1)
df2 <- read.csv( file = f2)

df2 <- df2[ df2$Mean_Ds > 0 & df2$Mean_Dn > 0, ]

genes <- intersect( rownames(df1), df2$Gene)
ng <- length(genes)

afunc <- function(gene){
  # extract
  gdf1 <- df1[ gene, 'omega']
  gdf2 <- df2[ df2$Gene == gene, ]
  # subset
  gsp <- gdf2$Species
  nsp <- length(gsp)
  ftdf <- as.data.frame( matrix( NA, nsp, 6, dimnames = list( 
    gsp, c('gene','species','pN','pS','pN/pS','dN/dS'))))
  ftdf$gene <- gene
  ftdf$species <- gsp
  for( j in 1:nsp){
    ftdf[ j, 3:6] <- c( gdf2$Mean_pN[j], gdf2$Mean_pS[j], gdf2$pN_over_pS[j], gdf1)
  }
  ftdf
}

ftls <- lapply( genes, function(g) tryCatch(expr = afunc(g),
                                            error = function(e) NA))
ftls <- do.call( rbind, ftls)
ftls <- ftls[ !apply( ftls, 1, anyNA), ]
df <- as.data.frame(ftls)
