setwd('~/ASU Dropbox/Ravina Telkar/intersp/final_315_samples/mito')

f1 <- 'dnds/ancestral/ancestral_dnds_highprob_filterallsites_0.60.csv'
f2 <- 'pnps/pnps_filtered_sites/pnps_mean_high_prob_0.6_sites.csv'
f3 <- 'dnds/pairwise/closest_pairs_dS.csv'

fout <- 'gene_species_anc_DnDs_NI.csv'

df1 <- read.csv( file = f1)
df2 <- read.csv( file = f2)
df3 <- read.csv( file = f3)

## select genes and species with pairwise dS < 1
df3 <- df3[ df3$dS < 1, ]

df1 <- df1[ apply( df1, 1, function(x) 
  any( df3$gene == x[1] & (df3$sp1 == x[2] | df3$sp2 == x[2]))), ]

# fix species names in df2
df2$Species <- sapply( strsplit( df2$Species, split = '', fixed = T),
                       function(x) paste(x[1:4],collapse=''))

species <- intersect( df1$Species, df2$Species)
nsp <- length(species)
df1 <- df1[ df1$Species %in% species, ]
df2 <- df2[ df2$Species %in% species, ]

df1 <- df1[ order(df1$Species), ]
df2 <- df2[ order(df2$Species), ]

genes <- intersect( df1$Gene, df2$Gene)
ng <- length(genes)

afunc <- function(gene){
  # extract
  gdf1 <- df1[ df1$Gene == gene, ]
  gdf2 <- df2[ df2$Gene == gene, ]
  # filter 
  gdf1 <- gdf1[ gdf1$ds > 0 & gdf1$dn > 0, ]
  gdf2 <- gdf2[ gdf2$Mean_Ds > 0 & gdf2$Mean_Dn > 0, ]
  # subset
  gsp <- intersect( gdf1$Species, gdf2$Species)
  nsp <- length(gsp)
  gdf1 <- gdf1[ match( gsp, gdf1$Species), ]  
  gdf2 <- gdf2[ match( gsp, gdf2$Species), ]
  ftdf <- as.data.frame( matrix( NA, nsp, 13, dimnames = list( 
    gsp, c('gene','species','pN','pS','dN','dS','pN/pS','dN/dS','obs NI','est NI','low','high','p'))))
  ftdf$gene <- gene
  ftdf$species <- gsp
  for( j in 1:nsp){
    dvec <- c( gdf2$Mean_Dn[j], gdf2$Mean_Ds[j], gdf1$dn[j], gdf1$ds[j])
    dmat <- matrix( dvec, 2, 2, byrow = T)
    pnps <- dvec[1]/dvec[2]
    dnds <- dvec[3]/dvec[4]
    ni <- pnps/dnds
    ft <- fisher.test( dmat)
    ftdf[ j, 3:13] <- c( dvec, pnps, dnds, ni, ft$estimate, ft$conf.int, ft$p.value)
  }
  ftdf
}

ftls <- lapply( genes, function(g) tryCatch(expr = afunc(g),
           error = function(e) NA))
ftls <- do.call( rbind, ftls)
ftls <- ftls[ !apply( ftls, 1, anyNA), ]
df <- as.data.frame(ftls)
# adjust p-values for multiple testing
df$adjp <- p.adjust( df$p, method = 'holm')

write.csv( x = df, file = fout, row.names = F, quote = F)

spdfs <- split( df, df$species)
