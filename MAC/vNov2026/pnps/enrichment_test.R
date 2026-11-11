ffn <- 'gene_cat_class.tsv'
fni <- 'gene_sp_ni_pairwise.csv'
fout <- 'enrichment_test_out.csv'

# minimum no. of entries in a category
mindata <- 15

# maximum dS for a comparison
maxds <- 0.8

dfn <- read.csv(file = ffn, sep = '\t')
dni <- read.csv(file = fni)

# apply max dS threshold
dni <- dni[ dni$dS <= maxds, ]

# species with valid NI 
sps <- unique( c(dni$focal_species,dni$closest_species))
cat( sprintf("Species with valid NI = %s\n", paste(sps,collapse = ', ')))

# list of NI for each gene
ls_gene_ni <- split( dni, dni$gene)

# total no. of genes
ngenes <- length(ls_gene_ni)
cat(sprintf("# genes = %i\n",ngenes))


# genes with all NI < 1
posg <- names(ls_gene_ni)[sapply( ls_gene_ni, function(X) all(X$NI<1))]

## limit data to genes with categories
dfn <- dfn[ dfn$category != "uncategorized", ]

dfn <- dfn[ dfn$gene %in% dni$gene, ]

cat_genes <- split( dfn$gene, dfn$category)

cat_ni <- lapply( cat_genes, function(x) 
  do.call( rbind, ls_gene_ni[x]))
cat_ni <- cat_ni[ sapply( cat_ni, nrow) >= mindata]
ncat <- length(cat_ni)
cat(sprintf("# Gene categories with at least %i comparisons = %i\n", mindata, ncat))

# no. of values in each func category
cat_size <- sapply( cat_ni, nrow)

# total no. of values
tot_size <- sum(cat_size)

# no. of values that are positive
tot_pos <- sum( sapply( cat_ni, function(X) sum(X$NI<1)))

out <- matrix( NA, ncat, 2)
rownames(out) <- names(cat_ni)
colnames(out) <- c('pos','neg')

out <- as.data.frame(out)

for(i in 1:ncat){
  cdf <- cat_ni[[i]]
  # NI < 1
  npos <- sum(cdf$NI < 1)
  out[ i, 1] <- npos
  # NI > 1
  out[i,2] <- cat_size[i] - npos
}

# positive and negative in others
out$opos <- tot_pos - out$pos
out$oneg <- (tot_size - tot_pos) - out$neg

# odds ratio
out$OR <- (out$pos/out$opos)/(out$neg/out$oneg)

# Fisher's exact test
out$pval <- apply( out[,1:4], 1, function(x) 
  fisher.test(x = matrix( x, 2, 2, byrow = T), alternative = 'g')$p.value)
