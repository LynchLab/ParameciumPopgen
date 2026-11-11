# infer genes presence/absence based on depth of coverage for each gene region

# command line arguments
args <- commandArgs(trailingOnly = T)
dr = args[1] # sample dir with cov.bed

# paths
fin <- file.path( dr, 'cov.bed')
fgp <- file.path( dr, 'genes_present.txt')
fmg <- file.path( dr, 'missing_genes.txt')
fex <- file.path( dr, 'excess_depth_genes.txt')

# minimum coverage threshold,
mbt = 0.6 # fraction of gene length overlapping at least one read

df <- read.csv( file = fin, header = F, sep = '\t', stringsAsFactors = F)

# ID
df$id = df$V4

# gene length
df$glen = df$V3-df$V2

# breadth of coverage
df$bcov = df$V8/df$glen

# per-base depth of coverage
df$dcov = round(df$V7/df$glen)

# genes with zero breadth of coverage
gcb0 <- df[ df$bcov == 0, ]

# genes with min bcov
gpb <- df[ df$bcov >= mbt, ]

# genes with zero depth of coverage
gcd0 <- df[ df$dcov == 0, ]

# find the depth cut-off to subset dcd
# first extract non-zero values
dp <- df$dcov[ df$dcov > 0]
# find their mean
md <- mean(dp)
# further, keep only values smaller than twice of the above mean
dp <- dp[ dp < 2*md]
# recalculate mean
md <- mean(dp) 
# get standard deviation
std <- sd(dp)
# minimum depth threshold, for a gaussian distribution 3SD covers 99.7% data
mdt = md - 3*std
# extract genes satsifying above threshold
gpd <- df[ df$dcov >= mdt & df$dcov > 0, ]
# genes with excess depth of coverage
E <- df[ df$dcov >= 2*md, ]

# take intersection of genes selected for breadth & depth of coverage
gp <- df[ intersect( rownames(gpb), rownames(gpd)), ]
# take union of genes with either zero breadth or zero depth
g0 <- df[ union( rownames(gcb0), rownames(gcd0)), ]

# write above vector of genes to files
write.table( x = gp, file = fgp, sep='\t', row.names = F, quote = F)
write.table( x = g0, file = fmg, sep='\t', row.names = F, quote = F)
write.table( x = E, file = fex, sep='\t', row.names = F, quote = F)
