setwd('~/ASU Dropbox/Ravina Telkar/intersp/final_315_samples/mito')

f1 <- 'gene_species_anc_DnDs_NI.csv'
f2 <- 'gene_species_anc_DnDs_NI_filtered_0.6.csv'

df1 <- read.csv(file = f1)
df2 <- read.csv(file = f2)

df1$name <- apply( df1[ , c('gene','species')], 1, paste, collapse='_')
df2$name <- apply( df2[ , c('gene','species')], 1, paste, collapse='_')

comnames <- intersect( df1$name, df2$name)

df1 <- df1[ match( comnames, df1$name), ]
df2 <- df2[ match( comnames, df2$name), ]

plot( df1$obs.NI, df2$obs.NI, log = 'xy')
abline(0,1)
abline(h=1)
abline(v=1)

wilcox.test( df1$obs.NI, df2$obs.NI)

plot( df1$p, df2$p)
abline(0,1)
abline(h=0.05)
abline(v=0.05)