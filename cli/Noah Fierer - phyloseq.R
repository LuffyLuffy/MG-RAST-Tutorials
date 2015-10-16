source("http://bioconductor.org/biocLite.R")
biocLite("phyloseq")

install.packages("ape")
install.packages("ade4")
install.packages("doParallel")
install.packages("foreach")
install.packages("ggplot2")
install.packages("igraph0")
install.packages("picante")
install.packages("plyr")
install.packages("RJSONIO")
install.packages("scales")
install.packages("testthat")
install.packages("vegan")

library(ggplot2)
library("biom")
library("phyloseq")

# reading json/biom file
amplicon.biom = read_biom("/Users/Andi/Development/tmp/data/phyloseq/amplicon.greengenes.biom")
wgs.biom      = read_biom("/Users/Andi/Development/tmp/data/phyloseq/wgs.greengenes.biom")
merged.biom      = read_biom("/Users/Andi/Development/tmp/data/phyloseq/merged.greengenes.biom")
#refseq    = read_biom("~/Development/tmp/data/phyloseq/test4.biom")

### Create amplicon object
# OTU matrix
otumatrix = as(biom_data(amplicon.biom), "matrix")
OTU = otu_table(otumatrix, taxa_are_rows=TRUE)

# import metadata
metadata_file="/Users/Andi/Development/tmp/data/phyloseq/amplicon.metadata"
META = import_qiime_sample_data(metadata_file)

# import taxa
taxa_file="/Users/Andi/Development/tmp/data/phyloseq/amplicon.taxa"
taxa.frame <- read.csv(taxa_file , header = TRUE , sep = "\t")

taxmat = as.matrix(taxa.frame)
rownames(taxmat) <- taxmat[,1]
taxmat <- taxmat[,-1]
#colnames(taxmat2) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", 
#                      "Species")
TAX = tax_table(taxmat)

amplicon = phyloseq(OTU, TAX , META)


### Create wgs object
# OTU matrix
otumatrix = as(biom_data(wgs.biom), "matrix")
OTU = otu_table(otumatrix, taxa_are_rows=TRUE)

# import metadata
metadata_file="/Users/Andi/Development/tmp/data/phyloseq/wgs.metadata"
META = import_qiime_sample_data(metadata_file)

# import taxa
taxa_file="/Users/Andi/Development/tmp/data/phyloseq/wgs.taxa"
taxa.frame <- read.csv(taxa_file , header = TRUE , sep = "\t")

taxmat = as.matrix(taxa.frame)
rownames(taxmat) <- taxmat[,1]
taxmat <- taxmat[,-1]
#colnames(taxmat2) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", 
#                      "Species")
TAX = tax_table(taxmat)

wgs = phyloseq(OTU, TAX , META)


### Create merged 16S object
# OTU matrix
otumatrix = as(biom_data(merged.biom), "matrix")
OTU = otu_table(otumatrix, taxa_are_rows=TRUE)

# import metadata
metadata_file="/Users/Andi/Development/tmp/data/phyloseq/merged.metadata"
META = import_qiime_sample_data(metadata_file)

# import taxa
taxa_file="/Users/Andi/Development/tmp/data/phyloseq/merged.taxa"
taxa.frame <- read.csv(taxa_file , header = TRUE , sep = "\t")

taxmat = as.matrix(taxa.frame)
rownames(taxmat) <- taxmat[,1]
taxmat <- taxmat[,-1]
#colnames(taxmat2) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", 
#                      "Species")
TAX = tax_table(taxmat)

merged = phyloseq(OTU, TAX , META)







# Plots Alpha Diversity Amplicon
plot_richness(amplicon) + geom_boxplot()
plot_ordination(amplicon , ordinate(amplicon) ) + geom_point(size=10)
plot_ordination(amplicon , ordinate(amplicon, "MDS") ) + geom_point(size=10)
plot_ordination(amplicon , ordinate(amplicon, "PCoA") , color = "sample.data.biome") + geom_point(size=5)
plot_heatmap(amplicon , title = "Cross-Site" , taxa.label = "Genus")
plot_bar(amplicon, fill = "Family")


# Plots Alpha Diversity WGS
plot_richness(wgs) + geom_boxplot()
plot_ordination(wgs , ordinate(wgs) ) + geom_point(size=10)
plot_ordination(wgs , ordinate(wgs, "MDS") ) + geom_point(size=10)
plot_ordination(wgs , ordinate(wgs, "PCoA") , color = "sample.data.biome") + geom_point(size=5)
plot_heatmap(wgs , title = "Cross-Site (WGS)" , taxa.label = "Genus")
plot_bar(wgs, fill = "Family")


# Plots Alpha Diversity Merged
plot_richness(merged) + geom_boxplot()
plot_ordination(merged , ordinate(merged) ) + geom_point(size=10)
plot_ordination(merged , ordinate(merged, "MDS") ) + geom_point(size=10)
plot_ordination(merged , ordinate(merged, "PCoA") , color = "sample.data.biome") + geom_point(size=5)
plot_heatmap(merged , title = "Cross-Site (WGS)" , taxa.label = "Genus")
plot_bar(merged, fill = "Family")

plot_ordination(wgs , ordinate(wgs, "PCoA") , color = "library.type") + geom_point(size=5)
plot_ordination(amplicon , ordinate(amplicon, "PCoA") , color = "library.type") + geom_point(size=5)
plot_ordination(merged , ordinate(merged, "PCoA") , color = "library.type") + geom_point(size=5)
