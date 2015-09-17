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
library("phyloseq")
library("biom")

# reading json/biom file
#greengenes = read_biom("~/Development/data/matrix.greengenes.biom")
greengenes = read_biom("~/Development/data/matrix.test.biom")

# columns
observation_metadata(greengenes)
# there are 5 levels
colnames(observation_metadata(greengenes))

sample_metadata(greengenes)

# Create pyloseq object
otumatrix = as(biom_data(greengenes), "matrix")
OTU = otu_table(otumatrix, taxa_are_rows=TRUE)
taxmat = as.matrix(observation_metadata(x), rownames.force=TRUE)
TAX = tax_table(taxmat)
META <- sample_metadata(greengenes)

physeq = phyloseq(OTU, TAX, META)


# Basic data access
taxa_names(physeq)

# Alpha Diversity
plot_richness(physeq) + geom_boxplot()


#
plot_ordination(physeq , ordinate(physeq) ) + geom_point(size=10)
# plot_ordination(physeq , ordinate(physeq, "MDS") ) + geom_point(size=10)

plot_heatmap(physeq , title = "Super Coooool")



 