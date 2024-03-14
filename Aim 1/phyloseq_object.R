library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(picante)

#### Load data ####
# Change file paths as necessary
metafp <- "halfvarson_metadata_wrangled.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "/Users/jamesforward/Desktop/MICB 475/475 project/IBD project/ibd_export/table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "/Users/jamesforward/Desktop/MICB 475/475 project/IBD project/ibd_export/tax_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "/Users/jamesforward/Desktop/MICB 475/475 project/IBD project/ibd_export/tree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# OTU tables should be a matrix
# with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$'sample-id'
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
ibd <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(ibd)
sample_data(ibd)
tax_table(ibd)
phy_tree(ibd)

rarecurve(t(as.data.frame(otu_table(ibd))), cex=0.1)
ibd_rare <- rarefy_even_depth(ibd, rngseed = 1, sample.size = 124392)

#### Beta diversity #####
unifrac <- UniFrac(ibd_rare, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)

# check which methods you can specify
?distance

pcoa_unifrac <- ordinate(ibd_rare, method="PCoA", distance=unifrac)

plot_ordination(ibd_rare, pcoa_unifrac, color = "days_post_transplant", shape="cage_id")

gg_pcoa <- plot_ordination(ibd_rare, pcoa_unifrac, color = "days_post_transplant", shape="cage_id") +
  scale_color_gradient(name = "Days after transplant", low = "yellow", high = "blue") +
  labs(pch="Cage ID", col = "Days after transplant")
gg_pcoa

ggsave("plot_pcoa.png"
       , gg_pcoa
       , height=4, width=5)

