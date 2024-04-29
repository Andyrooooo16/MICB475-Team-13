### Loading Packages ###

library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#### Load data ####
# Change file paths as necessary
metafp <- "hc_halfvarson_metadata_wrangled_new.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "HC_ibd_exports/table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "HC_ibd_exports/tax_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "HC_ibd_exports/tree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

meta <- meta %>%
  mutate(DiseaseSev_with_inflammation = case_when(
    disease_severity == "High" & inflammation == TRUE ~ "High with inflammation",
    disease_severity == "High" & inflammation == FALSE ~ "High no inflammation",
    disease_severity == "Medium" & inflammation == TRUE ~ "Medium with inflammation",
    disease_severity == "Medium" & inflammation == FALSE ~ "Medium no inflammation",
    disease_severity == "Low" & inflammation == TRUE ~ "Low with inflammation",
    disease_severity == "Low" & inflammation == FALSE ~ "Low no inflammation",
    disease_severity == "Healthy Control" & inflammation == FALSE ~ "Healthy Control",
    disease_severity == "Healthy Control" & inflammation == TRUE ~ "Healthy Control",
    
  )
  )
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
save(ibd, file = "ibd.RData")

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(ibd)
sample_data(ibd)
tax_table(ibd)
phy_tree(ibd)

#### Filtering ####
# Remove non-bacterial sequences, if any
ibd_filt <- subset_taxa(ibd,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
# Remove samples with less than 100 reads
ibd_filt <- prune_samples(sample_sums(ibd_filt)>105, ibd_filt)

save(ibd_filt, file = "ibd_filt.RData")

#### Rarefied Data ####
ibd_rare <- rarefy_even_depth(ibd_filt, rngseed = 1, sample.size = 124392)
save(ibd_rare, file = "ibd_rare.RData")

