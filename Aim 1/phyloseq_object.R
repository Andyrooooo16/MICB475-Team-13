library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(picante)

#### Load data ####
# Change file paths as necessary
metafp <- "../halfvarson_metadata_wrangled.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "../ibd_export/table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "../ibd_export/tax_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "../ibd_export/tree_export/tree.nwk"
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


#### Adding column to metadata file that shows different cases of disease severity and inflammation ####
meta_final <- meta %>%
      mutate(disease_and_inflammation = case_when(
          disease_severity == "Low" & inflammation == FALSE ~ "low_not_inflammed",
          disease_severity == "Low" & inflammation == TRUE ~ "low_inflammed",
          disease_severity == "Medium" & inflammation == FALSE ~ "medium_not_inflammed",
          disease_severity == "Medium" & inflammation == TRUE ~ "low_inflammed",
          disease_severity == "High" & inflammation == FALSE ~ "low_not_inflammed",
          disease_severity == "High" & inflammation == TRUE ~ "low_inflammed",
      ))


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

######### ANALYZE ##########
# Remove non-bacterial sequences, if any
ibd_filt <- subset_taxa(ibd,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
ibd_filt <- prune_samples(sample_sums(ibd_filt)>105, ibd_filt)

### INDICATOR SPECIES 

library(indicspecies)

#### Indicator Species/Taxa Analysis ####
# glom to Genus
ibd_genus <- tax_glom(ibd_filt, "Genus", NArm = TRUE)
ibd_genus_RA <- transform_sample_counts(ibd_genus, fun=function(x) x/sum(x))

#ISA
isa_ibd <- multipatt(t(otu_table(ibd_genus_RA)), cluster = sample_data(ibd_genus_RA)$`inflammation_with_surgery`)
summary(isa_ibd)
taxtable <- tax_table(ibd) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of 
# anything beyond the glomed taxa level
isa_ibd$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()

View(sample_data(ibd_genus_RA))

library(ggpubr)

#################
ibd_rare <- rarefy_even_depth(ibd_filt, rngseed = 1, sample.size = 124392)

#### Alpha diversity ######
plot_richness(ibd_rare) 

plot_richness(ibd_rare, measures = c("Shannon","Chao1")) 

gg_richness_inf_surg <- plot_richness(ibd_rare, x = "inflammation_with_surgery", measures = c("Shannon","Chao1")) +
  xlab("Inflammation with Surgery") +
  geom_boxplot()
gg_richness_inf_surg

gg_richness_dis_sev <- plot_richness(ibd_rare, x = "disease_severity", measures = c("Shannon","Chao1")) +
  xlab("Disease Severity") +
  geom_boxplot()
gg_richness_dis_sev

gg_richness_cd_loc <- plot_richness(ibd_rare, x = "cd_location", measures = c("Shannon","Chao1")) +
  xlab("Crohn's Disease Location") +
  geom_boxplot()
gg_richness_cd_loc

# phylogenetic diversity

# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(ibd_rare)), phy_tree(ibd_rare),
                 include.root=F) 
?pd

# add PD to metadata table
sample_data(ibd_rare)$PD <- phylo_dist$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(ibd_rare), aes(inflammation_with_surgery, PD)) + 
  geom_boxplot() +
  xlab("Subject ID") +
  ylab("Phylogenetic Diversity")

# view plot
plot.pd


#### Beta diversity #####
bc_dm <- distance(ibd_rare, method="bray")
# check which methods you can specify
?distance

pcoa_bc <- ordinate(ibd_rare, method="PCoA", distance=bc_dm)

# Calprotectin + Disease Severity
plot_ordination(ibd_rare, pcoa_bc, color = "calprotectin", shape="disease_severity")

gg_pcoa <- plot_ordination(ibd_rare, pcoa_bc, color = "calprotectin", shape="disease_severity")+
  scale_color_gradient(name = "Calprotectin Levels", low = "orange", high = "blue") +
  labs(pch="Disease Severity", col = "Calprotectin levels")
gg_pcoa

#### Taxonomy bar plots ####

# Plot bar plot of taxonomy
plot_bar(ibd_rare, fill="Phylum") 

# Convert to relative abundance
ibd_RA <- transform_sample_counts(ibd_rare, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first
ibd_phylum <- tax_glom(ibd_RA, taxrank = "Phylum", NArm=FALSE)

gg_taxa <- plot_bar(ibd_phylum, fill="Phylum") + 
  facet_wrap(.~inflammation_with_surgery, scales = "free_x")
gg_taxa


