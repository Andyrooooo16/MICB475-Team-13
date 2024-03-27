### Loading Packages ###

library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(picante)
library(ggplot2)
library(ggpubr)
library(rstatix)

custom_theme <- theme_minimal() + 
  theme(legend.position = "none", 
        text = element_text(family = "Helvetica"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.background = element_rect(fill = "white", color = NA), 
        panel.border = element_rect(color = "black", fill = NA, size = 0.5), 
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 14, face = "bold"))

# Set the custom theme as the default theme for all ggplot2 plots
theme_set(custom_theme)


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

meta <- meta %>%
  mutate(DiseaseSev_with_inflammation = case_when(
    disease_severity == "High" & inflammation == TRUE ~ "High with inflammation",
    disease_severity == "High" & inflammation == FALSE ~ "High no inflammation",
    disease_severity == "Medium" & inflammation == TRUE ~ "Medium with inflammation",
    disease_severity == "Medium" & inflammation == FALSE ~ "Medium no inflammation",
    disease_severity == "Low" & inflammation == TRUE ~ "Low with inflammation",
    disease_severity == "Low" & inflammation == FALSE ~ "Low no inflammation",
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
# Remove samples with less than 100 reads
ibd_filt <- prune_samples(sample_sums(ibd_filt)>105, ibd_filt)

save(ibd_filt, file = "ibd_filt.RData")

#################
ibd_rare <- rarefy_even_depth(ibd_filt, rngseed = 1, sample.size = 124392)

#### Alpha diversity ######
plot_richness(ibd_rare) 

plot_richness(ibd_rare, measures = c("Shannon","Chao1")) 

#
#
#
# Calculate diversity indices
diversity_data <- estimate_richness(ibd_rare, measures = c("Shannon", "Chao1"))


# Inflammation with Surgery Alpha Diversity
diversity_data$inflammation_with_surgery <- sample_data(ibd_rare)$inflammation_with_surgery

diversity_data$inflammation_with_surgery <- factor(diversity_data$inflammation_with_surgery, 
                                                   levels = c("no_inflammation_no_surgery", "no_inflammation_with_surgery", 
                                                              "inflammation_no_surgery", "inflammation_with_surgery"), ordered = TRUE)

diversity_long_infsurg <- pivot_longer(diversity_data, c("Shannon", "Chao1"), 
                                       names_to = "Measure", values_to = "Value")

p_iws <- ggplot(diversity_long_infsurg, aes(x = inflammation_with_surgery, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#56B4E9", "Chao1" = "#E69F00")) +
  labs(x = "Inflammation with Surgery",
       y = "Diversity Index") +
  facet_wrap(~Measure, scales = "free_y", ncol = 2)

p_iws +
  geom_pwc(
    aes(group = inflammation_with_surgery), tip.length = 0,
    method = "t_test", label = "p.adj.format",
    bracket.nudge.y = 0.1
  )

#
#
#
#
# Disease Severity Alpha diversity

diversity_data$disease_severity <- factor(sample_data(ibd_rare)$disease_severity, 
                                          levels = c("Low", "Medium", "High"), ordered = TRUE)

diversity_long_severity <- pivot_longer(diversity_data, c("Shannon", "Chao1"), 
                                        names_to = "Measure", values_to = "Value")

diversity_long_severity$disease_severity <- factor(diversity_long_severity$disease_severity, 
                                                   levels = c("Low", "Medium", "High"), ordered = TRUE)




p_ds <- ggplot(diversity_long_severity, aes(x = disease_severity, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#56B4E9", "Chao1" = "#E69F00")) +
  labs(x = "Disease Severity",
       y = "Diversity Index") +
  facet_wrap(~Measure, scales = "free_y", ncol = 2)

p_ds +
  geom_pwc(
    aes(group = disease_severity), tip.length = 0,
    method = "t_test", label = "p.adj.format",
    bracket.nudge.y = 0.1
  )

#
#
#
# CD Location alpha diversity

diversity_data$cd_location <- sample_data(ibd_rare)$cd_location

filtered_diversity_data <- diversity_data %>%
  filter(cd_location != "Ileocolonic and Upper-GI (L3+L4)")

diversity_long_location_filtered <- pivot_longer(filtered_diversity_data, c("Shannon", "Chao1"), names_to = "Measure", values_to = "Value")

p_cdloc <- ggplot(diversity_long_location_filtered, aes(x = cd_location, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#56B4E9", "Chao1" = "#E69F00")) +
  labs(x = "CD Location",
       y = "Diversity Index") +
  facet_wrap(~Measure, scales = "free_y", ncol = 2) 

p_cdloc + 
  geom_pwc(
    aes(group = cd_location), tip.length = 0,
    method = "t_test", label = "p.adj.format",
    bracket.nudge.y = 0.1
  )

#
#
#
# Disease Severity with Inflammation Alpha Diversity
diversity_data$DiseaseSev_with_inflammation <- sample_data(ibd_rare)$DiseaseSev_with_inflammation

diversity_long_ds_inf <- pivot_longer(diversity_data, c("Shannon", "Chao1"), 
                                      names_to = "Measure", values_to = "Value")

p_dsinf <- ggplot(diversity_long_ds_inf, aes(x = DiseaseSev_with_inflammation, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#56B4E9", "Chao1" = "#E69F00")) +
  labs(x = "Disease Severity and Inflammation",
       y = "Diversity Index") +
  facet_wrap(~Measure, scales = "free_y", ncol = 2)

p_dsinf +
  geom_pwc(
    aes(group = DiseaseSev_with_inflammation), tip.length = 0,
    method = "t_test", label = "p.adj.format",
    bracket.nudge.y = 0.1
  )


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

ibd_rare$inflammation <- as.numeric(as.character(ibd_rare$inflammation))

bc_dm <- distance(ibd_rare, method="bray")
# check which methods you can specify
?distance

pcoa_bc <- ordinate(ibd_rare, method="PCoA", distance=bc_dm)

# Calprotectin + Disease Severity
plot_ordination(ibd_rare, pcoa_bc, color = "inflammation", shape="disease_severity")

gg_pcoa <- plot_ordination(ibd_rare, pcoa_bc, color = "inflammation", shape="disease_severity")+
  scale_color_gradient(name = "inflammation", low = "red", high = "blue") +
  labs(pch="Disease Severity", col = "Calprotectin levels")
gg_pcoa

#### Taxonomy bar plots ####

# Plot bar plot of taxonomy
plot_bar(ibd_rare, fill="Phylum") 

# Convert to relative abundance
ibd_RA <- transform_sample_counts(ibd_rare, function(x) x/sum(x))

ggplot(data = melted_df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Abundance of OTUs across Samples", x = "Sample", y = "Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# To remove black bars, "glom" by phylum first
ibd_phylum <- tax_glom(ibd_RA, taxrank = "Phylum", NArm=FALSE)

gg_taxa <- plot_bar(ibd_phylum, fill="Phylum") + 
  facet_wrap(.~inflammation_with_surgery, scales = "free_x")
gg_taxa

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
isa_ibd_table <- isa_ibd$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable, by="ASV") %>%
  filter(p.value < 0.05)

# Write the table to a CSV file
write.csv(isa_ibd_table, "isa_ibd_table.csv", row.names = FALSE)



