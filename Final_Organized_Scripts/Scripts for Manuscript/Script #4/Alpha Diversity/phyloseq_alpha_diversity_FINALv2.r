### Loading Packages ###

library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(picante)
library(ggplot2)
library(ggpubr)
library(rstatix) 
library(patchwork)

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
                                                   levels = c("Healthy Control", "no_inflammation_no_surgery", "no_inflammation_with_surgery", 
                                                              "inflammation_no_surgery", "inflammation_with_surgery"), ordered = TRUE)

diversity_long_infsurg <- pivot_longer(diversity_data, c("Shannon", "Chao1"), 
                                       names_to = "Measure", values_to = "Value")

shannon_data_iws <- diversity_long_infsurg %>%
  filter(Measure == "Shannon")

chao1_data_iws <- diversity_long_infsurg %>%
  filter(Measure == "Chao1")

compare_means_shannon <- compare_means(Shannon ~ inflammation_with_surgery,  data = diversity_data)
compare_means_chao1 <- compare_means(Chao1 ~ inflammation_with_surgery,  data = diversity_data)

significant_comparisons_shannon <-list(c("Healthy Control", "no_inflammation_with_surgery"), 
                                       c("no_inflammation_no_surgery", "no_inflammation_with_surgery"))

p_Shannon_iws <- ggplot(shannon_data_iws, aes(x = inflammation_with_surgery, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#56B4E9")) +
  labs(x = "Inflammation with Surgery",
       y = "") + stat_compare_means(method = "wilcox.test", comparisons = significant_comparisons_shannon, label = "p.signif")

significant_comparisons_chao1 <-list(c("Healthy Control", "no_inflammation_with_surgery"), 
                                     c("no_inflammation_no_surgery", "no_inflammation_with_surgery"), 
                                     c("Healthy Control", "inflammation_no_surgery"), 
                                     c("Healthy Control", "inflammation_with_surgery"), 
                                     c("no_inflammation_with_surgery", "inflammation_no_surgery"))

p_Chao1_iws <- ggplot(chao1_data_iws, aes(x = inflammation_with_surgery, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Chao1" = "#E69F00")) +
  labs(x = "Inflammation with Surgery",
       y = "Diversity Index") + stat_compare_means(method = "wilcox.test", comparisons = significant_comparisons_chao1, label = "p.signif")

p_iws <-  p_Chao1_iws | p_Shannon_iws
p_iws
ggsave("inflam_w_surg.png", p_iws)

#

#
#
#
# Disease Severity Alpha diversity

diversity_data$disease_severity <- factor(sample_data(ibd_rare)$disease_severity, 
                                          levels = c("Healthy Control", "Low", "Medium", "High"), ordered = TRUE)

diversity_long_severity <- pivot_longer(diversity_data, c("Shannon", "Chao1"), 
                                        names_to = "Measure", values_to = "Value")

diversity_long_severity$disease_severity[is.na(diversity_long_severity$disease_severity)] <- 'Healthy Control'

diversity_long_severity$disease_severity <- factor(diversity_long_severity$disease_severity, 
                                                   levels = c("Healthy Control","Low", "Medium", "High"), ordered = TRUE)

shannon_data_ds <- diversity_long_severity %>%
  filter(Measure == "Shannon")

chao1_data_ds <- diversity_long_severity %>%
  filter(Measure == "Chao1")

compare_means_shannon <- compare_means(Shannon ~ disease_severity,  data = diversity_data)
compare_means_chao1 <- compare_means(Chao1 ~ disease_severity,  data = diversity_data)

significant_comparisons_shannon_ds <-list(c("Healthy Control", "Medium"), 
                                          c("Healthy Control", "High"), 
                                          c("Low", "Medium"),
                                          c("Low", "High"))

significant_comparisons_chao1_ds <-list(c("Healthy Control", "Low"),
                                        c("Healthy Control", "Medium"), 
                                        c("Healthy Control", "High"), 
                                        c("Low", "Medium"),
                                        c("Low", "High"))

p_Shannon_ds <- ggplot(shannon_data_ds, aes(x = disease_severity, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#56B4E9", "Chao1" = "#E69F00")) +
  labs(x = "Disease Severity",
       y = "")  + stat_compare_means(method = "wilcox.test", comparisons = significant_comparisons_shannon_ds, label = "p.signif")

p_Chao1_ds <- ggplot(chao1_data_ds, aes(x = disease_severity, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#56B4E9", "Chao1" = "#E69F00")) +
  labs(x = "Disease Severity",
       y = "")  + stat_compare_means(method = "wilcox.test", comparisons = significant_comparisons_chao1_ds, label = "p.signif")

p_ds <-  p_Chao1_ds | p_Shannon_ds
p_ds
ggsave("dseverity.png", p_ds)

#
#
#
# Disease Severity with Inflammation Alpha Diversity
diversity_data$inflammation <- factor(sample_data(ibd_rare)$inflammation) 

diversity_data$DiseaseSev_with_inflammation <- sample_data(ibd_rare)$DiseaseSev_with_inflammation

diversity_long_ds_inf <- pivot_longer(diversity_data, c("Shannon", "Chao1"), 
                                      names_to = "Measure", values_to = "Value")

true_inflammation <- filter(diversity_long_ds_inf, inflammation == TRUE)
false_inflammation <- filter(diversity_long_ds_inf, inflammation == FALSE)

data_Shannon_ds_inf_true <- subset(true_inflammation, Measure == "Shannon")
data_Chao1_ds_inf_true <- subset(true_inflammation, Measure == "Chao1")

data_Shannon_ds_inf_false <- subset(false_inflammation, Measure == "Shannon")
data_Chao1_ds_inf_false <- subset(false_inflammation, Measure == "Chao1")

data_Shannon_ds_inf_true$DiseaseSev_with_inflammation <- factor(data_Shannon_ds_inf_true$DiseaseSev_with_inflammation, 
                                                                levels = c("Healthy Control","Low with inflammation", "Medium with inflammation", "High with inflammation"), ordered = TRUE)

data_Shannon_ds_inf_false$DiseaseSev_with_inflammation <- factor(data_Shannon_ds_inf_false$DiseaseSev_with_inflammation, 
                                                                 levels = c("Healthy Control","Low no inflammation", "Medium no inflammation", "High no inflammation"), ordered = TRUE)

data_Chao1_ds_inf_true$DiseaseSev_with_inflammation <- factor(data_Chao1_ds_inf_true$DiseaseSev_with_inflammation, 
                                                              levels = c("Healthy Control","Low with inflammation", "Medium with inflammation", "High with inflammation"), ordered = TRUE)

data_Chao1_ds_inf_false$DiseaseSev_with_inflammation <- factor(data_Chao1_ds_inf_false$DiseaseSev_with_inflammation, 
                                                               levels = c("Healthy Control","Low no inflammation", "Medium no inflammation", "High no inflammation"), ordered = TRUE)

compare_means_shannon <- compare_means(Value ~ DiseaseSev_with_inflammation,  data = data_Shannon_ds_inf_true)
compare_means_chao1 <- compare_means(Value ~ DiseaseSev_with_inflammation,  data = data_Chao1_ds_inf_true)

compare_means_shannon <- compare_means(Value ~ DiseaseSev_with_inflammation,  data = data_Shannon_ds_inf_false)
compare_means_chao1 <- compare_means(Value ~ DiseaseSev_with_inflammation,  data = data_Chao1_ds_inf_false)

#true inflammation for Shannon and Chao were ns so 
significant_comparisons_chao1_ds_inf_true <-list(c("Medium with inflammation","Healthy Control"),
                                                 c("Medium with inflammation","Low with inflammation"))

significant_comparisons_shannon_ds_inf_true <-list(c("Medium with inflammation","Healthy Control"),
                                                   c("Medium with inflammation","Low with inflammation"))

significant_comparisons_chao1_ds_inf_false <-list(c("Medium no inflammation","Healthy Control"),
                                                  c("Medium no inflammation","Low no inflammation"))

significant_comparisons_shannon_ds_inf_false <-list(c("Medium no inflammation","Healthy Control"),
                                                    c("Medium no inflammation","Low no inflammation"))

compare_means_shannon <- compare_means(Value ~ DiseaseSev_with_inflammation,  data = data_Shannon_ds_inf_false)
compare_means_chao1 <- compare_means(Value ~ DiseaseSev_with_inflammation,  data = data_Chao1_ds_inf_false)

p_Shannon_ds_inf_true <- ggplot(data_Shannon_ds_inf_true, aes(x = DiseaseSev_with_inflammation, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#AED6F1")) +
  labs(x = "Disease Severity and Inflammation", y = "Diversity Index") +
  facet_wrap(~Measure, scales = "free_y", ncol = 2) +
  coord_cartesian(ylim = c(0, 5)) + stat_compare_means(method = "wilcox.test", comparisons = significant_comparisons_shannon_ds_inf_true, label = "p.signif")
# ns but can remove commment + stat_compare_means(method = "wilcox.test", comparisons = significant_comparisons_shannon_ds_inf_true, label = "p.signif")

p_Chao1_ds_inf_true <- ggplot(data_Chao1_ds_inf_true, aes(x = DiseaseSev_with_inflammation, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Chao1" = "#F9E79F")) +
  labs(x = "", y = "Diversity Index") +
  facet_wrap(~Measure, scales = "free_y", ncol = 2) +
  coord_cartesian(ylim = c(0, 500)) + stat_compare_means(method = "wilcox.test", comparisons = significant_comparisons_chao1_ds_inf_true, label = "p.signif")

p_dsinf_true <-  p_Chao1_ds_inf_true | p_Shannon_ds_inf_true

false_inflammation$DiseaseSev_with_inflammation <- factor(false_inflammation$DiseaseSev_with_inflammation, 
                                                          levels = c("Healthy Control","Low no inflammation", "Medium no inflammation", "High no inflammation"), ordered = TRUE)

p_shannon_dsinf_false <- ggplot(data_Shannon_ds_inf_false, aes(x = DiseaseSev_with_inflammation, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Shannon" = "#56B4E9")) +
  labs(x = "Disease Severity and Inflammation",
       y = "") +
  facet_wrap(~Measure, scales = "free_y", ncol = 4) + stat_compare_means(method = "wilcox.test", comparisons = significant_comparisons_shannon_ds_inf_false, label = "p.signif")

p_chao1_dsinf_false <- ggplot(data_Chao1_ds_inf_false, aes(x = DiseaseSev_with_inflammation, y = Value, fill = Measure)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Chao1" = "#E69F00")) +
  labs(x = "Disease Severity and Inflammation",
       y = "Diversity Index") +
  facet_wrap(~Measure, scales = "free_y", ncol = 4) + stat_compare_means(method = "wilcox.test", comparisons = significant_comparisons_chao1_ds_inf_false, label = "p.signif")

p_dsinf_false <-  p_chao1_dsinf_false | p_shannon_dsinf_false


p_dsinf_true 
ggsave("dsinf_true.png", p_dsinf_true)
p_dsinf_false
ggsave("dsinf_false.png", p_dsinf_false)
