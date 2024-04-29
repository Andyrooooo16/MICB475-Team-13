### Loading Packages ###

library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(picante)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(patchwork)
library(indicspecies)



### Setting theme for plots ###
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
load("ibd_filt.RData")
load("ibd_rare.RData")
load("ibd.RData")

#### Taxonomy bar plots ####

# Plot bar plot of taxonomy
plot_bar(ibd_rare, fill="Phylum") 

# Convert to relative abundance
ibd_RA <- transform_sample_counts(ibd_rare, function(x) x/sum(x))

# ggplot(data = melted_df, aes(x = Sample, y = Abundance, fill = Phylum)) +
#   geom_bar(stat = "identity", position = "stack") +
#   theme_minimal() +
#   labs(title = "Abundance of OTUs across Samples", x = "Sample", y = "Abundance") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# To remove black bars, "glom" by phylum first
ibd_phylum <- tax_glom(ibd_RA, taxrank = "Phylum", NArm=FALSE)

gg_taxa <- plot_bar(ibd_phylum, fill="Phylum") + 
  facet_wrap(.~inflammation_with_surgery, scales = "free_x")
gg_taxa

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

isa_ibd_table_filt <- isa_ibd_table %>%
  filter(stat > 0.87)

ibd_RA_melt <- psmelt(ibd_RA)

summary_by_otu_mean <- ibd_RA_melt %>%
  group_by(OTU, inflammation_with_surgery) %>%
  summarize(
    mean_abundance = mean(Abundance, na.rm = TRUE)
  )

summary_by_otu_mean <- summary_by_otu_mean %>%
  rename(ASV = OTU)

filtered_unique_asv <- summary_by_otu_mean %>%
  filter(ASV %in% unique(isa_ibd_table_filt$ASV))

isa_ibd_table_filt_asv_genus <- isa_ibd_table_filt %>%
  select(ASV, Genus)

joined_unique_asv_data <- inner_join(filtered_unique_asv, isa_ibd_table_filt_asv_genus, by = "ASV")

bubble <- ggplot(joined_unique_asv_data, aes(x = inflammation_with_surgery, y = Genus)) + 
          geom_point(aes(size = mean_abundance, fill = Genus), alpha = 0.75, shape = 21) + 
          scale_size_continuous(limits = c(0.000001, 0.1), range = c(1,17), breaks = c(1,10,50,75)) + 
          labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Genus")

   
bubble

ggsave("isa_bubble.png", bubble)

