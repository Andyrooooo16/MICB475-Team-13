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
  theme(#legend.position = "none",
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

getwd()
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
  group_by(OTU, Genus,Family,Order,Class,Phylum, inflammation_with_surgery) %>%
  summarize(
    mean_abundance = mean(Abundance, na.rm = TRUE)
  )

summary_by_otu_mean <- summary_by_otu_mean %>%
  rename(ASV = OTU)

filtered_unique_asv <- summary_by_otu_mean %>%
  filter(ASV %in% unique(isa_ibd_table_filt$ASV))


# We want to add Family-level information to any genus labels that are unclear: 'UCG', 'uncultured'
to_change = c('g__UCG-002','g__UCG-003','g__UCG-005','g__UCG-010','g__uncultured')

# We can use mutate() and ifelse() to quickly modify just these rows.
# ifelse() is very helpful - the first argument should be a true/false test, the second is what happens if it's true, and the third is if it's false.
joined_unique_asv_data = filtered_unique_asv %>% 
  mutate(Genus = ifelse(Genus %in% to_change, paste(Family,Genus,sep='.'),Genus))

# Upon inspection UCG-010 is actually uncharacterized at the Genus level too! 
# Let's replace the family info with order info.
joined_unique_asv_data = joined_unique_asv_data %>% 
  mutate(Genus = ifelse(Genus == 'f__UCG-010.g__UCG-010', paste(Order,"g__UCG-010",sep='.'),Genus))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bubble <- ggplot(joined_unique_asv_data, aes(x = inflammation_with_surgery, y = Genus)) + 
          geom_point(aes(size = mean_abundance, fill = Genus), alpha = 0.75, shape = 21) + 
          scale_size_continuous(limits = c(0.000001, 0.1), range = c(1,17), breaks = c(1,10,50,75)) + 
          labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Genus")
   
bubble

ggsave("isa_bubble.png", bubble)

#### EDIT: ADDING COLOURS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Add a column that tells us whether we should colour the point something other than grey.
# Index indicates the type of indicator (eg. each different combination of 1s and 0s gets a different number).
isa_colour = isa_ibd_table_filt %>% select(ASV:s.no_inflammation_with_surgery,index) %>% 
  # All the 1s and 0s were in separate columns, so let's turn them into a single column to match joined_unique_asv_data.
  pivot_longer(cols = -c(ASV,index), names_to = 'inflammation_with_surgery', values_to = 'not_gray') %>% 
  # Now we can remove the 's.' prefix on inflammation_with_surgery to make it match our dataset
  mutate(inflammation_with_surgery = str_remove(inflammation_with_surgery,'s.'))

joined2 = joined_unique_asv_data %>% left_join(isa_colour)

# Let's assign specific colours to each group of interest.
# I think we discussed colouring by column, but after inspection of isa_ibd_table_filt, you actually only have four different types of indicators (see line 138). That's not an overwhelming number, let's use that!

# If not_gray is zero, let's change index to zero too. (Importantly, zero is not already a value in the index column) 
# That essentially creates a new value of index that we can assign the colour grey.
joined2 = joined2 %>% mutate(index = ifelse(not_gray==0,0,index)) %>% 
  # And let's make it a character so that it doesn't create a colour gradient
  mutate(index = as.character(index))
# index now can be used to denote a separate colour for each indicator type, as well as gray for any groups where the taxon of interest is not considered an indicator.

# Prep the colours we want using the scales package 
# (scales:: directly refers to the package, you don't need to use library())
# First, we generate n-1 colours (since we'll be assigning gray to one of the index values)
plot_cols = c('grey30',scales::hue_pal()(length(unique(joined2$index))-1))
# Let's match each colour up with our index values
names(plot_cols) = sort(unique(joined2$index)) # Sort it so that 0 comes first, since the gray colour comes first in plot_cols
plot_cols

# Finally, let's order our taxa according to the indicator type to make it easier to see trends. Remember to take out the zeros since that's true for multiple indicator types!
order = joined2 %>% filter(index !='0') %>% arrange(index) %>% pull(Genus) %>% unique()
joined2 = joined2 %>% mutate(Genus = factor(Genus,levels = order))

# Finally, the plot:
ggplot(joined2, aes(x = inflammation_with_surgery, y = Genus)) +
  # We are now filling by index. I also changed alpha back to 1 to make the colours pop more
  geom_point(aes(size = mean_abundance, fill = index), alpha = 1, shape = 21) + 
  # your original scale_size_continuous had limits, which caused some values to be removed. It also caused the legend to disappear. The below just scales the points to make them bigger without changing the dataset.
  scale_size(range = c(2,14)) +
  labs( x= "", y = "", size = "Relative\nAbundance (%)", fill = "Indicator\nType") +
  # This tells ggplot to use the colours we specified above.
  # Guide='none' removes this part of the legend (might be easier to just explain in your figure legend, as otherwise you should change the numbers to something more descriptive)
  scale_fill_manual(values = plot_cols, guide='none')

ggsave("isa_bubble_V2.png", height = 11, width = 9)
