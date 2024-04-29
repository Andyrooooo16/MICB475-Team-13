# Load libraries
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggVennDiagram)


##### Basic Core microbiome not used in final project #####
##### Scroll to "Final" section #####

# Load unrarefied phyloseq object
load("ibd_filt.RData")

# Transform sample counts to relative abundance
ibd_RA <- transform_sample_counts(ibd_filt, fun=function(x) x/sum(x))

# Filter dataset by surgery status
ibd_nosurg <- subset_samples(ibd_RA, `cd_resection` == "no")
ibd_surg <- subset_samples(ibd_RA, `cd_resection` == "yes")
                                  
# Filter dataset by inflammation status
ibd_noinf <- subset_samples(ibd_RA, `inflammation` == FALSE)
ibd_inf <- subset_samples(ibd_RA, `inflammation` == TRUE)     

# Filter dataset by CD location!! Need to add code for this

# Determine core ASVs in each group
# May need to adjust detection and prevalence values later
nosurg_ASVs <- core_members(ibd_nosurg, detection = 0, prevalence = 0.7)
surg_ASVs <- core_members(ibd_surg, detection = 0, prevalence = 0.7)

noinf_ASVs <- core_members(ibd_noinf, detection = 0, prevalence = 0.7)
inf_ASVs <- core_members(ibd_inf, detection = 0, prevalence = 0.7)

# Determine if ASV IDs meet criteria specified
nosurg_ASVs
surg_ASVs

noinf_ASVs
inf_ASVs

# Determine which taxa ASVs belong to
tax_table(prune_taxa(nosurg_ASVs, ibd_filt))
tax_table(prune_taxa(surg_ASVs, ibd_filt))
                                  
tax_table(prune_taxa(noinf_ASVs, ibd_filt))
tax_table(prune_taxa(inf_ASVs, ibd_filt))

# Plot relative abundance of ASVs
prune_taxa(nosurg_ASVs,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`cd_resection`, scales ="free")

prune_taxa(noinf_ASVs,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`inflammation`, scales ="free")

# Combine nosurg and surg into a list
surg_list_full <- list(No_Surgery = nosurg_ASVs, Surgery = surg_ASVs)
surg_list_full

# Combine noinf and inf into a list
inf_list_full <- list(No_Inflammation = noinf_ASVs, Inflammation = inf_ASVs)
inf_list_full

# Generate Venn diagram for surgery & save
surg_venn <- ggVennDiagram(x = surg_list_full)
surg_venn

ggsave("venn_surg.png", surg_venn)  

# Generate Venn diagram for surgery & save
inf_venn <- ggVennDiagram(x = inf_list_full)
inf_venn

ggsave("venn_inf.png", inf_venn)  
                    

##### FINAL Core microbiome Inflammation and Surgery #####

# Load unrarefied phyloseq object
load("ibd_filt.RData")

# Transform sample counts to relative abundance
ibd_RA <- transform_sample_counts(ibd_filt, fun=function(x) x/sum(x))

# Filter dataset by surgery status
ibd_nosurg_noinflam <- subset_samples(ibd_RA, `inflammation_with_surgery` == "no_inflammation_no_surgery")
ibd_surg_noinflam <- subset_samples(ibd_RA, `inflammation_with_surgery` == "no_inflammation_with_surgery")
ibd_nosurg_inflam <- subset_samples(ibd_RA, `inflammation_with_surgery` == "inflammation_no_surgery")
ibd_surg_inflam <- subset_samples(ibd_RA, `inflammation_with_surgery` == "inflammation_with_surgery")


# Filter dataset by CD location!! Need to add code for this

# Determine core ASVs in each group
# May need to adjust detection and prevalence values later
ibd_nosurg_noinflam_ASVs <- core_members(ibd_nosurg_noinflam, detection = 0, prevalence = 0.7)
ibd_surg_noinflam_ASVs <- core_members(ibd_surg_noinflam, detection = 0, prevalence = 0.7)
ibd_nosurg_inflam_ASVs <- core_members(ibd_nosurg_inflam, detection = 0, prevalence = 0.7)
ibd_surg_inflam_ASVs <- core_members(ibd_surg_inflam, detection = 0, prevalence = 0.7)

# Determine which taxa ASVs belong to
tax_table(prune_taxa(ibd_nosurg_noinflam_ASVs, ibd_filt))
tax_table(prune_taxa(ibd_surg_noinflam_ASVs, ibd_filt))
tax_table(prune_taxa(ibd_nosurg_inflam_ASVs, ibd_filt))
tax_table(prune_taxa(ibd_surg_inflam_ASVs, ibd_filt))

# Plot relative abundance of ASVs
prune_taxa(ibd_nosurg_noinflam_ASVs,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`inflammation_with_surgery`, scales ="free")

prune_taxa(ibd_surg_noinflam_ASVs,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`inflammation_with_surgery`, scales ="free")

prune_taxa(ibd_nosurg_inflam_ASVs,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`inflammation_with_surgery`, scales ="free")

prune_taxa(ibd_surg_inflam_ASVs,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`inflammation_with_surgery`, scales ="free")

# Combine nosurg and surg into a list
inf_surg_list_full <- list(Inflammation_Surgery = ibd_surg_inflam_ASVs, 
                           Inflammation_No_Surgery = ibd_nosurg_inflam_ASVs, 
                           No_Inflammation_No_Surgery = ibd_nosurg_noinflam_ASVs,
                           No_Inflammation_Surgery = ibd_surg_noinflam_ASVs)

# Combine noinf and inf into a list
inf_list_full <- list(No_Inflammation = noinf_ASVs, Inflammation = inf_ASVs)
inf_list_full

# Generate Venn diagram for surgery status and inflammation status
final_venn <- ggVennDiagram(x = inf_surg_list_full) +
    scale_fill_gradient(low = "white", high = "orange")
final_venn
ggsave("final_venn.png", final_venn) 


##### FINAL CD_LOCATION Not sure if we will use #####

# Load unrarefied phyloseq object
load("ibd_filt.RData")

# Transform sample counts to relative abundance
ibd_RA <- transform_sample_counts(ibd_filt, fun=function(x) x/sum(x))

# Filter dataset by surgery status
ibd_L1 <- subset_samples(ibd_RA, `cd_location` == "Ileal (L1)")
ibd_L2 <- subset_samples(ibd_RA, `cd_location` == "Colonic (L2)")
ibd_L3 <- subset_samples(ibd_RA, `cd_location` == "Ileocolonic (L3)")

# Determine core ASVs in each group
# May need to adjust detection and prevalence values later
ibd_L1_ASVs <- core_members(ibd_L1, detection = 0, prevalence = 0.7)
ibd_L2_ASVs <- core_members(ibd_L2, detection = 0, prevalence = 0.7)
ibd_L3_ASVs <- core_members(ibd_L3, detection = 0, prevalence = 0.7)

# Determine which taxa ASVs belong to
tax_table(prune_taxa(ibd_L1_ASVs, ibd_filt))
tax_table(prune_taxa(ibd_L2_ASVs, ibd_filt))
tax_table(prune_taxa(ibd_L3_ASVs, ibd_filt))


# Plot relative abundance of ASVs
prune_taxa(ibd_L1_ASVs,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`cd_location`, scales ="free")

prune_taxa(ibd_L2_ASVs,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`cd_location`, scales ="free")

prune_taxa(ibd_L3_ASVs,ibd_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`cd_location`, scales ="free")

# Combine nosurg and surg into a list
location_list_full <- list(L1_Ileal = ibd_L1_ASVs, 
                           L2_Colonic = ibd_L2_ASVs, 
                           L3_Ileal_Colonic = ibd_L3_ASVs)

# Generate Venn diagram for surgery status and inflammation status
location_venn <- ggVennDiagram(x = location_list_full) # color
location_venn
ggsave("venn_location.png", inf_surg_venn)

### Final venn w/ only species of interest