# Load libraries
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggVennDiagram)

# Load unrarefied phyloseq object
load(_______)

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
                               



