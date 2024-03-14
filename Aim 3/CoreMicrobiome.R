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

# Determine core ASVs in each group
# May need to adjust detection and prevalence values later
