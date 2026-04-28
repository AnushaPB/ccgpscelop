library(here)
library(tidyverse)
source(here("general_functions.R"))

coords <- read_table(here("data", "ccgp_data", "58-Sceloporus.coords.txt"), col_names = FALSE)
nrow(coords)
vcfids <- read_csv(here("data_processing", "vcf_samples.txt"), col_names = "SampleID")
x <- read_csv(here("data_processing", "sample_removal_decisions.csv")) %>% rename(SampleID = initial_all_samps)

"Sceocc_MNW16-025-S33" %in% vcfids
x %>% filter(SampleID == "Sceocc_MNW16-025-S33") %>% select(details)

x %>% count(details)
#Scelocci_CAS236233

setdiff(coords$X1, x$SampleID)

kept_inds <- 
  x %>% 
  filter(details == "Retained in final assembly")
removed_inds <-
  x %>% 
  filter(details != "Retained in final assembly")


# Samples that were removed before creating the VCF file
pre_vcf <- x %>% filter(!SampleID %in% vcfids$SampleID)
pre_vcf %>% select(SampleID, details) 
# Just look at ones that were removed for reasons other than the filter
pre_vcf %>% select(SampleID, details) %>% filter(!grepl("filter", details))
# 1 Scelocci_CHI1382_DAW5-46-21 Unknown provenance

# Samples that need to be removed after creating the VCF file
post_vcf <- x %>% filter(SampleID %in% setdiff(vcfids$SampleID, kept_inds$SampleID))
post_vcf %>% select(SampleID, details)

# SUMMARY:
message("Started with ", nrow(x), " samples in the initial sample list.")
# Started with 298 samples in the initial sample list.

message("Removed ", nrow(pre_vcf), " samples before creating the VCF file.")
# Removed 32 samples before creating the VCF file

message("Of these ", nrow(pre_vcf), " samples, ", sum(grepl("filter", pre_vcf$details)), " were removed due to filtering, and ", sum(!grepl("filter", pre_vcf$details)), " were removed due to unknown provenance.")
#Of these 32 samples, 31 were removed due to filtering, and 1 was removed due to unknown provenance.

message("So the VCF file was created with ", nrow(vcfids), " samples.")
# So the VCF file was created with 266 samples.
stopifnot(nrow(vcfids) == nrow(x) - nrow(pre_vcf)) # Checking this

message("After creating the VCF file, we removed an additional ", nrow(post_vcf), " samples.")
# After creating the VCF file, we removed an additional 9 samples.

message("Of these ", nrow(post_vcf), " samples, ", sum(grepl("beckii", post_vcf$details)), " were removed due to being the beckii species, ", sum(grepl("swapped", post_vcf$details)), " were removed due to potential sample swaps based on pop structure")
# Of these 9 samples, 7 were removed due to being the beckii species, and 2 were removed due to potential sample swaps based on pop structure

message("Resulting in a final count of ", nrow(kept_inds), " samples in the final VCF file.")
# Resulting in a final count of 257 samples in the final VCF file.
stopifnot(nrow(kept_inds) == nrow(vcfids) - nrow(post_vcf)) # Checking this
