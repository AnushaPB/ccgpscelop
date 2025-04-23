library(tidyverse)
library(here)
library(furrr)
outputs <- here("analysis", "gea", "outputs") 

nonsyn_frq <- read_table(here(outputs, "nonsyn.frq"))
nonsyn_vec <- nonsyn_frq$MAF
names(nonsyn_vec) <- nonsyn_frq$SNP
nonsyn_vec <- round(nonsyn_vec, 2)

neutral_frq <- read_table(here(outputs, "nogeanogenes.frq"))
neutral_vec <- neutral_frq$MAF
names(neutral_vec) <- neutral_frq$SNP
neutral_vec <- round(neutral_vec, 2)

allnonsyn <- read_table(here(outputs, "allnonsyn.frq"))
allnonsyn_vec <- allnonsyn$MAF
names(allnonsyn_vec) <- allnonsyn$SNP
allnonsyn_vec <- round(allnonsyn_vec, 2)

# Use progress bar to track progress
# library(progress)
# pb <- progress_bar$new(
#   format = "  Matching SNPs [:bar] :percent in :elapsed",
#   total = length(nonsyn_vec),
#   clear = FALSE,
#   width = 60
# )

# # Initialize a flag vector for neutral_vec to track used SNPs
# used_neutral <- rep(FALSE, length(neutral_vec))
# pulled_snps <- rep(NA, length(nonsyn_vec))

# for (i in seq_along(nonsyn_vec)) {
#   pb$tick()  # Update progress bar
  
#   # Find candidate indices in neutral_vec that match the current nonsyn SNP and haven't been used
#   candidate_indices <- which(neutral_vec == nonsyn_vec[i] & !used_neutral)
  
#   if (length(candidate_indices) == 0) {
#     warning(paste("No match found for SNP:", names(nonsyn_vec)[i]))
#     pulled_snps[i] <- NA
#   } else {
#     # Use the first available match
#     match_index <- candidate_indices[1]
#     pulled_snps[i] <- names(neutral_vec)[match_index]
#     used_neutral[match_index] <- TRUE  # Mark as used
#   }
# }

# Sample matching frequencies from sample of 10 Million
neutral_vec_1M <- sample(neutral_vec, 10000000)

plan(multisession, workers = 5)
match_frq <- future_map_chr(nonsyn_vec, ~{
  names(sample(which(neutral_vec_1M == .x), 1))
}, .progress = TRUE, .options = furrr_options(seed = TRUE))
plan(sequential)

unique_match <- unique(match_frq)
length(match_frq) - length(unique_match) # Check non-matches
length(unique_match)
writeLines(unique_match, here(outputs, "neutral_snp_frqmatch.txt"))

# Repeat for all non=synonymous SNPs
library(furrr)
plan(multisession, workers = 5)
match_frq2 <- future_map_chr(nonsyn_vec, ~{
  names(sample(which(allnonsyn_vec == .x), 1))
}, .progress = TRUE, .options = furrr_options(seed = TRUE))
plan(sequential)

unique_match2 <- unique(match_frq2)
length(match_frq2) - length(unique_match2) # Check non-matches
length(unique_match2)
writeLines(unique_match2, here(outputs, "allnonsyn_snp_frqmatch.txt"))

