# Investigate RAM usage for RDA

# There are three values for each of the analyses: 
#     (1) reading/processing input data
#     (2) running the do_everything function
#     (3) exporting the resulting files

# For each of these, peakRAM() calculates three things: 
#     (1) elapsed time (secs)
#     (2) total RAM used (MiB)
#     (3) peak RAM used (MiB)

# Let's read in and bind together all the values we've run so far,
# also retrieving the dataset sizes (which should somewhat correspond
# to computational time / RAM usage).

library(tidyverse)
library(cowplot)
library(here)
theme_set(theme_cowplot())

# 151 scaffolds
filenames <- list.files(path = here("outputs/RDA"), pattern = "*peakRAM*", recursive = TRUE, full.names = TRUE)
shortfiles <- list.files(path = here("outputs/RDA"), pattern = "*peakRAM*", recursive = TRUE, full.names = FALSE)

# TODO only 24 scaffolds have results (due to scaffold length filters?)
dat <-
    1:length(shortfiles) %>% 
    lapply(function(x) {
        y <- readr::read_csv(filenames[x])
        y <- y %>% dplyr::mutate(filename = paste0(shortfiles[x]))
        return(y)
    }) %>% 
    dplyr::bind_rows() %>% 
    tidyr::separate_wider_delim(cols = filename, delim = "/", names = c("scaffold", "temp")) %>% 
    dplyr::select(-temp)
  
dat <- dat %>% 
    dplyr::mutate("Elapsed_Time_mins" = Elapsed_Time_sec/60) %>% 
    dplyr::select(-Elapsed_Time_sec) %>% 
    tidyr::pivot_longer(cols = c("Elapsed_Time_mins", "Total_RAM_Used_MiB", "Peak_RAM_Used_MiB"), names_to = "metric", values_to = "value")

p <- dat %>% 
    ggplot(aes(x = scaffold, y = value, fill = fxn)) +
    geom_bar(position = "stack", stat = "identity") +
    facet_wrap(~metric, scales = "free") +
    # geom_point(aes(x = scaffold, y = snps/100)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
p
ggsave(here("analysis/anne/outputs/RDA_RAM_usage.pdf"), width = 11, height = 6)

# To add on secondary panel that has numbers of SNPs (getting some NAs, not quite working)
rdadapt <- read_csv(here("outputs/RDA/58-Sceloporus_RDA_outliers_full_rdadapt.csv"), col_names = TRUE)

site_nos <- rdadapt %>% 
    group_by(scaff) %>% 
    summarize(snps = n()) %>% 
    rename(scaffold = scaff)

dat <- left_join(dat, site_nos)

p2 <- dat %>% 
    ggplot(aes(x = scaffold, y = snps)) +
    geom_point() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

plot_grid(p, p2, nrow = 1)
ggsave(here("analysis/anne/outputs/RDA_RAM_usage_composite.pdf"), width = 11, height = 6)