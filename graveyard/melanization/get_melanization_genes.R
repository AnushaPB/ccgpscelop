library(tidyverse)
source(here("general_functions.R"))
source(here("analysis", "gea", "functions_selection_stats.R"))
genes_bed <- get_all_genes_bed()


melanization_genes <- c(
  "TYR", "TYRP1", "DCT",
  "PMEL", "OCA2", "SLC45A2", "SLC24A5",
  "MITF", "SOX10", "PAX3", "LEF1",
  "MC1R", "ASIP", "POMC", "ADCY3", "ADCY5", "PRKACA",
  "MYO5A", "RAB27A", "MLPH",
  "GPNMB", "TFRC", "SOD1", "SOD2",
  "HPS1", "HPS2", "HPS3", "HPS4", "HPS5", "HPS6",
  "BLOC1S1", "BLOC1S2", "BLOC1S3", "BLOC1S4", "BLOC1S5", "BLOC1S6"
)

top_genes <- MC1R
Master switch for eumelanin vs pheomelanin. Classic background-matching gene; often shifts with darker substrates (ash/char).

ASIP
Antagonist of MC1R. Changes here can flip melanization without touching MC1R itself (very common in natural systems).

MITF
Master transcription factor controlling the whole melanocyte program. Fire-linked stressors (UV, heat) often act through MITF regulation.

TYR
Rate-limiting enzyme for melanin synthesis. If pigment amount differs, TYR often shows dosage/expression shifts.

TYRP1
Tyrosinase-related protein 1. Involved in melanin production; mutations can affect pigmentation.

top_genes <- c("MC1R", "ASIP", "MITF", "TYR", "TYRP1")

melanization_genes_bed <- 
  genes_bed %>%
  filter(str_detect(full_name, paste0("\\b(", paste(top_genes, collapse = "|"), ")\\b")))

bed_out <- melanization_genes_bed %>%
  transmute(scaffold,
            start = as.integer(start) - 1L,
            end   = as.integer(end)) %>%
  distinct() %>%
  arrange(scaffold, start, end)

write_tsv(bed_out, here("analysis", "melanization", "outputs", "melanization_genes.bed"), col_names = FALSE)
