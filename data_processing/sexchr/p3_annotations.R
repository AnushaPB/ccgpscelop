library(tidyverse)
library(here)

gff <- 
  read_tsv(
    here("data", "genome", "annotation", "complete.genomic.gff"),
    comment = "#",
    col_names = FALSE,
    col_types = cols(.default = "c")
  ) %>%
  select(chr = X1, full_name = X9) %>%
  filter(grepl("ID=gene", full_name)) %>%
  mutate(description = str_extract(full_name, "(?<=description=)[^;]+")) %>%
  mutate(gene_name = str_extract(full_name, "(?<=Name=)[^;]+")) %>%
  select(chr, gene_name, description)

sex_scaff <- 
  gff %>% 
  filter(chr %in% c("sex_linked_1", "sex_linked_2")) 

sex_scaff %>%
  filter(grepl("Y-linked", description))
  
sex_scaff %>%
  filter(grepl("X-linked", description))

# =============================
# Panel of Y-linked and X-linked genes
# =============================

# --- Y-linked markers ---
# Core single-copy or well-established markers
y_core <- c("SRY","ZFY","KDM5D","SMCY","UTY","DDX3Y","DBY",
            "RPS4Y1","RPS4Y2","EIF1AY","USP9Y","PRKY","AMELY")

# Multi-copy families often used in panels
y_fams <- c("RBMY","PRY","BPY2","BPY","CDY1","CDY2",
            "HSFY","XKRY","VCY","TSPY","TTTY")

# Other Y-linked protein-coding genes
y_other <- c("PCDH11Y","NLGN4Y","TBL1Y","TGIF2LY",
             "TMSB4Y","TXLNGY","TSPYL")

# Mouse-specific Y markers (optional, cross-reference)
y_mouse <- c("EIF2S3Y","Eif2s3y","Sry","Uty","Ddx3y","Kdm5d")

# Unified Y panel
y_markers <- unique(c(y_core, y_fams, y_other, y_mouse))


# --- X-linked markers (controls) ---
# Human X homologs of Y genes
x_markers <- c("ZFX","KDM5C","DDX3X","USP9X","RPS4X",
               "EIF1AX","PRKX","AMELX")

# - word-ish boundaries via non-alnum guards
# - allow family members (e.g., RBMY1A1)
# - avoid DAZAP1 when searching DAZ (negative lookahead)
rx_y  <- paste0("(^|[^A-Z0-9])(",
                paste(c(y_markers, "(?:(?<![A-Z])DAZ(?!AP))"), collapse="|"),
                ")([^A-Z0-9]|$)")
rx_x  <- paste0("(^|[^A-Z0-9])(", paste(x_markers, collapse="|"), ")([^A-Z0-9]|$)")

sex_lizards <- c("DMRT1",   # Z-linked master gene in birds; candidate in reptiles
              "SOX9",    # testis pathway activator
              "AMH",     # Anti-Müllerian Hormone, male pathway
              "FOXL2",   # ovarian pathway gene
              "CYP19A1", "Foxl2") # aromatase, estrogen synthesis, ovary development
rx_lizards <- paste0("(^|[^A-Z0-9])(", paste(sex_lizards, collapse="|"), ")([^A-Z0-9]|$)")

tagged <- sex_scaff %>%
  mutate(
    marker_class = case_when(
      str_detect(description, regex(rx_y, ignore_case = TRUE)) ~ "Y-candidate",
      str_detect(description, regex(rx_x, ignore_case = TRUE)) ~ "X-candidate",
      str_detect(description, regex(rx_lizards, ignore_case = TRUE)) ~ "Lizard-candidate",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(marker_class))

# Quick look at what matched
tagged %>% distinct(description, marker_class, chr) %>% arrange(marker_class) %>% print(n = 100)

tagged %>% group_by(chr) %>% count(marker_class)
