library(tidyverse)
library(here)
library(tidyverse)
library(tidyverse)
library(here)

# NCBI-flavored GFF3 -> Ensembl-ish IDs for bcftools csq
# - Keeps gene_biotype / transcript_biotype if present (NCBI style)
# - If transcript has CDS but no biotype, infers protein_coding
# - Emits: gene, transcript (from mRNA), and child features (exon/CDS/UTRs)
convert_gff_ncbi <- function(input_file, output_file) {
  message("Converting GFF3 from NCBI style to Ensembl style")

  cols <- c("chrom","source","feature_type","start","end","score","strand","phase","attributes")

  gff <- readr::read_tsv(
    input_file,
    comment = "#",
    col_names = cols,
    col_types = readr::cols(.default = readr::col_character()),
    progress = FALSE
  ) %>%
    mutate(
      .row   = row_number(),
      ID     = stringr::str_match(attributes, "(^|;)ID=([^;]+)")[,3],
      Parent = stringr::str_match(attributes, "(^|;)Parent=([^;]+)")[,3],
      gene_biotype       = stringr::str_match(attributes, "(^|;)gene_biotype=([^;]+)")[,3],
      transcript_biotype = stringr::str_match(attributes, "(^|;)transcript_biotype=([^;]+)")[,3]
    )

  # Gene map: original gene ID -> "gene:<ID>"
  gene_map <- gff %>%
    filter(feature_type == "gene") %>%
    transmute(gene_key = ID, new_gene_id = paste0("gene:", ID))

  # Transcripts (from mRNA), attach parent gene
  transcripts <- gff %>%
    filter(feature_type == "mRNA") %>%
    mutate(new_transcript_id = paste0("transcript:", ID)) %>%
    left_join(gene_map, by = c("Parent" = "gene_key")) %>%
    mutate(new_parent_gene = coalesce(new_gene_id, paste0("gene:", Parent)))

  # Mark transcripts that have at least one CDS
  tx_with_cds <- gff %>%
    filter(feature_type == "CDS", !is.na(Parent)) %>%
    distinct(Parent) %>%
    transmute(tx_id = Parent, has_cds = TRUE)

  transcripts <- transcripts %>%
    left_join(tx_with_cds, by = c("ID" = "tx_id")) %>%
    mutate(
      has_cds = replace_na(has_cds, FALSE),
      biotype_final = coalesce(
        transcript_biotype,                                # prefer transcript biotype
        gff$gene_biotype[match(Parent, gff$ID)],           # else gene's biotype
        if_else(has_cds, "protein_coding", NA_character_)  # else infer if CDS
      )
    )

  # Children: exon/CDS/UTRs -> Parent should be transcript:<ID>
  child_feats <- gff %>%
    filter(feature_type %in% c("exon", "CDS", "three_prime_UTR", "five_prime_UTR")) %>%
    left_join(transcripts %>% transmute(ID, new_transcript_id),
              by = c("Parent" = "ID")) %>%
    mutate(new_parent_tx = coalesce(new_transcript_id, paste0("transcript:", Parent)))

  # --- minimal vectorized helpers to preserve attributes ---
  rewrite_id <- function(attr, new_id) {
    stringr::str_replace(attr, "(^|;)ID=[^;]+", paste0("\\1ID=", new_id))
  }
  rewrite_parent <- function(attr, new_parent) {
    has_parent <- stringr::str_detect(attr, "(^|;)Parent=[^;]+")
    replaced   <- stringr::str_replace(attr, "(^|;)Parent=[^;]+", paste0("\\1Parent=", new_parent))
    dplyr::if_else(has_parent, replaced, paste0(attr, ";Parent=", new_parent))
  }
  append_biotype <- function(attr, biotype_val) {
    has_bio <- stringr::str_detect(attr, "(^|;)biotype=[^;]+")
    to_add  <- dplyr::if_else(is.na(biotype_val) | biotype_val == "", "", paste0(";biotype=", biotype_val))
    paste0(attr, dplyr::if_else(!has_bio & to_add != "", to_add, ""))
  }
  append_phase_if_cds <- function(attr, is_cds, phase_val) {
    paste0(attr, dplyr::if_else(is_cds, paste0(";Phase=", phase_val), ""))
  }

  # ---------- Build output rows ----------
  genes_out <- gff %>%
    filter(feature_type == "gene") %>%
    mutate(
      attributes = rewrite_id(attributes, paste0("gene:", ID)),
      attributes = append_biotype(attributes, gene_biotype)
    ) %>%
    select(chrom, source, feature_type, start, end, score, strand, phase, attributes, .row)

  transcripts_out <- transcripts %>%
    mutate(
      attributes = rewrite_id(attributes, new_transcript_id),
      attributes = rewrite_parent(attributes, new_parent_gene),
      attributes = append_biotype(attributes, biotype_final)
    ) %>%
    transmute(
      chrom, source,
      feature_type = "transcript",
      start, end, score, strand, phase,
      attributes,
      .row
    )

  children_out <- child_feats %>%
    mutate(
      attributes = rewrite_parent(attributes, new_parent_tx),
      attributes = append_phase_if_cds(attributes, feature_type == "CDS", phase)
    ) %>%
    select(chrom, source, feature_type, start, end, score, strand, phase, attributes, .row)

  out <- bind_rows(genes_out, transcripts_out, children_out) %>%
    arrange(.row) %>%
    select(-.row)

  # Perform conversion check
  conversion_check(original = gff, converted = out)

  message("Writing converted GFF3 to: ", output_file)
  readr::write_lines("##gff-version 3", output_file)
  readr::write_tsv(out, output_file, append = TRUE, col_names = FALSE)
  invisible(output_file)
}

conversion_check <- function(original, converted) {
  n_genes_orig <- sum(original$feature_type == "gene")
  n_genes_conv <- sum(converted$feature_type == "gene")
  if (n_genes_orig != n_genes_conv) {
    stop("Gene count mismatch after conversion: ", n_genes_orig, " vs ", n_genes_conv)
  }

  n_tx_orig <- sum(original$feature_type == "mRNA")
  n_tx_conv <- sum(converted$feature_type == "transcript")
  if (n_tx_orig != n_tx_conv) {
    stop("Transcript count mismatch after conversion: ", n_tx_orig, " vs ", n_tx_conv)
  }

  # Spot check a random set of genes
  set.seed(91196)
  message("Spot checking 100 genes for conversion accuracy...")
  map(1:100, ~{
    check_gene1 <- sample(original$ID[original$feature_type == "gene"], 1)
    original_gene1 <- original %>% filter(feature_type == "gene", ID == check_gene1) 
    converted_gene1 <- converted %>% filter(feature_type == "gene", str_detect(attributes, paste0("ID=gene:", check_gene1)))

    # Check for all columns matching (except attributes which may be reordered)
    direct_comparison <- colnames(converted_gene1)
    direct_comparison <- direct_comparison[direct_comparison != "attributes"]
    if (!all(original_gene1[, direct_comparison] == converted_gene1[, direct_comparison])){
      stop("Gene record mismatch after conversion for gene ID: ", check_gene1)
    }

    # Check attributes content
    converted_gene1_attrs <- 
      str_split(converted_gene1$attributes, ";")[[1]] %>% 
      str_split_fixed("=", n = 2) %>%
      as_tibble(.name_repair = "minimal") %>%
      set_names(c("key", "value")) %>%
      arrange(key)

    original_gene1_attrs <- 
      str_split(original_gene1$attributes, ";")[[1]] %>% 
      str_split_fixed("=", n = 2) %>%
      as_tibble(.name_repair = "minimal") %>%
      set_names(c("key", "value")) %>%
      arrange(key) 

    # Check that all original attributes (except ID) are present and matching
    attrs_check <-
      converted_gene1_attrs %>%
      filter(key %in% original_gene1_attrs$key) %>%
      # Replace `:` with `-` for key=ID to match original for comparison
      mutate(value = if_else(key == "ID", gsub("gene:", "", value), value)) %>%
      pull(value) == original_gene1_attrs$value

    if (!all(attrs_check)) {
      stop("Gene record mismatch after conversion for gene ID: ", check_gene1)
    }
  }, .progress = TRUE)

  message("Conversion check passed: ", n_genes_conv, " genes and ", n_tx_conv, " transcripts.")
}


input_file <- here("data", "genome", "annotation", "complete.genomic.gff")
output_file <- here("data", "genome", "annotation", "converted_annotation.gff")
convert_gff_ncbi(input_file, output_file)
