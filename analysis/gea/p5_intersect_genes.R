library(tidyverse)
library(here)
source(here("general_functions.R"))
source(here("analysis", "gea", "functions_selection_stats.R"))

# Function to get original locus names from GEA SNPs in genes
intersect_genes <- function(prefix){
  # Outpath
  outpath <- here("analysis", "gea", "outputs")
  
  # Read in GEA genes
  path <- here(outpath, paste0(prefix, "_gea_genes.bed"))
  message("Reading in GEA genes from ", path)
  genes <- read.table(path, sep="\t", quote="", fill=TRUE, stringsAsFactors=FALSE)
  
  # Get positions for SNPs in genes 
  mrna_pos <- 
    genes %>%
    dplyr::rename(
      scaffold  = V1,
      start     = V2,
      end       = V3,
      cds_start = V5,
      cds_end   = V6,
      description = V13
    ) %>%
    select(scaffold, start, end, cds_start, cds_end, description) %>%
    # strip any leading line numbers (defensive)
    mutate(attrs = stringr::str_replace(description, "^[0-9]+\\s+", "")) %>%
    # extract locus_tag for joining to gene_info
    mutate(
      locus_tag = stringr::str_match(attrs, "(?:^|;)locus_tag=([^;]+)")[,2]
    ) %>%
    distinct()
  
  # Link to full gene info via locus_tag
  gene_info <- 
    get_gene_structure() %>%
    dplyr::select(scaffold, locus_tag, full_name, gene_start, gene_end, cds_start, cds_end) %>%
    distinct()

  # Join SNP positions to genes
  gene_pos <- 
    mrna_pos %>%
    dplyr::left_join(
      gene_info,
      by = c("scaffold", "locus_tag", "cds_start", "cds_end")
    ) %>%
    # keep the gene-level fields you need
    dplyr::select(scaffold, start, end, gene_start, gene_end, full_name) %>%
    distinct()

  # Make sure all genes have names
  stopifnot(all(complete.cases(gene_pos$full_name)))

  # Get RDA positions (these have SNP names)
  rda <-
    read_csv(here(outpath, paste0(prefix, "_rda_ids.csv"))) %>%
    # Convert to 0-based start/end bed format
    mutate(start = position - 1, end = position)

  # Intersect the two to get SNP names for the SNPs in genes
  gene_snp <- left_join(gene_pos, rda, by = c("scaffold", "start", "end"))

  # gene_snp will have duplicate snp rows because snps can fall in multiple overlapping genes
  # write out a distinct list
  gene_snp_distinct <- gene_snp %>% distinct(locus) %>% pull(locus)

  # Write out just the locus names for plink (heterozygosity calcs)
  write.table(gene_snp_distinct, here(outpath, paste0(prefix, "_gea_genes_ids.txt")), quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Write out the full file for GO analysis
  write.csv(gene_snp, here(outpath, paste0(prefix, "_gea_genes_snp.csv")), row.names = FALSE)

  # Clean gene names
  gene_info <- unique(gene_snp$full_name)
  genes_clean <- clean_gene_annotations(gene_info)
  genes_org <- 
    genes_clean %>% 
    # Remove rows without uniprot_id
    drop_na(uniprot_id) %>% 
    # Edit ID to remove GNX-
    mutate(ID = gsub("GNX-", "", ID))  %>% 
    # remove version from uniprot_id
    mutate(uniprot_id = gsub("\\..*", "", uniprot_id)) 

  write_csv(genes_org, here(outpath, paste0(prefix, "_genes_list.csv")))

  # Calculate statistics
  summary_stats <- data.frame(
    Statistic = c(
      "Number of RDA + Linked SNPs", 
      "Number of SNPs in genes", 
      "Number of unique genes", 
      "Number of unique genes with UniProtID (GO genes)"),
    Value = c(
      length(unique(rda$locus)), 
      nrow(distinct(gene_pos, scaffold,  start, end)), 
      gene_pos$full_name %>% unique() %>% length(),
      genes_org$full_name %>% unique() %>% length()
    )
  )

  # Print statistics
  cat(
    "Summary Statistics:\n",
    "-------------------\n",
    sprintf("Number of RDA + Linked SNPs: %d\n", summary_stats$Value[1]),
    sprintf("Number of SNPs in genes: %d\n", summary_stats$Value[2]),
    sprintf("Number of unique genes: %d\n", summary_stats$Value[3]),
    sprintf("Number of unique genes with UniProtID (GO genes): %d\n", summary_stats$Value[4])
  )

  # Save statistics to a file
  write_csv(summary_stats, here(outpath, paste0(prefix, "_summary_stats.csv")))

  return(genes_org)
}

# Clean gene annotations
clean_gene_annotations <- function(gene_info) {
  # Function to extract information from each entry
  extract_info <- function(entry) {
    # Extract Gene ID
    gene_id <- gsub(".*ID=([^;]+)(;.*|$)", "\\1", entry)
    # Extract gene name
    gene_name <- gsub(".*Name=([^;]+)(;.*|$)", "\\1", entry)
    if (grepl("egapxtmp_", gene_name)) gene_name <- NA
    # Extract description
    description <- gsub(".*description=([^;]+)(;.*|$)", "\\1", entry)
    if (!grepl("description=", entry))   description <- NA

    return(data.frame(ID = gene_id, gene_name = gene_name, description = description, full_name = entry))
  }
  
  # Apply the extraction function to each row
  result <- map(gene_info, extract_info, .progress = TRUE) 
  
  # Convert to data frame with proper column names
  result_df <- bind_rows(result)

  return(result_df)
}


# Get positions for GEA genes
#intersect_genes("pca")
genes <- intersect_genes("bio1ndvi")
# Summary Statistics:                                                           
#  -------------------
#  Number of RDA + Linked SNPs: 1871451
#  Number of SNPs in genes: 23051
#  Number of unique genes: 9076
#  Number of unique genes with UniProtID (GO genes): 9076

# Get positions for all genes
all_genes <- get_all_genes_bed()
  
# Clean gene names
all_gene_info <- unique(all_genes$full_name)
length(all_gene_info) # 20,638 (note old had 74,808 genes)
all_genes_clean <- clean_gene_annotations(all_gene_info)
nrow(all_genes_clean) # Confirming 20,638 genes

# Check named
all_genes_clean %>% filter(!is.na(gene_name)) %>% nrow()

# All named genes have descriptions
all_genes_clean %>% filter(!is.na(gene_name), is.na(description)) %>% nrow()

# But not all genes with descriptions have names
missing_gene_name <- all_genes_clean %>% filter(is.na(gene_name), !is.na(description)) 
missing_gene_name %>% nrow() # 5535 genes with descriptions but no gene_name

# This is a problem because some genes of interest have descriptions but no gene_name
all_genes_clean %>% filter(grepl("shock 70", description)) %>% head()

all_genes_org <- 
  all_genes_clean %>% 
  # Remove rows without gene name or description
  filter(!is.na(gene_name) | !is.na(description)) 

nrow(all_genes_org) # 18670
nrow(drop_na(all_genes_org, gene_name)) #13135

# Write out full list
write_csv(all_genes_org, here("analysis", "gea", "outputs", "all_genes_list.csv"))
