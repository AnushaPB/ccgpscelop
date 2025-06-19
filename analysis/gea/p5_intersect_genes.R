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
    dplyr::rename(scaffold = V1, start = V2, end = V3, cds_start = V5, cds_end = V6, mrna_name = V13) %>%
    # Remove duplicate genes
    # (if you don't you will get an error later with the joins)
    # (needed for GEA genes since multiple SNPs can fall in the same gene so genes will be duplicated)
    select(scaffold, start, end, cds_start, cds_end, mrna_name) %>%
    # Extract mrna ID from the full name
    mutate(mrna = str_extract(mrna_name, "Parent=mrna-([0-9]+)")) %>%
    mutate(mrna = str_extract(mrna, "[0-9]+")) %>%
    distinct()

  # Load information about genes to link mrna to full gene name
  gene_info <- 
    get_gene_structure() %>% 
    dplyr::select(scaffold, mrna, full_name, gene_start, gene_end, cds_start, cds_end) %>%
    distinct()

  # Join gene positions with gene info to get full gene names
  gene_pos <- 
    mrna_pos %>% 
    left_join(gene_info, by = c("mrna", "cds_start", "cds_end", "scaffold")) %>%
    # Drop CDS info and just keep the gene 
    # (i.e., there are multiple coding regions in a gene, but we just want the gene name)
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
        gene_id <- gsub(".*ID=(.*?);.*", "\\1", entry)
        # Extract gene name and organism
        name_org <- gsub(".*Name=(.*?);.*", "\\1", entry)
        # Extract organism name from square brackets
        organism <- gsub(".*\\[(.*?)\\].*", "\\1", name_org)
        # Get the UniProt ID (everything before the :)
        uniprot_id <- gsub("^(.*?):.*", "\\1", organism)
        # Clean up organism name to remove any UniProt IDs
        organism <- gsub("^.*?:", "", organism)
        # Extract gene name (everything before the square bracket)
        gene_name <- gsub("\\s*\\[.*\\].*$", "", name_org)
        
        return(data.frame(ID = gene_id, gene_name = gene_name, organism = organism, name_org = name_org, uniprot_id = uniprot_id, full_name = entry))
    }
    
    # Apply the extraction function to each row
    result <- map(gene_info, extract_info, .progress = TRUE) 
    
    # Convert to data frame with proper column names
    result_df <- bind_rows(result)

    # Replace rows where ID is not present with NA
    result_df$organism <- ifelse(grepl("ID=", result_df$organism), NA, result_df$organism)
    result_df$uniprot_id <- ifelse(grepl("ID=", result_df$uniprot_id), NA, result_df$uniprot_id)
    
    return(result_df)
}


# Get positions for GEA genes
#intersect_genes("pca")
genes <- intersect_genes("bio1ndvi")

# Summary Statistics:                                                           
#  -------------------
#  Number of RDA + Linked SNPs: 1543125
#  Number of SNPs in genes: 34409 (all mutation types)
#  Number of unique genes: 17672
#  Number of unique genes with UniProtID (GO genes): 7988

# Get positions for all genes
all_genes <- get_all_genes_bed()
  
# Clean gene names
all_gene_info <- unique(all_gene_pos$full_name)
length(all_gene_info) # 74,808 genes
all_genes_clean <- clean_gene_annotations(all_gene_info)
nrow(all_genes_clean) # Confirming 74,808 genes

all_genes_org <- 
  all_genes_clean %>% 
  # Remove rows without uniprot_id
  drop_na(uniprot_id) %>% 
  # Edit ID to remove GNX-
  mutate(ID = gsub("GNX-", "", ID))  %>% 
  # remove version from uniprot_id
  mutate(uniprot_id = gsub("\\..*", "", uniprot_id)) 

nrow(all_genes_org) # 30,349 genes with uniprotids
nrow(distinct(all_genes_org, full_name)) # 30,349 unique original entries
nrow(distinct(all_genes_org, uniprot_id)) # 18,357 unique uniprot IDs

# Get repeated uniprot_ids
repeated_uniprot_ids <- 
  all_genes_org %>%
  group_by(uniprot_id) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  ungroup() %>%
  arrange(desc(n))

# Check original entries for repeated uniprot_ids
# (e.g. each has it's own GNX ID, but the name/uniprotID is the same)
repeated_uniprot_ids %>% filter(uniprot_id == "O75132") %>% dplyr::select(full_name)

# Write out full list
write_csv(all_genes_org, here("analysis", "gea", "outputs", "all_genes_list.csv"))

# Write uniprot IDs to text file
writeLines(all_genes_org$uniprot_id, here("analysis", "gea", "outputs", "all_genes_uniprot_ids.txt"))

# Get all genes not in GEA genes
all_genes_not_in_gea <- all_gene_pos %>% filter(!full_name %in% genes$full_name)
write_csv(all_genes_not_in_gea, here("analysis", "gea", "outputs", "all_genes_not_in_gea.csv"))
