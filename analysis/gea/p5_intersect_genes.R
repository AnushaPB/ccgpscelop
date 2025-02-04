library(tidyverse)
library(here)

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
        
        return(data.frame(ID = gene_id, gene_name = gene_name, organism= organism, full_name= name_org, uniprot_id = uniprot_id, original_entry = entry))
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

intersect_genes <- function(prefix){
  # Outpath
  outpath <- here("analysis", "gea", "outputs")
  
  # Read in GEA genes
  genes <- read.table(here(outpath, paste0(prefix, "_gea_genes.bed")), sep="\t", quote="", fill=TRUE, stringsAsFactors=FALSE)
  
  # Get positions for SNPs in genes
  gene_pos <- 
    genes %>% 
    dplyr::select(V1, V2, V3, V12) %>% 
    dplyr::rename(scaffold = V1, start = V2, end = V3, full_name = V12) %>%
    # Remove duplicates (if you don't you will get an error later with the joins)
    distinct()

  # Get RDA positions which have SNP names
  rda <- read_csv(here(outpath, paste0(prefix, "_rda_ids.csv")))

  # Intersect the two to get SNP names for the SNPs in genes
  gene_snp <- left_join(gene_pos, rda, by = c("scaffold", "start", "end"))

  # Write out just the locus names for plink (heterozygosity calcs)
  write.table(gene_snp$locus, here(outpath, paste0(prefix, "_gea_gene_ids.txt")), quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Write out the full file for GO analysis
  write.csv(gene_snp, here(outpath, paste0(prefix, "_gea_gene_snp.csv")), row.names = FALSE)

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
  cat(
    "Summary Statistics:\n",
    "-------------------\n",
    sprintf("Number of RDA + Linked SNPs: %d\n", nrow(rda)),
    sprintf("Number of SNPs in genes: %d\n", nrow(gene_pos)),
    sprintf("Number of unique genes: %d\n", gene_pos$full_name %>% unique() %>% length()),
    sprintf("Number of unique genes with UniProtID (GO genes): %d\n", genes_org$full_name %>% unique() %>% length())
  )
}


# Get positions for GEA genes
intersect_genes("pca")
intersect_genes("bio1ndvi")

# Get positions for all genes
all_genes <- read.table(here("analysis", "gea", "outputs", "all_genes.bed"), sep="\t", quote="", fill=TRUE, stringsAsFactors=FALSE)
all_gene_pos <- 
  genes %>% 
  dplyr::select(V1, V4, V5, V9) %>% 
  dplyr::rename(scaffold = V1, start = V4, end = V5, full_name = V9) %>%
  # Remove duplicates (if you don't you will get an error later with the joins)
  distinct()
  
# Clean gene names
all_gene_info <- unique(all_gene_pos$full_name)
length(all_gene_info)
all_genes_clean <- clean_gene_annotations(gene_info)

all_genes_org <- 
  all_genes_clean %>% 
  # Remove rows without uniprot_id
  drop_na(uniprot_id) %>% 
  # Edit ID to remove GNX-
  mutate(ID = gsub("GNX-", "", ID))  %>% 
  # remove version from uniprot_id
  mutate(uniprot_id = gsub("\\..*", "", uniprot_id)) 

write_csv(all_genes_org, here("analysis", "gea", "outputs", "all_genes_list.csv"))
