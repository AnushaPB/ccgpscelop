library(tidyverse)
library(here)

# Get positions for SNPs in genes
genes <- read.table(here("analysis", "gea", "outputs", "gea_genes.bed"), sep="\t", quote="", fill=TRUE, stringsAsFactors=FALSE)
gene_pos <- 
  genes %>% 
  dplyr::select(V1, V2, V3, V12) %>% 
  dplyr::rename(scaffold = V1, start = V2, end = V3, full_name = V12) %>%
  # Remove duplicates (if you don't you will get an error later with the joins)
  distinct()

# Get RDA positions which have SNP names
rda <- 
  read_csv(here("analysis", "gea", "outputs", "rda_sig_p01.csv")) 

# Intersect the two to get SNP names for the SNPs in genes
gene_snp <- left_join(gene_pos, rda, by = c("scaffold", "start", "end"))

# Write out just the locus names for plink
write.table(gene_snp$locus, here("analysis", "gea", "outputs", "gea_gene_ids.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

write.csv(gene_snp, here("analysis", "gea", "outputs", "gea_gene_snp.csv"), row.names = FALSE)

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

write_csv(genes_org, here("analysis", "gea", "outputs", "genes_list.csv"))

# Calculate statistics
cat(
  "Summary Statistics:\n",
  "-------------------\n",
  sprintf("Number of RDA SNPs: %d\n", nrow(rda)),
  sprintf("Number of SNPs in genes: %d\n", nrow(gene_pos)),
  sprintf("Number of unique genes: %d\n", gene_pos$full_name %>% unique() %>% length()),
  sprintf("Number of unique genes with UniProtID (GO genes): %d\n", genes_org$full_name %>% unique() %>% length())
)
