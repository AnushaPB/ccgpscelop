get_pi <- function(){
  # Read in windowed pi
  wpi <- read_table(here("analysis", "genetic_diversity", "outputs", "58-Sceloporus_10kb_windowpi.windowed.pi"))

  # Pull out long scaffolds (chr and Scaffold_13)
  i  <- 
    wpi %>% 
    count(CHROM) %>% 
    arrange(desc(n)) %>% 
    filter(n > 1000) %>% 
    pull(CHROM)

  # Filter and add a cumulative sum column
  df <-
    wpi %>%
    filter(CHROM %in% i) %>%
    mutate(CUMSUM = cumsum(BIN_START)) %>%
    mutate(CHROM = factor(CHROM, levels = i)) %>%
    rename(scaffold = CHROM)
  
  return(df)
}

get_tajimad <- function(windowkb = 10){
  pattern <- paste0(windowkb, "kb_tajimad.*\\.Tajima\\.D$")
  tjd_files <- list.files(here("analysis", "gea", "outputs"), pattern = pattern, full.names = TRUE)
  tjd_pops <- str_extract(tjd_files, "pop[0-9]+") %>% unique()
  names(tjd_files) <- tjd_pops

  tjd <- 
    map(tjd_files, read_table) %>% 
    bind_rows(.id = "pop") %>%
    mutate(TajimaD = as.numeric(TajimaD)) %>%
    dplyr::rename(scaffold = CHROM) %>%
    filter(grepl("chr", scaffold)) %>%
    mutate(scaffold = factor(scaffold, levels = paste0("chr", 1:11))) %>%
    group_by(pop) %>%
    mutate(top5 = if_else(TajimaD > quantile(TajimaD, 0.95, na.rm = TRUE), 1, 0), bottom5 = if_else(TajimaD < quantile(TajimaD, 0.05, na.rm = TRUE), 1, 0)) %>%
    mutate(outliers = case_when(top5 == 1 ~ "Top 5%", bottom5 == 1 ~ "Bottom 5%", TRUE ~ NA)) %>%
    ungroup() %>% 
    mutate(BIN_END = BIN_START + windowkb * 1000)
  
  return(tjd)
}

get_chromenv <- function(){
  path <- here("analysis", "chromenv", "outputs", "chromenv_sem.csv")
  message("Reading in ", path)
  read_csv(path)
}

get_nonsyn <- function(){
  path <- here("analysis", "gea", "outputs", "nonsynonymous.csv")
  message("Reading in ", path)
  read_csv(path)
}

get_syn <- function(){
  path <- here("analysis", "gea", "outputs", "synonymous.csv")
  message("Reading in ", path)
  read_csv(path)
}

get_all_genes_bed <- function(){
  path <- here("analysis", "gea", "outputs", "all_genes.bed")
  message("Reading in ", path)
  genes <- 
    read.table(path, sep="\t", quote="", fill=TRUE, stringsAsFactors=FALSE) %>%
    dplyr::rename(scaffold = V1, start = V2, end = V3, full_name = V10) %>%
    dplyr::select(scaffold, start, end, full_name)
  return(genes)
}

get_gene_structure <- function(){
  genes <- 
    get_all_genes_bed() %>%
    mutate(mrna = str_extract(full_name, "ID=GNX-([0-9]+)")) %>%
    mutate(mrna = str_extract(mrna, "[0-9]+")) %>%
    rename(gene_start = start, gene_end = end) %>%
    mutate(gene_name = str_extract(full_name, "Name=([^;]+)")) %>%
    mutate(gene_name = str_replace(gene_name, "Name=", "")) %>%
    mutate(gene_name = str_remove(gene_name, "\\[.*\\]"))

  path <- here("analysis", "gea", "outputs", "all_exons.bed")
  message("Reading in ", path)

  exons <-
    read.table(path, sep = "\t", quote = "", fill = TRUE, stringsAsFactors = FALSE)  %>%
    dplyr::rename(scaffold = V1, exon_start = V2, exon_end = V3, type = V8, info = V10) %>%
    mutate(mrna = str_extract(info, "Parent=mrna-([0-9]+)")) %>%
    mutate(mrna = str_extract(mrna, "[0-9]+")) %>%
    dplyr::select(scaffold, exon_start, exon_end, mrna) %>%
    left_join(genes, by = c("scaffold", "mrna"))

  path <- here("analysis", "gea", "outputs", "all_cds.bed")
  message("Reading in ", path)

  cds <-
    read.table(path, sep = "\t", quote = "", fill = TRUE, stringsAsFactors = FALSE) %>%
    rename(scaffold = V1, cds_start = V2, cds_end = V3, info = V10) %>%
    mutate(mrna = str_extract(info, "Parent=mrna-([0-9]+)")) %>%
    mutate(mrna = str_extract(mrna, "[0-9]+")) %>%
    dplyr::select(scaffold, cds_start, cds_end, mrna) %>%
    # Relationship is many to many because exons can have multiple CDS
    left_join(exons, by = c("scaffold", "mrna"), relationship = "many-to-many") 

  head(cds)

  # Check that matching worked
  stopifnot(all(complete.cases(cds$full_name)))

  return(cds)
}

get_all_genes_uniprotid <- function(){
  path <- here("analysis", "gea", "outputs", "all_genes_list.csv")
  message("Reading in ", path)
  df <- read_csv(path)
  return(df)
}

get_gea_genes <- function(nonsyn = FALSE){
  path <- here("analysis", "gea", "outputs", "bio1ndvi_gea_genes_snp.csv")
  message("Reading in ", path)
  gea <- 
    read_csv(path) %>%
    mutate(gene_name = str_extract(full_name, "Name=([^;]+)")) %>%
    mutate(gene_name = str_replace(gene_name, "Name=", "")) %>%
    mutate(gene_name = str_remove(gene_name, "\\[.*\\]"))

  message("Number of GEA SNPs in genes: ", nrow(gea))
  message("Number of unique genes with GEA SNPs: ", length(unique(gea$full_name)))

  if (nonsyn){
    outpath <- here("analysis", "gea", "outputs")
    path2 <-  here(outpath, "bio1ndvi_gea_genes_nonsyn.csv")
    message("Getting nonsynonymous genes from: ", path2)
    nonsyn <- 
      read_csv(path2) %>% 
      # Join with gea to get cleaned gene IDs
      left_join(gea, by = c("scaffold", "start", "end", "gene_start", "gene_end", "full_name", "locus", "position"))
    message("Number of nonsynonymous GEA SNPs in genes: ", nrow(nonsyn))
    message("Number of unique nonsynonymous genes with GEA SNPs: ", length(unique(nonsyn$full_name))) 
    return(nonsyn)
  }

  return(gea)
}

get_goi_names <- function(){
  goi_names <- c(
    "AllHSP70" = "Heat shock 70",                     # General HSP70 group
    "HSPA13" = "Heat shock 70 kDa protein 13", 
    "HSP30C" = "Heat shock protein 30C", 
    "HSPA12B" = "Heat shock 70 kDa protein 12B", 
    "HSFY" = "Heat shock transcription factor, Y-linked", 
    "HSPB2" = "Heat shock protein beta-2", 
    "HSPA1B" = "Heat shock 70 kDa protein 1B", 
    "HSPA1A" = "Heat shock 70 kDa protein 1A",
    "ITPR1" = "Inositol 1,4,5-trisphosphate receptor type 1", 
    "ITPR2" = "Inositol 1,4,5-trisphosphate receptor type 2", 
    "NPFFR1" = "Neuropeptide FF receptor 1", 
    "CALCR" = "Calcitonin receptor", 
    "GHR" = "Growth hormone receptor", 
    "TRPV3" = "Transient receptor potential cation channel subfamily V member 3", 
    "TRPA1" = "Transient receptor potential cation channel subfamily A member 1", 
    "TRPM8" = "TRPM8", 
    "EPAS1" = "Endothelial PAS domain-containing protein 1",
    "DNAJB6" = "DnaJ homolog subfamily B member 6"
  )

  return(goi_names)
}

get_goi <- function(type = "genes"){
  gea <- get_gea_genes()

  goi_names <- get_goi_names()
  message("Genes of interest: ", paste(goi_names, collapse = ", "))

  # Get GEA SNPs that fall within GOI
  goi <- 
    gea %>% 
    filter(purrr::reduce(map(goi_names, ~str_detect(full_name, .x)), `|`)) %>%
    filter(grepl("chr", scaffold)) %>%
    dplyr::rename(CHROM = scaffold)

  if (type == "snps") {message("Returning SNPs, to return genes set type to genes"); return(goi)}

  # Get full gene start/end
  all_genes <- get_all_genes_bed()
  goi_genes <- all_genes %>% filter(full_name %in% goi$full_name)

  if (type == "genes") {message("Returning genes, to return SNPs, set type to snps"); return(goi_genes)}
  stop("Invalid type. Choose 'genes' or 'snps'.")
}

get_moddf <- function(){
  path <- here("analysis", "genetic_diversity", "outputs", "model_df.csv")
  message("Reading in ", path)
  read_csv(path)
}

get_cor_snps_info <- function(labels = FALSE){
  path <- here("analysis", "gea", "outputs", "bio1ndvi_rda_linked_snps_info.csv")
  df <- read_csv(path) 

  # Get all outliers with linked SNPs
  linked <- df %>% drop_na(linked_locus)

  # Get non-synonymous SNPs in genes
  nonsyn <- get_gea_genes(nonsyn = TRUE)
  nonsyn_info <- nonsyn %>% dplyr::select(locus, full_name)

  # Join with linked information
  linked_info <- 
    linked %>% 
    # Get genes associated with linked SNPs
    left_join(nonsyn_info, by = c("linked_locus" = "locus")) %>%
    rename(linked_gene = full_name) %>%
    # Get genes associated with outlier SNPs
    left_join(nonsyn_info, by = c("outlier_locus" = "locus")) %>%
    rename(outlier_gene = full_name) %>%
    # Drop any where both linked_gene and outlier_gene are NA
    filter(!is.na(linked_gene) | !is.na(outlier_gene)) 

  # Collapse into just any gene associated with the outlier locus either directly or indirectly
  outlier_info <-
    linked_info %>%
    dplyr::select(scaffold, outlier_locus, linked_gene, outlier_gene) %>%
    pivot_longer(c(linked_gene, outlier_gene), names_to = "gene_link", values_to = "full_name") %>%
    # Remove any NA rows
    drop_na(full_name) %>%
    # Get unique genes
    distinct(scaffold, outlier_locus, full_name) %>%
    ungroup() %>%
    # Rejoin with nonsyn to get more info
    rename(locus = outlier_locus) %>%
    # Many-to-many relationship because one gene can be linked to multiple SNPs and one SNP can be linked to multiple genes
    # Remove locus and scaffold columns from nonsyn since the correct ones are in the main df
    left_join(dplyr::select(nonsyn, -locus, -scaffold), by = "full_name", relationship = "many-to-many")

  # Count how many loci are linked to multiple genes
  outlier_info %>%
    group_by(locus) %>%
    summarise(n_genes = n()) %>%
    count(n_genes > 1)
  # Count how many genes are linked to multiple loci
  outlier_info %>%
    group_by(full_name) %>%
    summarise(n_loci = n()) %>%
    count(n_loci > 1)

  # ADD BACK LOCI WITHOUT LINKED SNPS
  # Check to make sure columns are the same
  stopifnot(length(setdiff(names(nonsyn), names(outlier_info))) == 0)
  outlier_info_joined <- 
    outlier_info %>%
    bind_rows(nonsyn) %>%
    distinct()

  if (labels){
    # For labeling, make a combined label with all gene names associated with an outlier
    labels <-
      outlier_info %>%
      group_by(locus) %>%
      summarise(
        gene_name = paste(unique(gene_name), collapse = "/")
      ) %>%
      # DROP ANY GENES THAT DON'T HAVE A NAME
      filter(gene_name != "NA")
    return(labels)
  }

  return(outlier_info_joined)
}
