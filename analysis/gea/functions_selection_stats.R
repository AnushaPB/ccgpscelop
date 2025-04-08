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

get_tajimad <- function(windowkb = 100){
  pattern <- paste0(windowkb, "kb_tajimad.*\\.Tajima\\.D$")
  tjd_files <- list.files(here("analysis", "gea", "outputs"), pattern = pattern, full.names = TRUE)
  tjd_pops <- str_extract(tjd_files, "pop[0-9]+") %>% unique()
  names(tjd_files) <- tjd_pops

  tjd <- 
    map(tjd_files, read_table) %>% 
    bind_rows(.id = "pop") %>%
    mutate(TajimaD = as.numeric(TajimaD)) %>%
    rename(scaffold = CHROM) %>%
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


get_genes <- function(){
  path <- here("analysis", "gea", "outputs", "all_genes.bed")
  message("Reading in ", path)
  read.table(path, sep="\t", quote="", fill=TRUE, stringsAsFactors=FALSE) %>%
  rename(scaffold = V1, start = V4, end = V5, full_name = V9)
}

get_gea_genes <- function(nonsyn = FALSE){
  path <- here("analysis", "gea", "outputs", "bio1ndvi_gea_gene_snp.csv")
  message("Reading in ", path)
  gea <- read_csv(path)

  gea <- 
    gea %>% 
    mutate(
      id = str_extract(full_name, "GNX-\\d+"),
      name = str_extract(full_name, "(?<=Name=)[^\\[]+"),
      uniprotid = str_extract(full_name, "(?<=\\[)[^:]+")
    )

  message("Number of GEA genes: ", nrow(gea))

  if (nonsyn){
    path2 <- here("analysis", "gea", "outputs", "bio1ndvi_gea_gene_nonsyn_ids.txt")
    message("Getting nonsynonymous genes from: ", path2)
    nonsyn_genes <- readLines(path2)
    nonsyn_gea <- gea %>% filter(locus %in% nonsyn_genes)
    message("Number of nonsynonymous GEA genes: ", nrow(nonsyn_gea))
    return(nonsyn_gea)
  }

  return(gea)
}

get_goi_names <- function(){
  goi_names <- c("Heat shock", "Inositol 1,4,5-trisphosphate receptor type 1", "Inositol 1,4,5-trisphosphate receptor type 2", "Neuropeptide FF receptor 1", "Calcitonin receptor", "Growth hormone receptor", "Transient receptor potential cation channel subfamily V member 3", "Transient receptor potential cation channel subfamily A member 1", "TRPM8", "Endothelial PAS domain-containing protein 1")
  names(goi_names) <- c(
    "HSP",              # Heat shock → HSP (Heat Shock Protein)
    "IP3R1",            # Inositol 1,4,5-trisphosphate receptor type 1
    "IP3R2",            # Inositol 1,4,5-trisphosphate receptor type 2
    "NPFFR1",           # Neuropeptide FF receptor 1
    "CALCR",            # Calcitonin receptor
    "GHR",              # Growth hormone receptor
    "TRPV3",            # TRP cation channel subfamily V member 3
    "TRPA1",            # TRP cation channel subfamily A member 1
    "TRPM8",            # Already abbreviated
    "EPAS1"             # Endothelial PAS domain-containing protein 1
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
    filter(reduce(map(goi_names, ~str_detect(full_name, .x)), `|`)) %>%
    filter(grepl("chr", scaffold)) %>%
    rename(CHROM = scaffold)

  if (type == "snps") {message("Returning SNPs, to return genes set type to genes"); return(goi)}

  # Get full gene start/end
  all_genes <- get_genes()
  goi_genes <- all_genes %>% filter(full_name %in% goi$full_name)

  if (type == "genes") {message("Returning genes, to return SNPs, set type to snps"); return(goi_genes)}
  stop("Invalid type. Choose 'genes' or 'snps'.")
}

get_moddf <- function(){
  path <- here("analysis", "genetic_diversity", "outputs", "model_df.csv")
  message("Reading in ", path)
  read_csv(path)
}
