

get_ortho <- function(genes_org, target_organism, file_prefix, cache = TRUE){

  path <- here("analysis", "gea", "outputs", paste0(file_prefix, "_gene_orthologs.csv"))

  if (file.exists(path) & cache == TRUE) {
    message("Using cached file", path)
    return(read_csv(path))
  }

  organisms <- unique(genes_org$organism)
  names(organisms) <- organisms

  ortho_result <- 
    map(
      organisms,
      ~{
        convert_ortho(
          query = genes_org %>% filter(organism == .x) %>% pull(uniprot_id),
          organism = .x,
          target_organism = target_organism
        )
      }, .progress = TRUE
    )

  # Combine orthologs and original human genes to get the final query names
  organism_name <- 
    case_when(
      target_organism == "hsapiens" ~ "Homo sapiens", 
      target_organism == "ggallus" ~ "Gallus gallus",
      TRUE ~ NA)

  if (is.na(organism_name)) stop ("Organism must be ggallus or hsapiens")

  original <- filter(genes_org, grepl(organism_name, organism)) %>% dplyr::select(uniprot_id)
  genes_and_orthos <- 
    ortho_result %>%
    compact() %>%
    bind_rows() %>%
    dplyr::select(ortholog_name, uniprot_id = query) %>%
    bind_rows(original) %>%
    mutate(query = case_when(is.na(ortholog_name) ~ uniprot_id, TRUE ~ ortholog_name)) 

  # Count uniprotid duplicates (uniprot_id with more than one ortholog)
  dup <- genes_and_orthos %>% count(uniprot_id) %>% filter(n > 1) %>% nrow()
  if (dup > 0) warning(dup, " uniprot_ids with more than one ortholog, taking first ortholog")

  # Take the first ortholog for each uniprot_id so no genes are overrepresented because they have more orthologs 
  genes_and_orthos <- 
    genes_and_orthos %>%
    group_by(uniprot_id) %>%
    dplyr::slice(1) %>%
    ungroup()

  # Write out the orthologs
  message("Writing out orthologs to ", path)
  write_csv(genes_and_orthos, path)

  return(genes_and_orthos)
}

convert_ortho <- function(query, organism, target_organism){
  # Convert to gprofiler format
  source_organism <- convert_to_orgcode(organism)

  # Check if the source and target organisms are the same, if so, return NULL
  if (source_organism == target_organism) return(NULL)

  # Get the orthologs
  safe_gorth <- possibly(gorth)
  ortho <- safe_gorth(query = query, source_organism = source_organism, target_organism = target_organism)

  # Return NULL if no orthologs are found
  if (is.null(ortho)) return(NULL)

  df <- data.frame(ortholog_name = ortho$ortholog_name, query = ortho$input, source_organism = source_organism, target_organism = target_organism, description = ortho$description)

  return(df)
}

convert_to_orgcode <- function(organism) {
  # Extract the genus and species
  genus_species <- str_extract(organism, "^[A-Za-z]+ [a-z]+")
  
  # Extract the first letter of the genus and the full species name
  genus <- substr(genus_species, 1, 1)
  species <- str_extract(genus_species, " [a-z]+")
  
  # Concatenate the first letter of the genus and the species name
  short_code <- tolower(paste0(genus, species))
  
  # Remove the space
  short_code <- gsub(" ", "", short_code)
  
  return(short_code)
}


run_go <- function(query, custom_bg, organism, correction_method = "fdr") {
  gost(
    query = query, 
    organism = organism, 
    ordered_query = FALSE, 
    multi_query = FALSE, 
    significant = FALSE, 
    exclude_iea = FALSE, 
    measure_underrepresentation = FALSE, 
    evcodes = TRUE, 
    user_threshold = 0.05, 
    correction_method = correction_method,
    domain_scope = "custom_annotated",
    custom_bg = custom_bg, 
    numeric_ns = "", 
    sources = NULL, 
    as_short_link = FALSE, 
    highlight = FALSE
  )
}

# go_correction <- function(go_result, all_genes){
#   gea_genes <- 
#     go_result$result %>% 
#     as_tibble() %>%
#     # Just keep GO terms
#     filter(grepl("GO", source)) %>%
#     dplyr::select(intersection_size, query_size, term_size, term_name, source, term_id, p_value)

#   # Note: if you get an error here it means that there are NAs in query_size, which should not be the case - make sure no filtering was done on all_genes (all genes should not be filtered by p)
#   comparison <- 
#     right_join(all_genes, gea_genes, by = "term_id") %>%
#     rowwise() %>%
#     mutate(fish_p = go_fish(query_size, intersection_size, query_size_all, intersection_size_all))
    
#   # Term size all and term size should be the same
#   stopifnot(comparison$term_size_all == comparison$term_size)

#   # Get significant terms 
#   # *precision:* The precision of the enrichment, calculated as the ratio of intersection_size to query_size. Precision measures the proportion of relevant genes (i.e., genes associated with the GO term) in the query set.
#   # *recall:* The recall of the enrichment, calculated as the ratio of intersection_size to term_size. Recall measures the proportion of genes associated with the GO term that are present in the reference set.
#   sig <- 
#     comparison %>% 
#     # Correct p-values using fdr
#     mutate(fish_p = p.adjust(fish_p, method = "fdr")) %>%
#     filter(fish_p < 0.05) %>% 
#     arrange(fish_p) %>% 
#     mutate(precision = intersection_size / query_size,
#            recall = intersection_size / term_size) %>%
#     dplyr::select(term_name, term_id, source, fish_p, query_size_all, intersection_size, precision, recall)

#   return(sig)
# }

# go_fish <- function(query_size, intersection_size, query_size_all, intersection_size_all) {
#   # Calculate the counts for the contingency table
#   a <- intersection_size  # Number of GEA genes in the category
#   b <- query_size - intersection_size  # Number of GEA genes not in the category
#   c <- intersection_size_all - intersection_size  # Number of non-GEA Sceloporus genes in the category
#   d <- query_size_all - intersection_size_all  # Number of non-GEA Sceloporus genes not in the category
  
#   # Create the contingency table
#   contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
#   # Perform Fisher's exact test (one-sided for enrichment)
#   result <- fisher.test(contingency_table, alternative = "greater")
  
#   # Return the p-value
#   return(result$p.value)
# }



go_table <- function(df, n = 5) {
  # Rename columns
  df <- 
    df %>%
    dplyr::rename(
      `p` = fish_p,
      `Precision` = precision,
      `Recall` = recall,
      `Term name` = term_name,
      `Term ID` = term_id,
      `No. genes` = intersection_size
    ) %>%
    mutate(`Term Name` = str_to_sentence(`Term name`)) %>%
    head(n) %>%
    dplyr::select(`Term name`, `Term ID`, Precision, Recall, `No. genes`, p) 
  
  # Create gt table
  gt_table <- 
    df %>%
    gt() %>%
    fmt_scientific(
      columns = p,
      decimals = 2
    ) %>%
    fmt_number(
      columns = c(`Precision`, `Recall`),
      decimals = 2
    ) %>%
    data_color(
      columns = c(`Precision`, `Recall`),
      fn = scales::col_numeric(
        palette = RColorBrewer::brewer.pal(9, "YlGnBu"),
        domain = c(min(df$Precision, df$Recall), max(df$Precision, df$Recall))
      )
    )
  
  return(gt_table)
}

go_sig <- function(x, category, p = 0.05){
  x %>% 
    filter(source == category) %>%
    dplyr::select(term_name, intersection_size, precision, recall, p_value) %>%
    filter(p_value < p) %>%
    arrange(p_value)
}

go_sigbp <- function(x, p = 0.05) {
  go_sig(x, "GO:BP", p)
}

go_sigmf <- function(x, p = 0.05) {
  go_sig(x, "GO:MF", p)
}

go_top <- function(x, category, n = 10, arrange_by = "intersection_size", descending = TRUE, precision_cutoff = 0) {
  message("Filtering ", category, " terms with precision >= ", precision_cutoff)
  message("Arranging by ", arrange_by, " in ", if (descending) "descending" else "ascending", " order")
  message("Returning top ", n, " terms")
  x %>%
    filter(source == category) %>%
    filter(precision >= precision_cutoff) %>%
    dplyr::select(term_name, intersection_size, precision, recall, p_value) %>%
    arrange(if (descending) dplyr::desc(!!sym(arrange_by)) else !!sym(arrange_by)) %>%
    head(n)
}

go_topbp <- function(x, n = 10, arrange_by = "precision", descending = TRUE, precision_cutoff = 0) {
  go_top(x, "GO:BP", n, arrange_by, descending, precision_cutoff)
}

go_topmf <- function(x, n = 5, arrange_by = "precision", descending = TRUE, precision_cutoff = 0) {
  go_top(x, "GO:MF", n, arrange_by, descending, precision_cutoff)
}


