
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


run_go <- function(query, org_key){  
    go_result <-
        gost(
            query = query, 
            organism = org_key, ordered_query = FALSE, 
            multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
            measure_underrepresentation = FALSE, evcodes = TRUE, 
            user_threshold = 0.05, correction_method = "g_SCS", 
            domain_scope = "annotated", custom_bg = NULL, 
            numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE
        )

    return(go_result)
}


go_fish <- function(query_size, term_size, intersection_size, query_size_all, intersection_size_all) {
  # Calculate the counts for the contingency table
  a <- intersection_size  # Number of GEA genes in the category
  b <- query_size - intersection_size  # Number of GEA genes not in the category
  c <- intersection_size_all - intersection_size  # Number of non-GEA Sceloporus genes in the category
  d <- query_size_all - intersection_size_all  # Number of non-GEA Sceloporus genes not in the category
  
  # Create the contingency table
  contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  # Perform Fisher's exact test
  result <- fisher.test(contingency_table)
  
  # Return the p-value
  return(result$p.value)
}
