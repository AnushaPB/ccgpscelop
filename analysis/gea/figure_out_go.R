
library(gprofiler2)
gorth(query = c("Q8BLA8"), source_organism = "mmusculus",
target_organism = "hsapiens")
gorth(query = c("Q8K424"), source_organism = "mmusculus",
target_organism = "hsapiens", numeric_ns = "UNIPROT_GN")

ids <- c("Q8BLA8", "Q8K424")

gorth(query = c("Q8K424"), source_organism = "mmusculus",
target_organism = "hsapiens", numeric_ns = "UNIPROT_GN")

ids <- "TRPV3"
go_result <-
  gost(
      query = ids, 
      organism = "hsapiens", 
      ordered_query = FALSE, 
      multi_query = FALSE, 
      significant = TRUE,
      exclude_iea = FALSE, 
      measure_underrepresentation = FALSE, 
      evcodes = TRUE, 
      user_threshold = 0.05, 
      correction_method = "g_SCS", 
      domain_scope = "annotated", 
      custom_bg = NULL, 
      numeric_ns = "", 
      sources = NULL, 
      as_short_link = FALSE, 
      highlight = TRUE
  )

a <- go_result$result

genes_org$full_name[11514]
gorth(query = c("Q4KLT3"), source_organism = "xlaevis", target_organism = "hsapiens")
