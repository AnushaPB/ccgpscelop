
library(gprofiler2)

# Get two gene IDs for mouse/human
gorth(query = c("Q8BLA8"), source_organism = "mmusculus", target_organism = "hsapiens")
gorth(query = c("Q8K424"), source_organism = "mmusculus", target_organism = "hsapiens")
ids <- c("Q8BLA8", "Q8K424")

# When you run gorth, you will get nothing because the mouse IDs are used instead of the human ids
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

# You can convert to humans like this
gorth(query = c("Q8BLA8"), source_organism = "mmusculus", target_organism = "hsapiens")
gorth(query = c("Q8K424"), source_organism = "mmusculus", target_organism = "hsapiens")
# THe gorth says the ortholog name is this, which can be used to query gprofiler2
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
go_result$result %>% head()
