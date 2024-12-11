format_dist <- function(){
  path <- here("analysis", "ibdibe", "outputs")
  d <- read.table(here(path,"58-Sceloporus_annotated_pruned_0.6.mdist"))
  id <- read.table(here(path,"58-Sceloporus_annotated_pruned_0.6.mdist.id"))[,2]
  rownames(d) <- colnames(d) <- id
  path <- here(path, "58-Sceloporus_dist.csv")
  write.csv(d, path)
  message("wrote dist file to:", path)
}

get_gendist <- function(){
  path <- here("analysis", "ibdibe", "outputs")
  gendist <- read.csv(here(path, "58-Sceloporus_dist.csv"), row.names = 1)
  colnames(gendist) <- row.names(gendist)
  return(gendist)
}

unfold <- function(X, scale = TRUE) {
  x <- vector()
  for (i in 2:nrow(X)) x <- c(x, X[i, 1:i - 1])
  if (scale == TRUE) x <- scale(x, center = TRUE, scale = TRUE)
  return(x)
}

dist_to_df <- function(Y, X, stdz = TRUE){
  #Unfold X and Y
  y <- unfold(Y, scale = stdz)
  dfX <- purrr::map_dfc(X, unfold, scale = stdz) %>% purrr::map_dfc(as.numeric)

  # Make single variable dataframe
  df <- dfX %>%
    dplyr::mutate(Y = y) %>%
    tidyr::gather("var", "X", -Y)

  return(df)
}
