format_dist <- function(){
  d <- read.table(here("58-Sceloporus", paste0("58-Sceloporus_plinkdist_1ibs.mdist")))
  id <- read.table(here("58-Sceloporus", paste0("58-Sceloporus_plinkdist_1ibs.mdist.id")))[,2]
  rownames(d) <- colnames(d) <- id
  path <- here("data", paste0("58-Sceloporus_dist.csv"))
  write.csv(d, path)
  message("wrote dist file to:", path)
}

get_gendist <- function(){
  gendist <- read.csv(here("data", paste0("58-Sceloporus_dist.csv")), row.names = 1)
  colnames(gendist) <- row.names(gendist)
  return(gendist)
}
