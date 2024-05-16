format_dist <- function(){
  d <- read.table(here("data", "processed_data","58-Sceloporus_plinkdist_1ibs.mdist"))
  id <- read.table(here("data", "processed_data","58-Sceloporus_plinkdist_1ibs.mdist.id"))[,2]
  rownames(d) <- colnames(d) <- id
  path <- here("data", "processed_data", "58-Sceloporus_dist.csv")
  write.csv(d, path)
  message("wrote dist file to:", path)
}

get_gendist <- function(){
  gendist <- read.csv(here("data", "processed_data", "58-Sceloporus_dist.csv"), row.names = 1)
  colnames(gendist) <- row.names(gendist)
  return(gendist)
}
