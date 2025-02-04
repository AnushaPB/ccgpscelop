get_smcpp <- function(folder = "outputs", K = 9, i = 10){
  # returns NULL if csv doesn't exist
  safe_csv <- possibly(read.csv)

  combos <- expand.grid(K = 1:K, i = 1:i)

  get_file <- function(K, i){
    folder <- paste0("outputs", i)
    path <- here("analysis", "smcpp", "outputs", folder, paste0("pop", K), paste0("pop", K, ".csv"))
    df <- safe_csv(path)
    df$i <- i
    return(df)
  }

  result <- pmap(combos, get_file)
  
  # drops anything NULL
  result <- compact(result)

  # binds as data.frame
  result <- 
    result %>%
    bind_rows(result)

  # converts from generation times to years
  # note: I confirmed that this is the correct conversion by running this code and comparing the output csvs:
  # smc++ plot "$output_path/smcpp_plot.pdf" "$output_path/model.final.json" -c
  # smc++ plot "$output_path/smcpp_plot.pdf" "$output_path/model.final.json" -c -g 1
  # smc++ plot "$output_path/smcpp_plot.pdf" "$output_path/model.final.json" -c -g 2
  # (first two x are identical, third x is doubled, y stays the same in all cases)
  # generation time of 2 years for S. occidentalis (Jameson & Allison, 1976).
  result <- 
    result %>%
    mutate(years1 = x, years2 = x*2, years3 = x*3) %>% 
    pivot_longer(c(years1, years2, years3), names_to = "T", values_to = "years") %>%
    mutate(T = gsub("years", "", T)) %>%
    mutate(pop = gsub("pop", "",  label))
 
  return(result)
}
