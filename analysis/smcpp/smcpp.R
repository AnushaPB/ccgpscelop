get_smcpp <- function(folder = "outputs"){
  # returns NULL if csv doesn't exist
  safe_csv <- possibly(read.csv)
  
  result <- map(1:5, ~safe_csv(here("analysis", "smcpp", paste0("pop", .x), paste0("pop", .x, ".csv")))) 
  
  # drops anything that is not a data.frame
  # NULL cases will be lists
  result <- result[map_lgl(result, is.data.frame)]

  # binds as data.frame
  result <- 
    result %>%
    #FIX INDEX IS WRONG IF NOT ALL POPS ARE PRESENT
    imap(~mutate(.x, pop = .y)) %>%
    bind_rows(result)

  # converts from generation times to years
  # note: I confirmed that this is the correct conversion by running this code and comparing the output csvs:
  # smc++ plot "$output_path/smcpp_plot.pdf" "$output_path/model.final.json" -c
  # smc++ plot "$output_path/smcpp_plot.pdf" "$output_path/model.final.json" -c -g 1
  # smc++ plot "$output_path/smcpp_plot.pdf" "$output_path/model.final.json" -c -g 2
  # (first two x are identical, third x is doubled, y stays the same in all cases)
  result <- 
    result %>%
    mutate(years1 = x, years2 = x*2, years3 = x*3) %>% 
    pivot_longer(c(years1, years2, years3), names_to = "T", values_to = "years") %>%
    mutate(T = gsub("years", "", T))

  return(result)
}
