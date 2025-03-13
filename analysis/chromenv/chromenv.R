plotchrom_all <- function(data){
  print(plotchrom_c(data, "rsq") + ggtitle("R²")) 
  print(plotchrom_estimate(data, "csi_past") + ggtitle("Paleoclimate stability effect"))
  print(plotchrom_estimate(data, "tmean_dif")+ ggtitle("Contemporary climate change effect"))
  print(plotchrom_top(data) + ggtitle("Top variable"))
  print(plotchrom_top(filter(data, estimate < 0)) + ggtitle("Top variable (negative)"))
  print(plotchrom_unique(data) + ggtitle("Unique variable"))
  print(plotchrom_unique(filter(data, estimate < 0)) + ggtitle("Unique variable (negative)"))
}

plotchrom_c <- function(data, fill){
  scaff_rect <- get_scaff_rect(data)
  ggplot(data) +
    geom_rect(aes(xmin = BIN_START, xmax = BIN_END, ymin = 0, ymax = 1, fill = .data[[fill]])) +
    geom_rect(data = scaff_rect, aes(xmin = BIN_START, xmax = BIN_END), ymin = 0, ymax = 1, col = "black", fill = NA) +
    scaffold_theme() +
    scale_fill_viridis_c(option = "magma")
}

plotchrom_estimate <- function(data, fill){
  estimate_df <-
    data %>%
    dplyr::select(var, sig, BIN_START, BIN_END, scaffold) %>%
    pivot_wider(names_from = var, values_from = sig) 

  scaff_rect <- get_scaff_rect(data)
  plotchrom_div(estimate_df, fill)
}

plotchrom_div <- function(data, fill){
  ggplot(data) +
    geom_rect(aes(xmin = BIN_START, xmax = BIN_END, ymin = 0, ymax = 1, fill = .data[[fill]])) +
    geom_rect(data = scaff_rect, aes(xmin = BIN_START, xmax = BIN_END), ymin = 0, ymax = 1, col = "black", fill = NA) +
    scaffold_theme() +
    scale_fill_gradient2(midpoint = 0, low = "cyan3", mid = "white", high = "orange2", limits = c(-0.55, 0.55))
}

plotchrom_cat <- function(summary_data, data, fill){
  scaff_rect <- get_scaff_rect(data)
  colors <- get_colors()

  ggplot(summary_data) +
    geom_rect(data = data, aes(xmin = BIN_START, xmax = BIN_END, ymin = 0, ymax = 1), fill = "#e6e6e6") +
    geom_rect(aes(xmin = BIN_START, xmax = BIN_END, ymin = 0, ymax = 1, fill = .data[[fill]]), col = rgb(0,0,0,0)) +
    geom_rect(data = scaff_rect, aes(xmin = BIN_START, xmax = BIN_END), ymin = 0, ymax = 1, col = "black", fill = NA) +
    scaffold_theme() +
    scale_fill_manual(values = colors)
}

plotchrom_top <- function(data){
  top <-get_top(data)
  return(plotchrom_cat(top, data, "top_r"))
}

get_top <- function(data){
  top <-
    data %>%
    ungroup() %>%
    group_by(scaffold, BIN_START, BIN_END) %>%
    mutate(mag = abs(sig)) %>%
    dplyr::select(scaffold, BIN_START, BIN_END, var, mag) %>%
    mutate(top_r = case_when(all(is.na(mag)) ~ NA, mag == max(mag, na.rm = TRUE) ~ var, TRUE ~ NA)) %>%
    drop_na(top_r)

  return(top)
}

plotchrom_unique <- function(data){
  unique <- get_unique(data)
  return(plotchrom_cat(unique, data, "var"))
}

get_unique <- function(data){
  unique <-
    data %>%
    ungroup() %>%
    group_by(scaffold, BIN_START, BIN_END) %>%
    mutate(mag = abs(sig)) %>%
    dplyr::select(scaffold, BIN_START, BIN_END, var, mag) %>%
    # Filter for bins with only one significant correlation
    filter(sum(!is.na(mag)) == 1) %>%
    # Drop bins with no significant correlations so you only get the one significant var
    drop_na()

  return(unique)
}

get_scaff_rect <- function(data){
  scaff_rect <- data %>% group_by(scaffold) %>% summarise(BIN_START = min(BIN_START), BIN_END = max(BIN_END))
  return(scaff_rect)
}

get_colors <- function(){
  vars <- c("tmean_dif", "gHM", "bio1", "NDVI", "glacier", "csi_past", "Q")
  colors <- c("#A1C65D", "#FAC723", "#936FAC", "#E95E50", "#F29222", "#0CB2AF", "#e673e0")
  names(colors) <- vars
  return(colors)
}
