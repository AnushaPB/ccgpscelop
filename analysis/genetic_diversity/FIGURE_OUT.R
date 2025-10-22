
mod_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv"))
mod_sf <- 
  get_coords(sf = TRUE) %>% 
  right_join(mod_df) %>% 
  st_transform(3310)

# GLS can't handle 0 distances, so take the average of the coordinates
mod_sf_unique <- 
  mod_sf %>%
  group_by(x, y) %>%
  summarise_at(c(mod_vars, "Ho"), mean, na.rm = TRUE) %>%
  ungroup()

# Scale data
mod_sf_unique_scaled <-
  mod_sf_unique %>%
  mutate(
    across(
      all_of(setdiff(mod_vars, "glacier")),
      ~ as.numeric(scale(.x))   # <- drop matrix
    )
  )

# Confirm that all glacier values are 0 and 1
stopifnot(all(mod_sf_unique_scaled$glacier %in% c(0, 1)))

# Convert glacier to factor
mod_sf_unique_scaled$glacier_factor <- as.factor(mod_sf_unique$glacier)
mod_sf_unique$glacier_factor <- as.factor(mod_sf_unique$glacier)

f <- formula("Ho ~ NDVI + gHM +  csi_past + Q + tmean_dif + bio1")
final_model <- gls(f, correlation = corExp(form = ~ x + y, nugget = FALSE), data = mod_sf_unique_scaled)
summary(final_model)
