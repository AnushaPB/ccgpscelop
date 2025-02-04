library(spdep)
library(sf)
library(gstat)
library(nlme)
library(here)
library(tidyverse)
library(spdep)
library(spatialreg)

model_df <- read_csv(here("analysis/genetic_diversity/outputs/model_df.csv"))

# Convert to sf and sp objects (sp needed for semivariograms)
model_sf <- st_as_sf(model_df, coords = c("x", "y"), crs = 4326) %>% st_transform(3310)
sf_data_sp <- as(model_sf, "Spatial")

# Compute the semivariogram for heterozygosity
variogram_model <- variogram(Ho ~ 1, sf_data_sp)
vgm_model <- fit.variogram(variogram_model, model = vgm("Exp"))

# Plot semivariogram
png("delete.png")
plot(variogram_model, model = vgm_model)
dev.off()

# Get the range value (the distance at which the semivariogram reaches its sill)
range_value <- vgm_model$range[2]
print(range_value)

# Make model df with coordinates
dframe <- bind_cols(model_df, st_coordinates(model_sf))

# Fit the OLS model
ols_model <- lm(Ho ~ cur_lgm + bio1 + tmean_dif, data = dframe)
summary(ols_model)
model_sf$residuals_ols <- residuals(ols_model)

# Fit the GLS model
gls_model <- gls(Ho ~ cur_lgm + bio1 + tmean_dif, data = dframe, correlation = corExp(form = ~x+y), method = "ML")
summary(gls_model)
model_sf$residuals_gls <- residuals(gls_model)

# Create neighbors list based on range value from semivariogram
coords <- st_coordinates(model_sf)
nb <- dnearneigh(coords, d1 = 0, d2 = range_value)
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
plot(nb, coords = st_coordinates(model_sf))

# Fit the SEM model
sem_model <- errorsarlm(Ho ~ cur_lgm + bio1 + tmean_dif, data = dframe, listw = listw, zero.policy = TRUE)
summary(sem_model)
model_sf$residuals_sem <- residuals(sem_model)

# Compute the variogram of the residuals for each model
ols_variogram <- variogram(residuals_ols ~ 1, data = model_sf)
ols_vgm <- fit.variogram(ols_variogram, model = vgm("Exp"))

gls_variogram <- variogram(residuals_gls ~ 1, data = model_sf)
gls_vgm <- fit.variogram(gls_variogram, model = vgm("Exp"))

sem_variogram <- variogram(residuals_sem ~ 1, data = model_sf)
sem_vgm <- fit.variogram(sem_variogram, model = vgm("Exp"))

# Plot the variograms
plot(ols_variogram, model = ols_vgm, main = "OLS Residuals Variogram")
plot(gls_variogram, model = gls_vgm, main = "GLS Residuals Variogram")
plot(sem_variogram, model = sem_vgm, main = "SEM Residuals Variogram")


# Calculate Moran's I 
moran.test(model_sf$residuals_ols, listw)$estimate
moran.test(model_sf$residuals_gls, listw)$estimate
moran.test(model_sf$residuals_sem, listw)$estimate

# CONCLUSION: SEM model reduces the Moran's I to  ~0 and the residuals show the least spatial autocorrelation based on the semivariogram

# Test a range of neighborhood sizes
sem_test <- 
  map(seq(10000, 100000, 5000), ~{
    nb <- dnearneigh(coords, d1 = 0, d2 = .x)
    listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

    # Fit the SEM model
    sem_model <- errorsarlm(Ho ~ cur_lgm + bio1 + tmean_dif, data = dframe, listw = listw, zero.policy = TRUE)
    model_sf$residuals_sem <- residuals(sem_model)

    # Calculate Moran's I
    moran_sem <- moran.test(model_sf$residuals_sem, listw)

    # Get model stats
    cur_lgm <- broom::tidy(sem_model) %>% filter(term == "cur_lgm")
    
    data.frame(p_cur_lgm = cur_lgm$p.value, beta_cur_lgm = abs(cur_lgm$estimate), d = .x, Moran = moran_sem$estimate[1])
  }, .progress = TRUE) %>%
  bind_rows()

gg_df <- sem_test %>% pivot_longer(-d)
hlines <- data.frame(name = c("Moran", "beta_cur_lgm", "p_cur_lgm"), value = c(NA, NA, 0.05))
lowest_M <- sem_test %>% filter(Moran == min(Moran)) %>% pull(d)
zero_M <- sem_test %>% filter(Moran < 0) %>% filter(Moran == max(Moran)) %>% pull(d)

ggplot(gg_df, aes(x = d, y = value)) +
  geom_point() +
  geom_line() +
  geom_hline(data = hlines, aes(yintercept = value, col = "alpha = 0.05"), linetype = "dashed") +
  geom_vline(linetype = "dashed", aes(xintercept = lowest_M, col = "Minimum M")) +
  geom_vline(linetype = "dashed", aes(xintercept = zero_M, col = "First M <= 0")) +
  facet_wrap(~name, scales = "free", nrow = 1) +
  labs(x = "Distance (m)", y = "Statistic") 
