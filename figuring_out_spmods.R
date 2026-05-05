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
mod_vars <- c("glacier", "tmean_dif", "gHM", "bio1", "lgm_lig", "NDVI", "Q")
f <- as.formula(paste("Ho ~", paste(mod_vars, collapse = " + ")))
ols_model <- lm(f, data = dframe)
summary(ols_model)
model_sf$residuals_ols <- residuals(ols_model)
sf_data_sp$residuals_ols <- residuals(ols_model)

# Plot residuals
ggplot(model_sf, aes(col = residuals_ols)) +
  geom_sf() +
  scale_color_viridis_c(option = "turbo") +
  theme_classic()

# Compute the semivariogram for OLS residuals
variogram_model_ols <- variogram(residuals_ols ~ 1, sf_data_sp)
vgm_model_ols <- fit.variogram(variogram_model_ols, model = vgm("Exp"))
vgm_model_ols$range # nugget is zero
range_ols <- vgm_model_ols$range[2]

# Fit the GLS model
gls_model <- gls(f, 
                 data = dframe, 
                 correlation = corExp(form = ~ X + Y, nugget = TRUE),
                 #correlation = corExp(value = range_ols, nugget = FALSE, form = ~ 1, metric = "euclidean", fixed = FALSE),
                 #correlation = corExp(value = range_value, nugget = FALSE, form = ~ 1, metric = "euclidean", fixed = FALSE),
                 method = "ML")
summary(gls_model)
model_sf$residuals_gls <- residuals(gls_model)

# Create neighbors list based on range value from semivariogram
coords <- st_coordinates(model_sf)

nb1 <- dnearneigh(coords, d1 = 0, d2 = range_value)
listw1 <- nb2listw(nb1, style = "W", zero.policy = TRUE)

nb2 <- dnearneigh(coords, d1 = 0, d2 = range_ols)
listw2 <- nb2listw(nb2, style = "W", zero.policy = TRUE)

# Fit the SEM model
sem_model1 <- errorsarlm(f, data = dframe, listw = listw1, zero.policy = TRUE)
summary(sem_model1)
model_sf$residuals_sem1 <- residuals(sem_model1)

sem_model2 <- errorsarlm(f, data = dframe, listw = listw2, zero.policy = TRUE)
summary(sem_model2)
model_sf$residuals_sem2 <- residuals(sem_model2)

# Compute the variogram of the residuals for each model
ols_variogram <- variogram(residuals_ols ~ 1, data = model_sf)
ols_vgm <- fit.variogram(ols_variogram, model = vgm("Exp"))

gls_variogram <- variogram(residuals_gls ~ 1, data = model_sf)
gls_vgm <- fit.variogram(gls_variogram, model = vgm("Exp"))

sem_variogram1 <- variogram(residuals_sem1 ~ 1, data = model_sf)
sem_vgm1 <- fit.variogram(sem_variogram1, model = vgm("Exp"))

sem_variogram2 <- variogram(residuals_sem2 ~ 1, data = model_sf)
sem_vgm2 <- fit.variogram(sem_variogram2, model = vgm("Exp"))

# Plot the variograms
plot(ols_variogram, model = ols_vgm, main = "OLS Residuals Variogram")
plot(gls_variogram, model = gls_vgm, main = "GLS Residuals Variogram")
plot(sem_variogram1, model = sem_vgm, main = "SEM 1 Residuals Variogram")
plot(sem_variogram2, model = sem_vgm2, main = "SEM 2 Residuals Variogram")

# Calculate Moran's I 
moran.test(model_sf$residuals_ols, listw1)$estimate
moran.test(model_sf$residuals_gls, listw1)$estimate
moran.test(model_sf$residuals_sem1, listw1)$estimate
moran.test(model_sf$residuals_sem2, listw1)$estimate

moran.test(model_sf$residuals_ols, listw2)$estimate
moran.test(model_sf$residuals_gls, listw2)$estimate
moran.test(model_sf$residuals_sem1, listw2)$estimate
moran.test(model_sf$residuals_sem2, listw2)$estimate

# I think SEM model 2 is the best
summary(sem_model2)
moran.test(model_sf$residuals_sem2, listw2)$estimate # Moran's I = -0.04
AIC(sem_model2)

# CONCLUSION: SEM model reduces the Moran's I to  ~0 and the residuals show the least spatial autocorrelation based on the semivariogram

# Test a range of neighborhood sizes
sem_test <- 
  map(seq(10000, range_value + 10000, 2000), ~{
    nb <- dnearneigh(coords, d1 = 0, d2 = .x)
    listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

    # Fit the SEM model
    mod_vars <- c("glacier", "tmean_dif", "gHM", "bio1", "lgm_lig", "NDVI", "Q")
    f <- as.formula(paste("Ho ~", paste(mod_vars, collapse = " + ")))
    dframe_scale <- dframe %>% mutate(across(c(Ho, glacier, tmean_dif, gHM, bio1, lgm_lig, bio1, NDVI, Q), scale))
    sem_model <- errorsarlm(f, data = dframe_scale, listw = listw, zero.policy = TRUE)
    model_sf$residuals_sem <- residuals(sem_model)

    # Calculate Moran's I of null model residuals
    ols_model <- lm(f, data = dframe_scale)
    moran_ols <- moran.test(residuals(ols_model), listw)

    # Calculate Moran's I
    moran_sem <- moran.test(model_sf$residuals_sem, listw)

    # Get model stats
    paleo_change <- broom::tidy(sem_model) %>% filter(term == "lgm_lig")
    
    data.frame(d = .x, p_paleo_change = paleo_change$p.value, beta_paleo_change = abs(paleo_change$estimate), AIC = AIC(sem_model), Moran_sem = moran_sem$estimate[1], Moran_ols = moran_ols$estimate[1])
  }, .progress = TRUE) %>%
  bind_rows()

gg_df <- 
  sem_test %>% 
  mutate(sig = case_when(p_paleo_change < 0.05 ~ "sig", TRUE ~ "not sig")) %>%
  pivot_longer(c(-d, -sig)) %>%
  mutate(name = factor(name, levels = c("Moran_sem", "Moran_ols", "AIC", "beta_paleo_change", "p_paleo_change")))

hlines <- 
  data.frame(name = c("Moran_sem", "Moran_ols", "AIC", "beta_paleo_change", "p_paleo_change"), value = c(NA, NA, NA, NA, 0.05)) %>%
  mutate(name = factor(name, levels = c("Moran_sem", "Moran_ols", "AIC", "beta_paleo_change", "p_paleo_change")))

lowest_M <- sem_test %>% filter(Moran_sem == min(Moran_sem)) %>% pull(d)
zero_M <- sem_test %>% filter(Moran_sem < 0) %>% arrange(d) %>% slice(1) %>% pull(d)
lowest_AIC <- sem_test %>% filter(AIC == min(AIC)) %>% pull(d)

plt <-
  ggplot(gg_df, aes(x = d, y = value)) +
  geom_vline(linetype = "dotted", col = "gray", xintercept = seq(20000, 100000, 20000)) +
  geom_hline(data = hlines, aes(yintercept = value, col = "alpha = 0.05"), linetype = "dotted", lwd = 1) +
  geom_vline(linetype = "dashed", aes(xintercept = lowest_M, col = "Minimum M"), lwd = 1) +
  geom_vline(linetype = "dashed", aes(xintercept = zero_M, col = "First M <= 0"), lwd = 1) +
  geom_vline(linetype = "solid", aes(xintercept = range_value, col = "Ho range"), lwd = 1) +
  geom_vline(linetype = "solid", aes(xintercept = range_ols, col = "OLS residual range"), lwd = 1) +
  geom_vline(linetype = "solid", aes(xintercept = lowest_AIC, col = "Minimum AIC"), lwd = 1) +
  geom_line() +
  geom_point(aes(fill = sig), pch = 21) +
  facet_wrap(~name, scales = "free", nrow = 2) +
  scale_x_continuous(breaks = seq(0, max(gg_df$d), by = 20000)) +
  labs(x = "Distance (m)", y = "Statistic") +
  theme_classic() +
  theme(strip.background = element_blank())

pdf("sem_test.pdf", height = 5, width = 15)
print(plt)
dev.off()


par(mfrow = c(1, 5))
plot(dnearneigh(coords, d1 = 0, d2 = 20000), coords = st_coordinates(model_sf))
title("20 km")
plot(dnearneigh(coords, d1 = 0, d2 = 40000), coords = st_coordinates(model_sf))
title("40 km")
plot(dnearneigh(coords, d1 = 0, d2 = 60000), coords = st_coordinates(model_sf))
title("60 km")
plot(dnearneigh(coords, d1 = 0, d2 = 80000), coords = st_coordinates(model_sf))
title("80 km")
plot(dnearneigh(coords, d1 = 0, d2 = 100000), coords = st_coordinates(model_sf))
title("100 km")

par(mfrow = c(1, 4))
test_sf <- model_sf 
plot(dnearneigh(test_sf, d1 = 0, d2 = 78000), coords = st_coordinates(test_sf))
title("78 km")
plot(dnearneigh(test_sf, d1 = 0, d2 = 80000), coords = st_coordinates(test_sf))
title("80 km")
plot(dnearneigh(test_sf, d1 = 0, d2 = 82000), coords = st_coordinates(test_sf))
title("82 km")
plot(dnearneigh(test_sf, d1 = 0, d2 = 84000), coords = st_coordinates(test_sf))
title("84 km")


par(mfrow = c(1, 4))
test_sf <- model_sf %>% filter(Ho > 0.0008)
nrow(test_sf) - nrow(model_sf)
plot(dnearneigh(test_sf, d1 = 0, d2 = 78000), coords = st_coordinates(test_sf))
title("78 km")
plot(dnearneigh(test_sf, d1 = 0, d2 = 80000), coords = st_coordinates(test_sf))
title("80 km")
plot(dnearneigh(test_sf, d1 = 0, d2 = 82000), coords = st_coordinates(test_sf))
title("82 km")
plot(dnearneigh(test_sf, d1 = 0, d2 = 84000), coords = st_coordinates(test_sf))
title("84 km")

# Test a range of neighborhood sizes
sem_test <- 
  map(seq(10000, range_value + 10000, 2000), ~{
    nb <- dnearneigh(test_sf, d1 = 0, d2 = .x)
    listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

    # Fit the SEM model
    dframe_scale <- test_sf %>% mutate(across(c(Ho, lgm_lig, bio1, tmean_dif, glacier, Q, NDVI), scale))
    mod_vars <- c("glacier", "tmean_dif", "gHM", "bio1", "lgm_lig", "NDVI")
    f <- as.formula(paste("Ho ~", paste(mod_vars, collapse = " + ")))
    sem_model <- errorsarlm(Ho ~ lgm_lig + Q, data = dframe_scale, listw = listw, zero.policy = TRUE)
    test_sf$residuals_sem <- residuals(sem_model)

    # Calculate Moran's I
    moran_sem <- moran.test(test_sf$residuals_sem, listw)

    # Get model stats
    paleo_change <- broom::tidy(sem_model) %>% filter(term == "lgm_lig")
    
    data.frame(d = .x, p_paleo_change = paleo_change$p.value, beta_paleo_change = abs(paleo_change$estimate), AIC = AIC(sem_model), Moran = moran_sem$estimate[1])
  }, .progress = TRUE) %>%
  bind_rows()

gg_df <- 
  sem_test %>% 
  mutate(sig = case_when(p_paleo_change < 0.05 ~ "sig", TRUE ~ "not sig")) %>%
  pivot_longer(c(-d, -sig)) %>%
  mutate(name = factor(name, levels = c("Moran", "AIC", "beta_paleo_change", "p_paleo_change")))

hlines <- 
  data.frame(name = c("Moran", "AIC", "beta_paleo_change", "p_paleo_change"), value = c(NA, NA, NA, 0.05)) %>%
  mutate(name = factor(name, levels = c("Moran", "AIC", "beta_paleo_change", "p_paleo_change")))

lowest_M <- sem_test %>% filter(Moran == min(Moran)) %>% pull(d)
zero_M <- sem_test %>% filter(Moran < 0) %>% arrange(d) %>% slice(1) %>% pull(d)
lowest_AIC <- sem_test %>% filter(AIC == min(AIC)) %>% pull(d)

new_ols_resids <- residuals(lm(Ho ~ lgm_lig, data = test_sf))
new_ols_range <- fit.variogram(variogram(new_ols_resids ~ 1, test_sf), model = vgm("Exp"))$range[2]

plt <-
  ggplot(gg_df, aes(x = d, y = value)) +
  geom_vline(linetype = "dotted", col = "gray", xintercept = seq(20000, 100000, 20000)) +
  geom_hline(data = hlines, aes(yintercept = value, col = "alpha = 0.05"), linetype = "dotted", lwd = 1) +
  geom_vline(linetype = "dashed", aes(xintercept = lowest_M, col = "Minimum M"), lwd = 1) +
  geom_vline(linetype = "dashed", aes(xintercept = zero_M, col = "First M <= 0"), lwd = 1) +
  geom_vline(linetype = "solid", aes(xintercept = range_value, col = "Ho range"), lwd = 1) +
  geom_vline(linetype = "solid", aes(xintercept = new_ols_range, col = "New OLS residual range"), lwd = 1) +
  geom_vline(linetype = "solid", aes(xintercept = lowest_AIC, col = "Minimum AIC"), lwd = 1) +
  geom_line() +
  geom_point(aes(fill = sig), pch = 21) +
  facet_wrap(~name, scales = "free", nrow = 2) +
  scale_x_continuous(breaks = seq(0, max(gg_df$d), by = 20000)) +
  labs(x = "Distance (m)", y = "Statistic") +
  theme_classic() +
  theme(strip.background = element_blank())


pdf("sem_test2.pdf", height = 5, width = 10)
print(plt)
dev.off()


ggplot(model_sf, aes(x = lgm_lig, y = Ho)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggpubr::stat_cor() 


# NO Q

# Test a range of neighborhood sizes
sem_test <- 
  map(seq(10000, range_value + 10000, 2000), ~{
    nb <- dnearneigh(coords, d1 = 0, d2 = .x)
    listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

    # Fit the SEM model
    dframe_scale <- dframe %>% mutate(across(c(Ho, glacier, tmean_dif, gHM, bio1, lgm_lig, bio1, NDVI), scale))
    # Remove Q from model formula
    mod_vars <- c("glacier", "tmean_dif", "gHM", "bio1", "lgm_lig", "NDVI")
    f <- as.formula(paste("Ho ~", paste(mod_vars, collapse = " + ")))
    sem_model <- errorsarlm(f, data = dframe_scale, listw = listw, zero.policy = TRUE)
    model_sf$residuals_sem <- residuals(sem_model)

    # Calculate Moran's I
    moran_sem <- moran.test(model_sf$residuals_sem, listw)

    # Get model stats
    paleo_change <- broom::tidy(sem_model) %>% filter(term == "lgm_lig")
    
    data.frame(d = .x, p_paleo_change = paleo_change$p.value, beta_paleo_change = abs(paleo_change$estimate), AIC = AIC(sem_model), Moran = moran_sem$estimate[1])
  }, .progress = TRUE) %>%
  bind_rows()

gg_df <- 
  sem_test %>% 
  mutate(sig = case_when(p_paleo_change < 0.05 ~ "sig", TRUE ~ "not sig")) %>%
  pivot_longer(c(-d, -sig)) %>%
  mutate(name = factor(name, levels = c("Moran", "AIC", "beta_paleo_change", "p_paleo_change")))

hlines <- 
  data.frame(name = c("Moran", "AIC", "beta_paleo_change", "p_paleo_change"), value = c(NA, NA, NA, 0.05)) %>%
  mutate(name = factor(name, levels = c("Moran", "AIC", "beta_paleo_change", "p_paleo_change")))

lowest_M <- sem_test %>% filter(Moran == min(Moran)) %>% pull(d)
zero_M <- sem_test %>% filter(Moran < 0) %>% arrange(d) %>% slice(1) %>% pull(d)
lowest_AIC <- sem_test %>% filter(AIC == min(AIC)) %>% pull(d)

plt <-
  ggplot(gg_df, aes(x = d, y = value)) +
  geom_vline(linetype = "dotted", col = "gray", xintercept = seq(20000, 100000, 20000)) +
  geom_hline(data = hlines, aes(yintercept = value, col = "alpha = 0.05"), linetype = "dotted", lwd = 1) +
  geom_vline(linetype = "dashed", aes(xintercept = lowest_M, col = "Minimum M"), lwd = 1) +
  geom_vline(linetype = "dashed", aes(xintercept = zero_M, col = "First M <= 0"), lwd = 1) +
  geom_vline(linetype = "solid", aes(xintercept = range_value, col = "Ho range"), lwd = 1) +
  geom_vline(linetype = "solid", aes(xintercept = range_ols, col = "OLS residual range"), lwd = 1) +
  geom_vline(linetype = "solid", aes(xintercept = lowest_AIC, col = "Minimum AIC"), lwd = 1) +
  geom_line() +
  geom_point(aes(fill = sig), pch = 21) +
  facet_wrap(~name, scales = "free", nrow = 2) +
  scale_x_continuous(breaks = seq(0, max(gg_df$d), by = 20000)) +
  labs(x = "Distance (m)", y = "Statistic") +
  theme_classic() +
  theme(strip.background = element_blank())

pdf("sem_test.pdf", height = 5, width = 10)
print(plt)
dev.off()

# Test a range of neighborhood sizes
sem_test <- 
  map(seq(10000, range_value + 10000, 2000), ~{
    nb <- dnearneigh(test_sf, d1 = 0, d2 = .x)
    listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

    # Fit the SEM model
    dframe_scale <- test_sf %>% mutate(across(c(Ho, lgm_lig, bio1, tmean_dif, glacier, Q), scale))
    mod_vars <- c("glacier", "tmean_dif", "gHM", "bio1", "lgm_lig")
    f <- as.formula(paste("Ho ~", paste(mod_vars, collapse = " + ")))
    sem_model <- errorsarlm(f, data = dframe_scale, listw = listw, zero.policy = TRUE)
    test_sf$residuals_sem <- residuals(sem_model)

    # Calculate Moran's I
    moran_sem <- moran.test(test_sf$residuals_sem, listw)

    # Get model stats
    paleo_change <- broom::tidy(sem_model) %>% filter(term == "lgm_lig")
    
    data.frame(d = .x, p_paleo_change = paleo_change$p.value, beta_paleo_change = abs(paleo_change$estimate), AIC = AIC(sem_model), Moran = moran_sem$estimate[1])
  }, .progress = TRUE) %>%
  bind_rows()

gg_df <- 
  sem_test %>% 
  mutate(sig = case_when(p_paleo_change < 0.05 ~ "sig", TRUE ~ "not sig")) %>%
  pivot_longer(c(-d, -sig)) %>%
  mutate(name = factor(name, levels = c("Moran", "AIC", "beta_paleo_change", "p_paleo_change")))

hlines <- 
  data.frame(name = c("Moran", "AIC", "beta_paleo_change", "p_paleo_change"), value = c(NA, NA, NA, 0.05)) %>%
  mutate(name = factor(name, levels = c("Moran", "AIC", "beta_paleo_change", "p_paleo_change")))

lowest_M <- sem_test %>% filter(Moran == min(Moran)) %>% pull(d)
zero_M <- sem_test %>% filter(Moran < 0) %>% arrange(d) %>% slice(1) %>% pull(d)
lowest_AIC <- sem_test %>% filter(AIC == min(AIC)) %>% pull(d)

new_ols_resids <- residuals(lm(Ho ~ lgm_lig, data = test_sf))
new_ols_range <- fit.variogram(variogram(new_ols_resids ~ 1, test_sf), model = vgm("Exp"))$range[2]

plt <-
  ggplot(gg_df, aes(x = d, y = value)) +
  geom_vline(linetype = "dotted", col = "gray", xintercept = seq(20000, 100000, 20000)) +
  geom_hline(data = hlines, aes(yintercept = value, col = "alpha = 0.05"), linetype = "dotted", lwd = 1) +
  geom_vline(linetype = "dashed", aes(xintercept = lowest_M, col = "Minimum M"), lwd = 1) +
  geom_vline(linetype = "dashed", aes(xintercept = zero_M, col = "First M <= 0"), lwd = 1) +
  geom_vline(linetype = "solid", aes(xintercept = range_value, col = "Ho range"), lwd = 1) +
  geom_vline(linetype = "solid", aes(xintercept = new_ols_range, col = "New OLS residual range"), lwd = 1) +
  geom_vline(linetype = "solid", aes(xintercept = lowest_AIC, col = "Minimum AIC"), lwd = 1) +
  geom_line() +
  geom_point(aes(fill = sig), pch = 21) +
  facet_wrap(~name, scales = "free", nrow = 2) +
  scale_x_continuous(breaks = seq(0, max(gg_df$d), by = 20000)) +
  labs(x = "Distance (m)", y = "Statistic") +
  theme_classic() +
  theme(strip.background = element_blank())


pdf("sem_test2.pdf", height = 5, width = 10)
print(plt)
dev.off()
