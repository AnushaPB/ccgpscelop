# Read in heterozygosity for Adaptive SNPs

# Read in heterozygosity for Neutral SNPs

# OR SUBSET
# NEED TO THINK ABOUT WHETHER OR NOT TO LD PRUNE FIRST (probably should test both...)

# TEMP:
library(vcfR)
library(algatr)
library(tidyverse)
library(sf)
library(terra)
library(here)
source(here("general_functions.R"))
source(here("analysis", "genetic_diversity", "genetic_diversity.R"))

# Read in het data
callable_sites <- read.table(here("data", "ccgp_data", "58-Sceloporus_callable_sites_nsites.txt"))[1,1]

neutral <- 
  format_het(here("analysis", "genetic_diversity", "outputs", "58-Sceloporus.het"), callable_sites = callable_sites) %>%
  mutate(neutral_ho = Ho) %>%
  select(IID, neutral_ho)


gea <- 
  format_het(here("analysis", "gea", "outputs", "gea.het"), callable_sites = 1)%>%
  mutate(gea_ho = 1 - `O(HOM)`/`N(NM)`) %>%
  select(IID, gea_ho)

genes <- 
  format_het(here("analysis", "gea", "outputs", "genes.het"), callable_sites = 1) %>%
  mutate(genes_ho = 1 - `O(HOM)`/`N(NM)`) %>%
  select(IID, genes_ho)

df <- 
  left_join(neutral, gea, by = "IID") %>%
  left_join(genes, by = "IID") 

cols = c("GEA" = "cyan", "GEA genes" = "magenta")

pdf(here("analysis", "gea", "plots", "gea_neutral_adaptive_ho.pdf"), width = 6, height = 5)
cols = c("GEA" = "#00d5d1", "GEA genes" = "#e800a6")
ggplot(df) +
  geom_point(aes(x = neutral_ho, y = gea_ho, color = "GEA", fill = "GEA"), alpha = 0.5) +
  geom_point(aes(x = neutral_ho, y = genes_ho, color = "GEA genes", fill = "GEA genes"), alpha = 0.5) +
  geom_smooth(aes(x = neutral_ho, y = gea_ho, color = "GEA", fill = "GEA"), method = "lm") +
  geom_smooth(aes(x = neutral_ho, y = genes_ho, color = "GEA genes", fill = "GEA genes"), method = "lm") +
  theme_classic() +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  # Remove fill legend
  #guides(fill = FALSE) +
  labs(col = "", fill = "") +
  xlab("Neutral heterozygosity") +
  ylab("Adaptive heterozygosity") 
dev.off()

coords <- get_coords(sf = TRUE) %>% rename(IID = SampleID)
coords <- coords %>% left_join(df, by = "IID") %>% drop_na(neutral_ho)

gea_mod <- lm(gea_ho ~ neutral_ho, data = coords)
genes_mod <- lm(genes_ho ~ neutral_ho, data = coords)
coords$gea_resid <- residuals(gea_mod)
coords$genes_resid <- residuals(genes_mod)


gg_df <- 
  coords %>% 
  pivot_longer(cols = c("gea_resid", "genes_resid"), names_to = "type", values_to = "Residuals") %>%
  mutate(type = case_when(
    type == "gea_resid" ~ "GEA",
    type == "genes_resid" ~ "GEA genes"
  ))

ca <- get_ca()
plt1 <- 
  ggplot(coords) +
  geom_point(aes(x = neutral_ho, y = gea_ho, color = "GEA", fill = gea_resid), pch = 21, cex = 2) +
  geom_point(aes(x = neutral_ho, y = genes_ho, color = "GEA genes", fill = genes_resid), pch = 21, cex = 2) +
  geom_smooth(aes(x = neutral_ho, y = gea_ho, color = "GEA"), method = "lm") +
  geom_smooth(aes(x = neutral_ho, y = genes_ho, color = "GEA genes"), method = "lm") +
  theme_classic() +
  scale_color_manual(values = cols) +
  scale_fill_gradient2(low = "#1434A4", mid = "#faf0e6", high = "#EC5800", midpoint = 0) +
  # Remove fill legend
  #guides(fill = FALSE) +
  labs(col = "", fill = "Residuals") +
  xlab("Neutral heterozygosity") +
  ylab("Adaptive heterozygosity") 

plt2 <- 
  ggplot(gg_df) +
  geom_sf(data = ca) +
  geom_sf(aes(fill = Residuals, col = type), pch = 21, col = "gray20", cex = 2) + 
  # divergent color scale
  scale_fill_gradient2(low = "#1434A4", mid = "#faf0e6", high = "#EC5800", midpoint = 0) +
  scale_color_manual(values = c("#00d5d1", "#e800a6")) +
  facet_wrap(~type) +
  theme_void()

plt3a <- 
  ggplot(coords) +
    geom_sf(data = ca) +
    geom_sf(aes(col = neutral_ho)) +
    scale_color_viridis_c(option = "magma") +
    labs(col = "Neutral Ho") +
    theme_void()
plt3b <- 
  ggplot(coords) +
    geom_sf(data = ca) +
    geom_sf(aes(col = gea_ho)) +
    scale_color_viridis_c(option = "magma") +
    labs(col = "GEA Ho") +
    theme_void()
plt3c <- 
  ggplot(coords) +
    geom_sf(data = ca) +
    geom_sf(aes(col = genes_ho)) +
    scale_color_viridis_c(option = "magma") +
    labs(col = "GEA genes Ho") +
    theme_void()


library(cowplot)
# Arrange plt3a, plt3b, plt3c in a single row
plt3 <- plot_grid(plt3a, plt3b, plt3c, nrow = 1)

# Arrange plt1 and plt2 in a single row
plt12 <- plot_grid(plt1, plt2, nrow = 1, rel_widths = c(2, 3), labels = c("B", "C"))

# Arrange plt3 and plt12 in two rows
final_plot <- plot_grid(plt3, plt12, nrow = 2, labels = c("A", "", ""))


pdf(here("analysis", "gea", "plots", "gea_neutral_adaptive_resid.pdf"), width = 12, height = 7)
final_plot
dev.off()

# Spatial model 1
mod_sf <- bind_cols(coords, st_coordinates(coords))

library(spdep)
library("spatialreg")
nb_test <- map(c(10000, 50000, 100000, 200000), ~{
  nb <- dnearneigh(mod_sf, d1 = 0, d2 = .x)
  listw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  slm_model <- lagsarlm(genes_ho ~ neutral_ho, data = mod_sf, listw = listw, na.action = "na.omit", zero.policy = TRUE)
  sem_model <- errorsarlm(genes_ho ~ neutral_ho, data = mod_sf, listw = listw, na.action = "na.omit", zero.policy = TRUE)
  return(data.frame(slm = AIC(slm_model), sem = AIC(sem_model), d2 = .x))
}, .progress = TRUE) %>% bind_rows()

ggplot(nb_test) +
  geom_line(aes(x = d2, y = slm, col = "SLM")) +
  geom_line(aes(x = d2, y = sem, col = "SEM")) +
  theme_classic() +
  labs(x = "Distance", y = "AIC") 

nb_test %>% filter(slm == min(slm) | sem == min(sem)) %>% head(1)

build_spatial_models <- function(nbdist) {
  # Assuming mod_sf, mod_df, and best_model are available in the environment
  # Create neighborhood using coordinates
  nb <- dnearneigh(mod_sf, d1 = 0, d2 = nbdist)
  listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

  # Calculate Moran's I on residuals of best mod
  mod_sf$residuals <- residuals(best_model)
  morans_i <- moran.test(mod_sf$residuals, listw)
  print(morans_i)

  # Ordinary least squares model
  ols_model <- lm(genes_ho ~ neutral_ho, data = mod_sf)
  # test for evidence of autocorrelation
  lm_tests <- lm.RStests(ols_model, listw, test = "all")
  print(lm_tests)

  # Spatial lag model
  slm_model <- lagsarlm(genes_ho ~ neutral_ho, data = mod_sf, listw = listw, na.action = "na.omit", zero.policy = TRUE)

  # Spatial error model
  sem_model <- errorsarlm(genes_ho ~ neutral_ho, data = mod_sf, listw = listw, na.action = "na.omit", zero.policy = TRUE)

  # Compare AIC
  writeLines(paste("AIC OLS:", round(AIC(ols_model), 2)))
  writeLines(paste("AIC SEM:", round(AIC(sem_model), 2)))
  writeLines(paste("AIC SLM:", round(AIC(slm_model), 2)))

  # Return a list of models
  return(list(ols_model = ols_model, sem_model = sem_model, slm_model = slm_model))
}

spmod <- build_spatial_models(10000)
ols_model <- spmod$ols_model
sem_model <- spmod$sem_model
slm_model <- spmod$slm_model

# Plot residuals
resids <- 
  mod_sf %>% 
  mutate(OLS = residuals(ols_model), SEM = residuals(sem_model), SLM = residuals(slm_model)) 

plt1 <- 
  ggplot(resids) + 
  geom_sf(data = ca) +
  geom_sf(aes(fill = OLS), col = "gray20", pch = 21, cex = 2) +
  scale_fill_gradient2(low = "#1434A4", mid = "#faf0e6", high = "#EC5800", midpoint = 0) +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5), strip.text.y = element_text(angle = 90))

plt2 <- 
  ggplot(resids) + 
  geom_sf(data = ca) +
  geom_sf(aes(fill = SEM), col = "gray20", pch = 21, cex = 2) +
  scale_fill_gradient2(low = "#1434A4", mid = "#faf0e6", high = "#EC5800", midpoint = 0) +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5), strip.text.y = element_text(angle = 90))


plt3 <- 
  ggplot(resids) + 
  geom_sf(data = ca) +
  geom_sf(aes(fill = SLM), col = "gray20", pch = 21, cex = 2) +
  scale_fill_gradient2(low = "#1434A4", mid = "#faf0e6", high = "#EC5800", midpoint = 0) +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5), strip.text.y = element_text(angle = 90))


pdf(here("analysis", "gea", "plots", "spatial_model_nbdist10k.pdf"), width = 10, height = 4)
spmod1 <- plot_grid(plt1, plt2, plt3, nrow = 1)
spmod1
dev.off()

# Spatial model 2

# Construct IDW matrix
# Calculate pairwise geographic distances
geoMat <- st_distance(coords)

# Convert to a numeric matrix while preserving the original matrix structure
geoMat <- matrix(as.numeric(units::drop_units(geoMat)), nrow = nrow(geoMat), ncol = ncol(geoMat))

# Ensure the matrix is symmetrical
geoMat[upper.tri(geoMat)] <- geoMat[lower.tri(geoMat)]

# Apply inverse distance weighting
spdata.IDW <- 1 / geoMat
diag(spdata.IDW) <- 0  # Set diagonal to 0 to avoid division by zero in IDW
dummy <- rep(1, nrow(coords))
dframe <- cbind(coords, st_coordinates(coords), dummy) # Add it to the data frame

library(nlme)
null.model <- lme(fixed = genes_ho ~ neutral_ho, data = dframe, random = ~ 1 | dummy, method = "ML")
summary(null.model)

# Residuals are spatially autocorrelated
null_resids <- residuals(null.model)
null_M <- ape::Moran.I(null_resids, spdata.IDW)
null_M

exp.sp <- update(null.model, correlation = corExp(1, form = ~ X + Y), method = "ML")
writeLines("Results for exponential spatial correlation:\n")
summary(exp.sp)

# False convergence:
#gau.sp <- update(null.model, correlation = corGaus(1, form = ~ X + Y), method = "ML")
#writeLines("\nResults for Gaussian spatial correlation:\n")
#summary(gau.sp)

sph.sp <- update(null.model, correlation = corSpher(1, form = ~ X + Y), method = "ML")
writeLines("\nResults for spherical spatial correlation:\n")
summary(sph.sp)

# Get the AIC for all models
aic <- AIC(null.model, exp.sp, sph.sp)
best <- rownames(aic)[which.min(aic[, "AIC"])]
best

# the best model is exp.sp
best_model <- get(best)

summary(best_model)

# Residuals are still spatially autocorrelated
resids <- residuals(best_model)
Msp <- ape::Moran.I(resids, spdata.IDW)
Msp$observed
null_M$observed

dframe$resids <- resids
dframe$null_resids <- null_resids

spmod2 <- 
  ggplot(dframe) + 
  geom_sf(data = ca) + 
  labs(fill = "Variogram") +
  geom_sf(aes(fill = resids), col = "gray20", pch = 21, cex = 2) +
  scale_fill_gradient2(low = "#1434A4", mid = "#faf0e6", high = "#EC5800", midpoint = 0) +
  theme_void()

pdf(here("analysis", "gea", "plots", "spatial_model.pdf"), width = 12, height = 3)
plot_grid(spmod1, spmod2, nrow = 1, rel_widths = c(3, 1))
dev.off()

# predicted vs observed
coords$sppredicted <- predict(best_model)
coords$nullpredicted <- predict(null.model)
plt1 <- 
  ggplot(coords) +
  geom_point(aes(x = genes_ho, y = nullpredicted, col = "Null"),cex = 2, alpha = 0.5) +
  geom_point(aes(x = genes_ho, y = sppredicted, col = "Spatial"), cex = 2, alpha = 0.5) +
  geom_smooth(aes(x = genes_ho, y = sppredicted, col = "Spatial"), method = "lm") +
  geom_smooth(aes(x = genes_ho, y = nullpredicted, col = "Null"), method = "lm") +
  geom_abline(lty = "dashed") +
  labs(fill = "", col = "") +
  theme_classic() +
  xlab("Observed") +
  ylab("Predicted")

plt2 <-  
  ggplot(dframe) + 
  geom_sf(data = ca) + 
  labs(fill = "Spatial residuals") +
  geom_sf(aes(fill = resids), col = "gray20", pch = 21, cex = 2) +
  scale_fill_gradient2(low = "#1434A4", mid = "#faf0e6", high = "#EC5800", midpoint = 0) +
  theme_void()
plt3 <-  
  ggplot(dframe) + 
  geom_sf(data = ca) + 
  labs(fill = "Null residuals") +
  geom_sf(aes(fill = null_resids), col = "gray20", pch = 21, cex = 2) +
  scale_fill_gradient2(low = "#1434A4", mid = "#faf0e6", high = "#EC5800", midpoint = 0) +
  theme_void()


pdf(here("analysis", "gea", "plots", "predicted_vs_observed.pdf"), width = 13, height = 4)
plot_grid(plt1, plt2, plt3, nrow = 1, rel_widths = c(5, 4, 4))
dev.off()

# PCA
pca_genes <- read.table(here("analysis", "gea", "outputs", "genes.eigenvec"))