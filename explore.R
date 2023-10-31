library(here)
library(tidyverse)
library(vcfR)
library(wingen)
library(sf)
library(algatr)
library(gdl)
library(terra)
library(tigris)
library(raster)

# Load the U.S. state boundaries data
states <- states(cb = TRUE)
# Extract the boundary of California (CA)
ca <- states[states$STUSPS == "CA", "STUSPS"]

# sample coords
coords <- read_table(here("data/58-Sceloporus.coords.txt"), col_names = FALSE)
colnames(coords) <- c("ID", "x", "y")

# WGS data
#vcfbig <- read.vcfR(here("58-Sceloporus/CCGP/58-Sceloporus_annotated_pruned_0.6.vcf.gz"))
vcf <- read.vcfR(here("data/58-Sceloporus_annotated_pruned_0.6_JALMGF010000010.1.vcf"))

# only include overlapping samples
vcf <- vcf[, c(1,which(colnames(vcf@gt) %in% coords$ID))]
coords <- coords %>% filter(ID %in% colnames(vcf@gt)[-1])

# check order is identical
stopifnot(all(coords$ID == colnames(vcf@gt)[-1]))

# project coords
coords_longlat <- st_as_sf(coords, coords = c("x", "y"), crs = 4326)
coords_proj <- st_transform(coords_longlat, 3085)
# Create layer for wingen
lyr <- coords_to_raster(coords_proj, res = 10000, buffer = 10)

# Run wingen
pg <- window_gd(vcf, coords_proj, lyr, stat = "pi", wdim = 11, fact = 0, rarify = FALSE)
plot_gd(pg)
pg_longlat <- terra::project(pg, "+proj=longlat")

#download states from tigris
states <- states(cb = TRUE)

# reproject into wgs84 to match coordinates
states_longlat <- st_transform(states, crs = "+proj=longlat")
states_proj <- st_transform(states, crs = 3085)

# subset out CONUS
conus <- states_longlat[-which(states_longlat$NAME %in% c("Alaska", "Hawaii", "Puerto Rico", "American Samoa", "Guam", "Commonwealth of the Northern Mariana Islands", "United States Virgin Islands")), "STUSPS"]
# convert to SPDF for plotting
conus <- as_Spatial(conus)

# subset out Northern US
NUS_proj <- states_proj[which(states_proj$NAME %in% c("California")), "STUSPS"]
NUS_longlat <- states_longlat[which(states_longlat$NAME %in% c("California")), "STUSPS"]

kg <- krig_gd(pg, index = 1, lyr, disagg_grd = 4)

mg <-
  terra::project(kg, terra::crs(NUS_longlat)) %>%
  terra::mask(NUS_longlat) %>%
  trim()


# See the SpatialKDE package for more options and details about using KDE
k <- kde(
  coords_proj,
  kernel = "quartic",
  band_width = 100000,
  decay = 1,
  grid = raster(kg),
)

par(pty = "s")
# Visualize KDE layer
plot_count(k, main = "KDE")

mg2 <- mask_gd(kg[[1]], terra::project(rast(k), terra::crs(kg)), minval = 1)

library(randomForest)

env <- raster::stack(here("data", "CLEANED_ENVDATA_NOCOR.tif"))
abs <- data.frame(layer = sampleRandom(env, 100), pa = 0)
pres <- data.frame(layer = raster::extract(env, coords), pa = 1)
df <- rbind(pres, abs)
df <- df[complete.cases(df),]
names(df) <- c(names(env), "pa")

mod <- randomForest(pa ~ ., df, mtry = 1)
sdm <- predict(env, mod)
plot(sdm)

writeRaster(sdm, "data/SDM.tif")

sdm <- project(sdm, crs(coords_proj))
sdm20 <- terra::aggregate(sdm, 30, na.rm = TRUE)
distmat <- get_resdist(coords, lyr = sdm30, ncores= 20)
#write.csv(distmat, "distmat30.csv")

rg <- resist_gd(vcf, coords_proj, lyr = sdm30, distmat = distmat,
                maxdist = quantile(distmat, 0.01, na.rm = TRUE))


# GDM -----

load_algatr_example()
gendist_id <- read_table("58-Sceloporus/CCGP/58-Sceloporus_annotated_pruned_0.6.dist.id", col_names = FALSE)
gendist <- read_table("58-Sceloporus/CCGP/58-Sceloporus_annotated_pruned_0.6.dist", col_names = FALSE)
gendist <- as.matrix(gendist)
gendist_id <- as.matrix(gendist_id)
row.names(gendist) <- colnames(gendist) <- gendist_id[,2]
gendist <- gendist[coords$ID, coords$ID]

gdm <- gdm_do_everything(gendist = gendist, coords = coords[,c("x", "y")], envlayers = CA_env)

# RDA ----
dos <- vcf_to_dosage(vcf)
#dos <- simple_impute(dos)
na <- apply(dos, 2, function(x) mean(is.na(x)))
doscc <- dos[,(na == 0)]
env <- terra::extract(CA_env, coords[,c("x", "y")])
cci <- complete.cases(env)
mod_full <- rda_run(doscc[cci,], env[cci,], model = "full")
rda_sig_p <- rda_getoutliers(mod_full, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.0000000000000000000001, plot = FALSE)
rda_snps <- rda_sig_p$rda_snps
rda_gen <- doscc[cci, rda_snps]
cor_df <- rda_cor(rda_gen, env[cci,])
cor_df <- cor_df %>%
  filter(p < 0.05) %>%
  dplyr::group_by(snp) %>%
  dplyr::filter(abs(r) == max(abs(r)))
cor_df %>% group_by(var) %>% summarize(count = n())

# Distinct maps ----
PC_snps <- cor_df %>% filter(var == "CA_rPCA1") %>% pull(snp)
rda_snps_i <- which(names(rda_sig_p$pvalues) %in% PC_snps)
future::plan("multisession", workers = 8)
#safe_window_p <- safely(window_p)
pstk1 <- window_p(vcf[rda_snps_i,], coords_proj, lyr, wdim = 11, fact = 0)
dpg1 <- window_gd(vcf[rda_snps_i,], coords_proj, lyr, stat = "pi", wdim = 11, fact = 0, rarify = FALSE)
future::plan("sequential")

PC_snps <- cor_df %>% filter(var == "CA_rPCA2") %>% pull(snp)
rda_snps_i <- which(names(rda_sig_p$pvalues) %in% PC_snps)
future::plan("multisession", workers = 10)
#safe_window_p <- safely(window_p)
pstk2 <- window_p(vcf[rda_snps_i,], coords_proj, lyr, wdim = 11, fact = 0, parallel_option = 2)
pg2 <- window_gd(vcf[rda_snps_i,], coords_proj, lyr, stat = "pi", wdim = 11, fact = 0, rarify = FALSE)
future::plan("sequential")


PC_snps <- cor_df %>% filter(var == "CA_rPCA3") %>% pull(snp)
rda_snps_i <- which(names(rda_sig_p$pvalues) %in% PC_snps)
future::plan("multisession", workers = 10)
#safe_window_p <- safely(window_p)
pstk3 <- window_p(vcf[rda_snps_i,], coords_proj, lyr, wdim = 11, fact = 0, parallel_option = 2)
pg3 <- window_gd(vcf[rda_snps_i,], coords_proj, lyr, stat = "pi", wdim = 11, fact = 0, rarify = FALSE)
future::plan("sequential")


dlstk <- divloss_p(pstk1, prop = TRUE)
dlavg <- mean(dlstk, na.rm = TRUE)
dlavgscl <- mean(scale(dlstk), na.rm = TRUE)
plot_gd(dlavg, lyr, col = viridis::inferno(100, direction = -1))
plot_gd(dlavgscl, lyr, col = viridis::inferno(100, direction = -1))

kdlavg <- krig_gd(dlavg, index = 1, lyr, disagg_grd = 4)
crs(kdlavg) <- terra::crs(NUS_proj)
mdlavg <-
  terra::project(kdlavg, terra::crs(NUS_longlat)) %>%
  terra::mask(NUS_longlat) %>%
  trim()

kg <- krig_gd(pg, index = 1, lyr, disagg_grd = 4)
crs(kg) <- terra::crs(NUS_proj)
mg <-
  terra::project(kg, terra::crs(NUS_longlat)) %>%
  terra::mask(NUS_longlat) %>%
  trim()

dlimp <- abs(dlavg)
col.matrix <- colmat(nquantiles = 5, xlab = "Distinct Diversity", ylab = " Genetic Diversity")
bivmap <- bivariate.map(rasterx = raster(dlimp),
                        rastery = raster(pg1),
                        colormatrix = col.matrix,
                        nquantiles = 5)
plot_gd(bivmap, col = as.vector(col.matrix), legend = FALSE, main = "Bivariate Map")



col.matrix <- colmat(nquantiles = 2, xlab = "Distinct Diversity", ylab = " Genetic Diversity")
bivmap <- bivariate.map(rasterx = raster(dlimp),
                        rastery = raster(mg),
                        colormatrix = col.matrix,
                        nquantiles = 2)
plot_gd(bivmap, col = as.vector(col.matrix), legend = FALSE, main = "Bivariate Map")

df <- data.frame(dlimp = values(dlimp), pg = values(pg))
mod <- lm(dlimp ~ pg, data = df)
predicted_values <- predict(mod, values(pg))

# Calculate residuals by subtracting predicted values from observed values
residuals <- your_data$dependent_variable - predicted_values

resid <- pg
resid[] <- residuals(mod)

# NO MISSING

# NICHE DISTANCE
env <- rast(here("data", "CLEANED_ENVDATA_NOCOR.tif"))
env <- terra::project(env, crs(coords_proj))
plot(env[[1]], col = viridis::mako(100))

df <- extract(env, coords_proj, xy = TRUE)
df <- df[,c("x", "y", "CA_rPCA1", "CA_rPCA2", "CA_rPCA3")]
center <- df %>% summarize_at(c("CA_rPCA1", "CA_rPCA2", "CA_rPCA3"), mean, na.rm = TRUE)
df <-
  df %>%
  mutate(envdist = sqrt((CA_rPCA1 - center$CA_rPCA1)^2 + (CA_rPCA2 - center$CA_rPCA2)^2 + (CA_rPCA3 - center$CA_rPCA3)^2),
         env1dist = sqrt((CA_rPCA1 - center$CA_rPCA1)^2))

env_df <- as.data.frame(env[[1:3]], xy = TRUE)
env_df <-
  env_df %>%
  mutate(envdist = sqrt((CA_rPCA1 - center$CA_rPCA1)^2 + (CA_rPCA2 - center$CA_rPCA2)^2 + (CA_rPCA3 - center$CA_rPCA3)^2),
         env1dist = sqrt((CA_rPCA1 - center$CA_rPCA1)^2)) %>%
  mutate(envdist = case_when(envdist > max(df$envdist, na.rm = TRUE) ~ NA, TRUE ~ envdist),
         env1dist = case_when(env1dist > max(df$env1dist, na.rm = TRUE) ~ NA, TRUE ~ env1dist))

df_sf <- df %>% st_as_sf(coords = c("x", "y"), crs = st_crs(coords_proj))

ggplot() +
  geom_raster(data = env_df, aes(x = x, y = y, fill = envdist)) +
  geom_sf(data = df_sf, aes(fill = envdist), cex = 3, pch = 21, col = "black") +
  scale_color_viridis_c(option = "plasma") +
  scale_fill_viridis_c(option = "plasma") +
  theme_void()


lyr <- coords_to_raster(coords_proj, res = 10000, buffer = 10)
winenv <- window_general(df$CA_rPCA1, coords = coords_proj, lyr = lyr, stat = mean, na.rm = TRUE, wdim = 11, fact = 0)
winenvvar <- window_general(df$CA_rPCA1, coords = coords_proj, lyr = lyr, stat = var, na.rm = TRUE, wdim = 11, fact = 0)
winenvdist <- window_general(df$env1dist, coords = coords_proj, lyr = lyr, stat = mean, na.rm = TRUE, wdim = 11, fact = 0)

ggplot_gd(winenv, bkg = NUS_proj) + ggtitle("PCA 1")
ggplot_gd(winenvvar, bkg = NUS_proj) + ggtitle("PCA 1 (Variance)")
ggplot_gd(log(winenvdist), bkg = NUS_proj) + ggtitle("PCA 1 log(Distance from Mean)")

