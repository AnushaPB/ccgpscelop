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
#vcf <- read.vcfR(here("data/58-Sceloporus_annotated_pruned_0.6_JALMGF010000010.1.vcf"))
vcf <- read.vcfR(here("data/small_vcf.vcf"))

# only include overlapping samples
vcf <- vcf[, c(1,which(colnames(vcf@gt) %in% coords$ID))]
coords <- coords %>% filter(ID %in% colnames(vcf@gt)[-1])
coords <- coords %>% mutate(ID = factor(ID, levels = colnames(vcf@gt)[-1])) %>% arrange(ID)

# check order is identical
coords <- coords %>% arrange(factor(ID, levels =colnames(vcf@gt)[-1]))
stopifnot(all(coords$ID == colnames(vcf@gt)[-1]))

# project coords
coords_longlat <- st_as_sf(coords, coords = c("x", "y"), crs = 4326)
coords_proj <- st_transform(coords_longlat, 3310)

# RDA ----
# convert to dosage
dos <- vcf_to_dosage(vcf)

# deal with NA values (for now just removing SNPs with NA values)
#dos <- simple_impute(dos)
na <- apply(dos, 2, function(x) mean(is.na(x)))
doscc <- dos[,(na == 0)]

# extract env data and only retain complete cases
# algatr example contains PCs for CA
load_algatr_example()
env <- terra::extract(CA_env, coords[,c("x", "y")])
cci <- complete.cases(env)

# run RDA
mod_full <- rda_run(doscc[cci,], env[cci,], model = "full")

# get outliers
# note: plotting will cause a crash because there are so many snps
rda_sig_p <- rda_getoutliers(mod_full, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.0000000000000000000001, plot = FALSE)

# pull out outliers from dosage matrix 
rda_snps <- rda_sig_p$rda_snps
rda_gen <- doscc[cci, rda_snps]

# calculate correlations to identify the associated env variable
cor_df <- 
  rda_cor(rda_gen, env[cci,]) %>%
  # retain only significant correlations
  filter(p < 0.05) %>%
  # keep value that is highest for each snp
  dplyr::group_by(snp) %>%
  dplyr::filter(abs(r) == max(abs(r)))

# get number of SNPs by variable
cor_df %>% group_by(var) %>% summarize(count = n())

# Distinct maps ----

# get SNPs associated with a each env variable as a list
snps_i_ls <- 
  unique(cor_df$var) %>%
  set_names() %>%
  map(~{
    snps <- cor_df %>% filter(var == .x) %>% pull(snp)
    # get indeces for rda_gen
    which(colnames(dos) %in% snps)
  })


lyr <- coords_to_raster(coords_proj, res = 10000, buffer = 10)

# Set up parallel backend
library(doParallel)

# setup parallel session
cl <- makeCluster(8) 
registerDoParallel(cl)

# Run PC1
res_PC1 <- run_windows(snps_i_ls$CA_rPCA1, dos, coords_proj, raster(lyr))
writeRaster(res_PC1$pstk, "pstk_pc1.tif", overwrite = TRUE)
writeRaster(res_PC1$dpg, "dpg_pc1.tif", overwrite = TRUE)

# repeat for PC2 and 3:
res_PC2 <- run_windows(snps_i_ls$CA_rPCA2, dos, coords_proj, raster(lyr))
writeRaster(res_PC2$pstk, "pstk_pc2.tif", overwrite = TRUE)
writeRaster(res_PC2$dpg, "dpg_pc2.tif", overwrite = TRUE)

# repeat for PC3: 
res_PC3 <- run_windows(snps_i_ls$CA_rPCA3, dos, coords_proj, raster(lyr))
writeRaster(res_PC3$pstk, "pstk_pc3.tif", overwrite = TRUE)
writeRaster(res_PC3$dpg, "dpg_pc3.tif", overwrite = TRUE)

stopCluster(cl)

# function to calculate distinct stacks
future::plan("multisession", workers = 10)
dlstk1 <- divloss_p(res_PC1$pstk, prop = TRUE)
writeRaster(dlstk1, "dlstk_pc1.tif", overwrite = TRUE)
dlstk2 <- divloss_p(res_PC2$pstk, prop = TRUE)
writeRaster(dlstk2, "dlstk_pc2.tif", overwrite = TRUE)
dlstk3 <- divloss_p(res_PC3$pstk, prop = TRUE)
writeRaster(dlstk3, "dlstk_pc3.tif", overwrite = TRUE)
future::plan("sequential")

# create stack
dlstk <- map(1:3, ~rast(paste0("dlstk_pc", .x, ".tif")))
names(dlstk) <- c("PC1", "PC2", "PC3")

dpg <- rast(map(1:3, ~rast(paste0("dpg_pc", .x, ".tif"))[[1]]))
names(dpg) <- c("PC1", "PC2", "PC3")

dlavg <- rast(map(dlstk, ~mean(.x, na.rm = TRUE)))
dlavgscl <- rast(map(dlstk, ~mean(scale(.x), na.rm = TRUE)))

CA_env <- project(rast(CA_env), crs(ca))
dlavg <- project(dlavg, crs(ca))
dlavgscl <- project(dlavgscl, crs(ca))
dpg <- project(dpg, crs(ca))

ggdf <-
  bind_rows(as.data.frame(dlavg, xy = TRUE), data.frame(as.data.frame(dlavgscl, xy = TRUE), scale = "Scaled")) %>%
  mutate(scale = case_when(is.na(scale) ~ "Not scaled", TRUE ~ scale)) %>% 
  pivot_longer(c(PC1, PC2, PC3)) %>%
  # TRANSFORM SO THAT HIGHER VALUES = MORE DISTINCT
  mutate(value = value * -1)

plt1 <-
  ggplot(filter(ggdf, scale == "Scaled")) +
  geom_sf(data = ca) +
  geom_raster(aes(x = x, y = y, fill = value)) +
  facet_grid(scale ~ name) +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = "DI (scaled)") +
  theme_void() +
  theme(strip.text.y = element_blank())

plt2 <-
  ggplot(filter(ggdf, scale == "Not scaled")) +
  geom_sf(data = ca) +
  geom_raster(aes(x = x, y = y, fill = value)) +
  facet_grid(scale ~ name) +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = "DI") +
  theme_void() +
  theme(strip.text.y = element_blank())

plt3 <- 
  dpg %>%
  as.data.frame(xy = TRUE) %>%
  pivot_longer(c(-x, -y)) %>%
  ggplot() +
  geom_sf(data = ca) +
  geom_raster(aes(x = x, y = y, fill = value)) +
  facet_grid(. ~ name) +
  scale_fill_viridis_c(option = "inferno", na.value = NA) +
  labs(fill = "pi") +
  theme_void() +
  theme(strip.text = element_blank())

plt4 <- 
  CA_env %>%
  #mask(resample(dlavg, CA_env)) %>%
  as.data.frame(xy = TRUE) %>%
  mutate_at(c("CA_rPCA1", "CA_rPCA2", "CA_rPCA3"), scale) %>%
  pivot_longer(c(-x, -y)) %>%
  ggplot() +
  geom_sf(data = ca) +
  geom_raster(aes(x = x, y = y, fill = value)) +
  facet_grid(. ~ name) +
  scale_fill_viridis_c(option = "turbo", na.value = NA) +
  labs(fill = "Env") +
  theme_void() +
  theme(strip.text = element_blank())

gridExtra::grid.arrange(plt2, plt1, plt3, plt4, nrow = 4)
gridExtra::grid.arrange(plt2, plt3, plt4, nrow = 3)

plot(dlavg, col = viridis::inferno(100, direction = -1))
plot_gd(dlavgscl, lyr, col = viridis::inferno(100, direction = -1))

# krig results for PC1
kdlavg <- krig_gd(dlavg[[1]], index = 1, lyr, disagg_grd = 4)
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

# plot bivariate map
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

# linear model to test whether distinct diversity is predicted by genetic diversity
df <- data.frame(dlimp = values(dlimp), pg = values(pg))
mod <- lm(dlimp ~ pg, data = df)
predicted_values <- predict(mod, values(pg))

# Calculate residuals by subtracting predicted values from observed values
# this tells us where geneetic diversity does not predict distinct diversity
residuals <- your_data$dependent_variable - predicted_values
resid <- pg
resid[] <- residuals(mod)

# NO MISSING

# wingen ---

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

# SDM/resist_gd ------------
library(randomForest)

if (!file.exists(here("data", "SDM.tif"))){
  env <- stack(here(wdir, "data", "CLEANED_ENVDATA_NOCOR.tif"))
  abs <- data.frame(layer = sampleRandom(env, 100), pa = 0)
  pres <- data.frame(layer = raster::extract(env, coords), pa = 1)
  df <- rbind(pres, abs)
  df <- df[complete.cases(df),]
  names(df) <- c(names(env), "pa")

  mod <- randomForest(pa ~ ., df, mtry = 1)
  sdm <- predict(env, mod)
  plot(sdm)

  writeRaster(sdm, "data/SDM.tif")
} else {
  sdm <- rast("data/SDM.tif")
}

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



# Distinct individuals ----
coordscc <- coords[cci,]
rdapc <- prcomp(rda_gen)$x[,1:3]
rdapc_mean <- colMeans(rdapc, na.rm = TRUE)
rda_df <- 
 data.frame(coordscc, pcdist = rdapc) %>% 
 pivot_longer(c(PC1, PC2, PC3))

ggplot(rda_df) +
geom_point(aes(x = x, y = y, col = value))+
facet_wrap(~name) +
coord_sf() +
scale_color_viridis_c(option = "turbo")

rdapc <- prcomp(doscc[cci, !(colnames(doscc) %in% rda_snps)])$x[,1:3]
rdapc_mean <- colMeans(rdapc, na.rm = TRUE)
rda_df2 <- 
 data.frame(coordscc, pcdist = rdapc) %>% 
 pivot_longer(c(PC1, PC2, PC3))

ggplot(rda_df2) +
geom_point(aes(x = x, y = y, col = value))+
facet_wrap(~name) +
coord_sf() +
scale_color_viridis_c(option = "turbo")

rda_df3 <-
  bind_rows(
    data.frame(rda_df, method = "adaptive"),
    data.frame(rda_df2, method = "neutral")
  ) %>%
  mutate(method = factor(method, levels = c("neutral", "adaptive")))

ggplot(rda_df3) +
geom_sf(data = NUS_longlat) +
geom_point(aes(x = x, y = y, col = value))+
facet_grid(method~name) +
coord_sf() +
scale_color_viridis_c(option = "turbo") +
theme(
  panel.background = element_rect(color = "gray", fill = NA),
  strip.background = element_rect(color = "gray", fill = "gray"),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank()
)

gendist <- as.matrix(dist(rda_gen))
gdm <- gdm_do_everything(gendist = gendist, coords = coords[cci,c("x", "y")], envlayers = CA_env, scale_gendist = TRUE, quiet = TRUE)
plotRGB(gdm$rast$pcaRastRGB)
