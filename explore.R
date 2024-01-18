library(here)
library(tidyverse)
library(vcfR)
library(wingen)
library(sf)
library(algatr)
library(gdl)
library(terra)
library(raster)

# sample coords
coords <- read_table(here("data/58-Sceloporus.coords.txt"), col_names = FALSE)
colnames(coords) <- c("ID", "x", "y")

# WGS data
#vcfbig <- read.vcfR(here("58-Sceloporus/CCGP/58-Sceloporus_annotated_pruned_0.6.vcf.gz"))
vcf <- read.vcfR(here("data/58-Sceloporus_annotated_pruned_0.6_JALMGF010000010.1.vcf"))

# only include overlapping samples
vcf <- vcf[, c(1,which(colnames(vcf@gt) %in% coords$ID))]
coords <- coords %>% filter(ID %in% colnames(vcf@gt)[-1])
coords <- coords %>% mutate(ID = factor(ID, levels = colnames(vcf@gt)[-1])) %>% arrange(ID)

# check order is identical
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
cl <- makeCluster(8) 

# setup parallel session

# Run PC1
registerDoParallel(cl)
res_PC1 <- run_windows(snps_i_ls$CA_rPCA1, dos, coords_proj, raster(lyr))
stopCluster(cl)

# repeat for PC2 and 3:
future::plan("multisession", workers = 10)
res_PC2 <- run_windows(snps_i_ls$CA_rPCA2, dos, coords_proj, lyr)
future::plan("sequential")

# repeat for PC3: 
future::plan("multisession", workers = 10)
res_PC1 <- run_windows(snps_i_ls$CA_rPCA3, dos, coords_proj, lyr)
future::plan("sequential")

# function to calculate distinct stacks
dlstk <- divloss_p(pstk1, prop = TRUE)
dlavg <- mean(dlstk, na.rm = TRUE)
dlavgscl <- mean(scale(dlstk), na.rm = TRUE)
plot_gd(dlavg, lyr, col = viridis::inferno(100, direction = -1))
plot_gd(dlavgscl, lyr, col = viridis::inferno(100, direction = -1))

# krig results
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

