

library(terra)
library(sf)
library(here)
library(tidyverse)
source(here("general_functions.R"))
ca <- get_ca()

#Name: CHELSA V2
#Source: Karger, D.N., Conrad, O., Böhner, J., Kawohl, T., Kreft, H., Soria-Auza, R.W., Zimmermann, N.E., Linder, P., Kessler, M. (2017): Climatologies at high resolution for the Earth land surface areas. Scientific Data. 4 170122. https://doi.org/10.1038/sdata.2017.122
#Link: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2F
#Downloaded using chelsa/get_chelsa.sh script

#Information from website:
#CHELSA (Climatologies at high resolution for the earth’s land surface areas) is a very high resolution (30 arc sec, ~1km) global downscaled climate data set currently hosted by the Swiss Federal Institute for Forest, Snow and Landscape Research WSL. It is built to provide free access to high resolution climate data for research and application, and is constantly updated and refined.

# Get chelsa files
chelsa_files <- list.files(here("data", "env", "chelsa"), full.names = TRUE)
chelsa <- rast(chelsa_files)

# Crop to California
chelsa_crop <- crop(chelsa, ca)
chelsa_mask <- mask(chelsa_crop, ca)

# Export to a single file
writeRaster(chelsa_mask, here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1.tif"), overwrite = TRUE)

# Create a PCA raster from the bioclimatic layers:
library(RStoolbox)

# Perform raster PCA on just the bioclimatic data
pca <- rasterPCA(scale(chelsa_mask))

writeRaster(pca$map[[1:3]], here("data", "env", "california_chelsa_bioclim_1981-2010_V.2.1_pca.tif"), overwrite = TRUE)

# Read in LGM data
lgm <- rast(here("data", "env", "chelsa_pmip", "CHELSA_PMIP_CCSM4_BIO_01.tif"))
lgm <- crop(lgm, ca)
lgm <- mask(lgm, ca)
# Convert from Kelvin/10 to Celsius
lgm <- lgm/10 - 273.15
writeRaster(lgm, here("data", "env", "california_chelsa_pmip_ccsm4.tif"), overwrite = TRUE)

# Read in TraCE21K data
files <- list.files(here("data", "env", "chelsa_trace21k"), full.names = TRUE, pattern = "CHELSA_TraCE21k_bio01")
trace <- rast(files)
trace <- crop(trace, ca)
trace <- mask(trace, ca)

writeRaster(here("data", "env", "california_chelsa_trace21k.tif"), overwrite = TRUE)

# Calculate variance
trace <- rast(here("data", "env", "california_chelsa_trace21k.tif"))
trace_var <- app(trace, var, na.rm = TRUE)
trace_sd <- app(trace, sd, na.rm = TRUE)
names(trace_var) <- c("trace21k_var")
names(trace_sd) <- c("trace21k_sd")
trace_stack <- c(trace_var, trace_sd)

png(here("data_processing", "env", "trace21k_var.png"))
plot(trace_stack, col = viridis::turbo(100))
dev.off()

coords <- get_coords(sf = TRUE)
mod_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv"))
time_key <- data.frame(timeid = rev(seq(-200, 20, 1)), timepoint = seq(0, 22, 0.1))
time_vals <- 
  terra::extract(trace, coords, ID = FALSE) %>%
  mutate(SampleID = coords$SampleID) %>% 
  pivot_longer(-SampleID, names_to = "variable", values_to = "value") %>%
  filter(grepl("_bio01_", variable)) %>%
  mutate(timeid = gsub("_V1.0", "", gsub("CHELSA_TraCE21k_bio01_", "", variable))) %>%
  mutate(timeid = as.numeric(timeid)) %>%
  # Convert to k BP
  left_join(time_key) %>%
  # Convert to BP
  mutate(yearbp = timepoint*1000) %>%
  left_join(mod_df) %>%
  drop_na(Ho)

cumsum <- 
  time_vals %>% 
  group_by(SampleID) %>%
  arrange(yearbp) %>%
  mutate(change = abs(value - lag(value, default = first(value)))) %>%
  #mutate(change = ifelse(change > 1, change, 0)) %>%
  summarise(cumsum = sum(change, na.rm = TRUE)) %>%
  ungroup()

df <- bind_cols(mod_df) %>% left_join(coords) %>% drop_na(Ho)
coords <- coords %>% filter(SampleID %in% df$SampleID)
df$trace_sd <- terra::extract(trace_sd, coords, ID = FALSE)[,1]
df$trace_var <- terra::extract(trace_var, coords, ID = FALSE)[,1]
df <- df %>% drop_na(trace_sd, trace_var, Ho)
summary(lm(scale(Ho) ~ scale(trace_sd), data = df))
summary(lm(scale(Ho) ~ scale(trace_var), data = df))

pdf(here("data_processing", "env", "trace21k_var.pdf"), width = 4, height = 4)
gglm("trace_sd", "Ho", df)
gglm("trace_var", "Ho", df)
gglm("cumsum", "Ho", df)
ggpartial("trace_sd", "Ho", c("trace_sd", "Q"), df)
ggpartial("lgm_lig", "Ho", c("lgm_lig", "Q"), df)
gglm("lgm_lig", "Ho", df)
ggplot(time_vals, aes(x = yearbp, y = value)) +
  geom_line(aes(group = SampleID)) +
  theme_bw()
dev.off()

