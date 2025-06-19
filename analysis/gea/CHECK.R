fut = envlayers$env_fut
fut = fut[["CHELSA_bio1_2071.2100_gfdl.esm4_ssp585_V.2.1"]]
cur = envlayers$env_pres$BIO1
pdf(here("analysis", "gea", "plots", "bio1_change.pdf"), width = 8, height = 5)

plot(fut, col = viridis::viridis(100), main = "Future")
plot(cur, col = viridis::viridis(100), main = "Current")
plot(fut - cur, col = viridis::viridis(100), main = "Change")
plot(offset$Proj_offset_RCP26$RDA1, col = viridis::viridis(100), main = "Offset 26")
plot(offset_85, col = viridis::viridis(100), main = "Offset 26")
plot(sum(offset_85), col = viridis::viridis(100), main = "Offset 26")
dev.off()


het <- get_het()

change = fut - cur
writeRaster(change, here("analysis", "gea", "outputs", "bio1_change.tif"), overwrite = TRUE)

csi <- rast(here("data", "env", "csi", "csi.tif"))

het <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv"))
coords <- get_coords(sf = TRUE)
coords$change <- terra::extract(change, coords)
coords$offset <- terra::extract(sum(offset_85), coords)
coords$past_csi <- terra::extract(csi[["Past"]], coords, ID = FALSE)[,1]
coords$csi <- terra::extract(csi[["SSP 3-7.0"]], coords, ID = FALSE)[,1]
df <- left_join(coords, het, by = "SampleID")
library(ggpubr)


pdf(here("analysis", "gea", "plots", "climate_change_vs_het.pdf"), width = 4, height = 4)
gglm("past_csi", "Ho", df) + xlab("Past climate stability")
gglm("csi", "Ho", df) + xlab("Future climate stability")


ggplot(df, aes(x = Ho, y = change)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Future climate change", y = "Genome-wide heterozygosity") +
  theme_classic() +
  stat_cor()

ggplot(df, aes(x = Ho, y = offset)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Future climate change", y = "Genome-wide heterozygosity") +
  theme_classic() +
  stat_cor()
dev.off()



pres <- rast(envlayers$env_pres)



# Assuming 'pres' is a raster stack or brick
pres <- stack(envlayers$env_pres)  # Create a raster stack

# Perform PCA using RasterPCA function
pca_result <- rasterPCA(pres)

# Extract the first two principal components (PCs)
pres <- pca_result$map[[1:2]]

rgb <- map(1:nlyr(pres), function(i) {
  layer <- pres[[i]]
  stats <- global(layer, fun = range, na.rm = TRUE)
  vmin <- stats[1, 1]
  vmax <- stats[1, 2]
  ((layer - vmin) / (vmax - vmin)) * 255
}) %>% rast()
rgb <- c(rgb, pres[[1]] * 0)
names(rgb) <- c("bio1", "bio2", "bio3")


pdf("bio1_ndvi_rgb.pdf", width = 8, height = 5)
plotRGB(rgb, r = 2, g = 3, b = 1)
dev.off()
