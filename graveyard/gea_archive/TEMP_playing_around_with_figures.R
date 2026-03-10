library(tidyverse)
library(here)
library(terra)
library(sf)
library(MetBrewer)
source(here("general_functions.R"))
source(here("analysis", "gea", "functions_selection_stats.R"))  
outpath <- here("analysis", "gea", "outputs_contaminated_genome")
plotpath <- here("analysis", "gea", "plots")

hsp70df <- readRDS(here(outpath, "hsp70df.rds"))
hsp70loadings <- readRDS(here(outpath, "hsp70loadings.rds"))

# Kriging
lyr <- wingen::coords_to_raster(ca_proj, res = 10000)

hsp70_pc1 <- rasterize(hsp70df, lyr, field = "PC1", fun = mean, na.rm = TRUE)
hsp70_pc1 <- krig_gd2(hsp70_pc1, nmax = 30)$prediction
names(hsp70_pc1) <- "PC1"
range_map <- get_range()
hsp70_pc1 <- mask(mask(hsp70_pc1, range_map), st_transform(ca, 3310))
hsp70_pc1 <- as.data.frame(hsp70_pc1, xy = TRUE) 

hsp70_pc2 <- rasterize(hsp70df, lyr, field = "PC2", fun = mean, na.rm = TRUE)
hsp70_pc2 <- krig_gd2(hsp70_pc2, nmax = 30)$prediction
names(hsp70_pc2) <- "PC2"
hsp70_pc2 <- mask(mask(hsp70_pc2, range_map), st_transform(ca, 3310))
hsp70_pc2 <- as.data.frame(hsp70_pc2, xy = TRUE)

# Make PCA plot
plt1 <- 
  ggplot(hsp70df) +
  geom_sf(data = ca, col = NA) +
  geom_tile(data = hsp70_pc1, aes(x = x, y = y, fill = PC1)) +
  geom_sf(aes(fill = PC1), col = "black", pch = 21, cex = 1.5) +
  scale_fill_gradientn(colors=MetBrewer::met.brewer("Hiroshige", direction = 1)) +
  labs(fill = "HSP70\nPC1") +
  theme_void()  +
  theme(
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"),
    legend.position = c(0.9, 0.55),  # (x, y) in npc coords (0-1)
    legend.justification = c(1, 0), 
    legend.key.height = unit(0.6, "cm"),
    legend.key.width = unit(0.6, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

plt2 <- 
  ggplot(hsp70df) +
  geom_hline(yintercept = 0, color = "gray", linetype = "solid") +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  annotate("text", x = max(hsp70df$PC1), y = 0.1, label = "PC1", hjust = 0.5, vjust = 0, color = "darkgray") +
  annotate("text", x = 0.1, y = max(hsp70df$PC2), label = "PC2", vjust = -5, hjust = 0.1, color = "darkgray") +
  geom_point(aes(x = PC1, y = PC2, col = PC1)) +
  geom_segment(data = hsp70loadings, aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5), 
         arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  ggrepel::geom_text_repel(
    data = hsp70loadings,
    aes(x = PC1 * 5, y = PC2 * 5, label = clean_locus),
    color = "black",
    max.overlaps = Inf
  ) +
  scale_color_gradientn(colors=MetBrewer::met.brewer("Hiroshige", direction = 1)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), legend.position = "none")


plt3 <- 
  ggplot(hsp70df) +
  geom_sf(data = ca, col = NA) +
  geom_tile(data = hsp70_pc2, aes(x = x, y = y, fill = PC2)) +
  geom_sf(aes(fill = PC2), col = "black", pch = 21, cex = 1.5) +
  scale_fill_gradientn(colors=MetBrewer::met.brewer("Hiroshige", direction = -1)) +
  labs(fill = "HSP70\nPC2") +
  theme_void()  +
  theme(
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"),
    legend.position = c(0.9, 0.55),  # (x, y) in npc coords (0-1)
    legend.justification = c(1, 0), 
    legend.key.height = unit(0.6, "cm"),
    legend.key.width = unit(0.6, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

plt4 <- 
  ggplot(hsp70df) +
  geom_hline(yintercept = 0, color = "gray", linetype = "solid") +
  geom_vline(xintercept = 0, color = "gray", linetype = "solid") +
  annotate("text", x = max(hsp70df$PC1), y = 0.1, label = "PC1", hjust = 0.5, vjust = 0, color = "darkgray") +
  annotate("text", x = 0.1, y = max(hsp70df$PC2), label = "PC2", vjust = -5, hjust = 0.1, color = "darkgray") +
  geom_point(aes(x = PC1, y = PC2, col = PC2)) +
  geom_segment(data = hsp70loadings, aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5), 
         arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  ggrepel::geom_text_repel(
    data = hsp70loadings,
    aes(x = PC1 * 5, y = PC2 * 5, label = clean_locus),
    color = "black",
    max.overlaps = Inf
  ) +
  scale_color_gradientn(colors=MetBrewer::met.brewer("Hiroshige", direction = -1)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), legend.position = "none")

  
pdf(here(plotpath, "hsp70_pca.pdf"), width = 8, height = 4)
cowplot::plot_grid(plt1, plt2, nrow = 1, labels = c("A", "B"))
cowplot::plot_grid(plt3, plt4, nrow = 1, labels = c("C", "D"))
dev.off()



