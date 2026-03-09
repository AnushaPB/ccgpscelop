library(tidyverse)
library(here)
library(sf)
library(terra)
library(wingen)
source(here("general_functions.R"))
ca <- get_ca() %>% st_transform(3310)
df <- 
  read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv")) %>%
  st_as_sf(coords = c("x", "y"), crs = 4326) %>%
  st_transform(3310)

lyr <- coords_to_raster(df, res = 10000, buffer = 20)
lyr <- crop(lyr, ca)

png(here("TEST.png"))
plot(lyr)
lines(ca)
dev.off()

custom_mean <- function(x) mean(x, na.rm = TRUE)
wgd <- window_general(df$Ho, df, lyr, custom_mean, wdim = 5)

#-------------------------------
# helper: raster -> point df
#-------------------------------
.rast_to_pts <- function(r) {
  pts <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(pts)[3] <- "value"
  pts
}

#-------------------------------
# helper: prediction grid
#-------------------------------
.make_pred_df <- function(grd) {
  terra::as.data.frame(grd, xy = TRUE, na.rm = FALSE)[, c("x", "y")]
}

#-------------------------------
# helper: vector predictions -> raster
#-------------------------------
.vec_to_rast <- function(pred, grd, name = "pred") {
  out_df <- .make_pred_df(grd)
  out_df[[name]] <- pred
  out <- terra::rast(grd)
  out[] <- out_df[[name]]
  names(out) <- name
  out
}

smooth_gam <- function(r, grd = NULL, k = 100) {
  if (is.null(grd)) grd <- r
  
  pts <- .rast_to_pts(r)
  
  fit <- mgcv::gam(
    value ~ s(x, y, bs = "tp", k = k),
    data = pts,
    method = "REML"
  )
  
  pred_xy <- .make_pred_df(grd)
  pred <- predict(fit, newdata = pred_xy)
  
  .vec_to_rast(pred, grd, name = names(r))
}

sgd <- smooth_gam(wgd[[1]], k = 100)
kgd <- wkrig_gd(wgd[[1]], weight_r = wgd[[2]])
kgd_bad <- wkrig_gd(wgd[[1]])
rangemap <- get_range()
mgd <- mask(sgd, rangemap) 
mkgd <- mask(kgd, rangemap)

png(here("TEST.png"), bg = "transparent")
ggplot_gd(mkgd, bkg = ca) + theme(
  panel.background = element_rect(fill = "transparent", color = NA),
  plot.background = element_rect(fill = "transparent", color = NA),
  legend.position = "none"
) 
dev.off()

png(here("TEST.png"), bg = "transparent")
ggplot_gd(mgd, bkg = ca, col = viridis::plasma(100)) + theme(
  panel.background = element_rect(fill = "transparent", color = NA),
  plot.background = element_rect(fill = "transparent", color = NA),
  legend.position = "right"
) 
dev.off()


sgd <- smooth_gam(wgd[[1]], k = 100)
highlyr <- disagg(lyr, 10)
sgd <- smooth_gam(wgd[[1]], highlyr, k = 20)
png(here("TEST.png"), width = 30, height = 40, units = "in", res = 300)
ggplot_gd(sgd) + theme_void() + theme(legend.position = "none")
dev.off()


png(here("TEST.png"), width = 3.5, height = 3, units = "in", res = 300)
ggplot(df, aes(x = csi_past, y = Ho)) +
  geom_point(aes(col = Ho), size = 1) +
  # Geom smooth with quadratic formula
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "black") + 
  theme_classic() +
  xlab("Paleoclimate stability") +
  ylab("Genetic Diversity") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  scale_color_viridis_c(option = "plasma")
dev.off()

df2 <- 
  df %>%
  mutate(
    fire_frq = gsub(" Frequency", "", fire_frq),
    fire_frq = factor(fire_frq, levels = c("No Fire/Low", "Intermediate", "High")),
    fire_recent = factor(fire_recent, levels = c("Burned", "Unburned"))
  )

png(here("TEST.png"), width = 4.5, height = 3.5, units = "in", res = 300)
ggplot(df2, aes(x = fire_frq, y = Ho)) +
  geom_boxplot(aes(fill = fire_recent), position = position_dodge(width = 0.8)) +
  geom_point(
    aes(color = Ho, fill = fire_recent),
    position = position_jitterdodge(
      jitter.width = 0.5,
      jitter.height = 0,
      dodge.width = 0.8
    ),
    size = 1
  ) +
  theme_classic() +
  xlab("Historical fire frequency") +
  ylab("Genetic Diversity") +
  labs(fill = "Recent fire") +
  # Remove color legend
  guides(color = "none") +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.18, 0.8),
    legend.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  ) +
  scale_color_viridis_c(option = "plasma") +
  scale_fill_manual(values = c("darkgray", "white")) 
dev.off()


# Just plot plasma legend only
png(here("TEST.png"), width = 3, height = 3, units = "in", res = 300)
ggplot(df, aes(x = csi_past, y = Ho)) +
  geom_point(aes(col = Ho), size = 1) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  ) +
  scale_color_viridis_c(option = "plasma", name = "Genetic Diversity")
dev.off()


