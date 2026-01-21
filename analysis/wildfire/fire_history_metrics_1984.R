library(tidyverse)
library(terra)
library(here)
library(sf)

files <- list.files("analysis/wildfire/conus_fire_history_metrics_1984_2022", full.names = TRUE, pattern = ".tif$")
fire_stack <- rast(files)

fire_values <- extract(fire_stack, mod_sf_unique_scaled, ID = FALSE)
mod_sf_unique_scaled$fire_recent <- scale(log1p(fire_values$conus_1984_2022_FRQ))
sum(is.na(mod_sf_unique_scaled$fire_count)) # 0

gls(Ho ~ bio1 + csi_past + tmean_dif + gHM + NDVI + glacier + fire_severity + groupveg + vdep + fire_recent + Q, data = mod_sf_unique_scaled,
    correlation = corExp(form = ~ x + y, nugget = FALSE)) %>% broom::tidy() %>% filter(p.value < 0.1) %>% mutate(p.value = round(p.value, 3)) %>% arrange(desc(abs(estimate)))

gls(Ho ~ bio1 + csi_past + tmean_dif + gHM + NDVI + glacier + fire_severity + vdep + groupveg + fire_recent,
    data = mod_sf_unique_scaled,
    correlation = corExp(form = ~ x + y, nugget = FALSE)) %>% broom::tidy() %>% filter(p.value < 0.1) %>% mutate(p.value = round(p.value, 3)) %>% arrange(desc(abs(estimate)))

lm(Q ~ bio1 + csi_past + tmean_dif + gHM + NDVI + glacier + fire_severity + vdep + groupveg + fire_recent,
    data = mod_sf_unique_scaled) %>% broom::tidy() %>% filter(p.value < 0.1) %>% mutate(p.value = round(p.value, 3)) %>% arrange(desc(abs(estimate)))


mod_Q <- 
  gls(Ho ~ bio1 + csi_past + tmean_dif + gHM + NDVI + glacier + fire_severity + groupveg + vdep + fire_recent + Q, data = mod_sf_unique_scaled, method = "ML", 
  correlation = corExp(form = ~ x + y, nugget = FALSE)) 

mod_noQ <- 
  gls(Ho ~ bio1 + csi_past + tmean_dif + gHM + NDVI + glacier + fire_severity + vdep + groupveg + fire_recent,
  data = mod_sf_unique_scaled, method = "ML", 
  correlation = corExp(form = ~ x + y, nugget = FALSE)) 


# Multicolinearity check with performance package 
performance::check_collinearity(mod_Q) # high correlation for fire_severity and group_veg (VIF > 3 for both)
performance::check_collinearity(mod_noQ) # VIF for fire severity moderate (2.36) and for groupveg it is very high (4.13)

# Post hoc test for groupveg
library(emmeans)
emm <- emmeans::emmeans(mod_noQ, ~ groupveg, mode = "appx-satterthwaite")
pairs(emm) %>% data.frame() %>% mutate(p.value = round(p.value, 3))


AIC(mod_noQ, mod_Q)
anova(mod_noQ, mod_Q)

# Psuedo R2
rsq <- function (x, y) cor(x, y) ^ 2
rsq(mod_sf_unique_scaled$Ho, predict(mod_Q))
rsq(mod_sf_unique_scaled$Ho, predict(mod_noQ))

pdf(here("analysis", "wildfire", "plots", "wildfire_vs_het_1984_2022.pdf"), width = 4, height = 4)
ggplot(mod_sf_unique_scaled, aes(x = log1p(fire_values$conus_1984_2022_FRQ), y = Ho)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  labs(x = "log(recent burn count + 1)", y = "Ho") +
  ggpubr::stat_cor() +
  theme_classic()

ggplot(mod_sf_unique_scaled, aes(x = predict(mod_Q), y = Ho)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  labs(x = "Predicted Ho (Q model)", y = "Observed Ho") +
  ggpubr::stat_cor() +
  theme_classic()

ggplot(mod_sf_unique_scaled, aes(x = predict(mod_noQ), y = Ho)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  labs(x = "Predicted Ho (no Q model)", y = "Observed Ho") +
  ggpubr::stat_cor() +
  theme_classic()

ggplot(mod_sf_unique_scaled, aes(x = predict(mod_noQ), y = Ho)) +
  geom_point(size = 1, aes(col = Q)) +
  geom_smooth(method = "lm") +
  labs(x = "Predicted Ho (no Q model)", y = "Observed Ho") +
  ggpubr::stat_cor() +
  scale_color_viridis_c(option = "plasma") + 
  theme_classic()

outlier <- which(mod_sf_unique_scaled$Ho > 0.002 & predict(mod_noQ) < 0.001)
ggplot(mod_sf_unique_scaled[-outlier,], aes(x = predict(mod_noQ)[-outlier], y = Ho)) +
  geom_point(size = 1) +
  geom_smooth(method = "lm") +
  labs(x = "Predicted Ho (no Q model)", y = "Observed Ho") +
  ggpubr::stat_cor() +
  theme_classic()
dev.off()


mod_nest <- lm(residuals(mod_noQ) ~ Q, data = mod_sf_unique_scaled)
summary(mod_nest)
pdf(here("analysis", "wildfire", "plots", "wildfire_vs_het_1984_2022_resid.pdf"), width = 5, height = 4)
ggplot(mod_sf_unique_scaled) +
  geom_sf(data = ca) +
  geom_sf(aes(col = resid(mod_noQ)), size = 1) +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  theme_void()

ggplot(mod_sf_unique_scaled) +
  geom_sf(data = ca) +
  geom_sf(aes(col = resid(mod_nest)), size = 1) +
  scale_color_viridis_c(option = "plasma", direction = -1) +
  theme_void()
dev.off()
