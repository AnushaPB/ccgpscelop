manual_preds_Ho <- function(model, data,
                            tmean_vals = NULL,
                            bio1_vals = NULL,
                            evt_level = NULL,
                            fire_recent_level = NULL) {
  data <- as.data.frame(data)
  b <- coef(model)

  # defaults
  if (is.null(tmean_vals)) {
    tmean_vals <- seq(min(data$tmean_dif, na.rm = TRUE),
                      max(data$tmean_dif, na.rm = TRUE),
                      length.out = 200)
  }

  if (is.null(bio1_vals)) {
    bio1_vals <- as.numeric(quantile(data$bio1, c(0.25, 0.5, 0.75), na.rm = TRUE))
  }

  # choose reference / fixed values
  glacier0   <- mean(data$glacier, na.rm = TRUE)
  gHM0       <- mean(data$gHM, na.rm = TRUE)
  csi0       <- mean(data$csi_past, na.rm = TRUE)
  fire_frq0  <- mean(data$fire_frq, na.rm = TRUE)

  if (is.null(evt_level)) {
    evt_level <- levels(data$evt)[1]
  }
  if (is.null(fire_recent_level)) {
    fire_recent_level <- levels(data$fire_recent)[1]
  }

  # helper to get coefficient if present, else 0
  getb <- function(name) if (name %in% names(b)) unname(b[name]) else 0

  # polynomial terms
  poly_obj <- poly(data$csi_past, 2, raw = TRUE)
  # for raw=TRUE, term 1 = x, term 2 = x^2
  csi1 <- csi0
  csi2 <- csi0^2

  out <- lapply(seq_along(bio1_vals), function(i) {
    bio <- bio1_vals[i]

    df <- data.frame(
      tmean_dif = tmean_vals,
      bio1 = bio,
      glacier = glacier0,
      gHM = gHM0,
      csi_past = csi0,
      evt = evt_level,
      fire_frq = fire_frq0,
      fire_recent = fire_recent_level
    )

    # start with intercept
    fit <- rep(getb("(Intercept)"), nrow(df))

    # numeric main effects
    fit <- fit + getb("glacier") * df$glacier
    fit <- fit + getb("gHM") * df$gHM
    fit <- fit + getb("tmean_dif") * df$tmean_dif
    fit <- fit + getb("bio1") * df$bio1
    fit <- fit + getb("fire_frq") * df$fire_frq

    # raw polynomial terms for csi_past
    fit <- fit + getb("poly(csi_past, 2, raw = TRUE)1") * csi1
    fit <- fit + getb("poly(csi_past, 2, raw = TRUE)2") * csi2

    # interaction: tmean_dif * bio1
    fit <- fit + getb("tmean_dif:bio1") * df$tmean_dif * df$bio1

    # factor: evt
    evt_coef_name <- paste0("evt", evt_level)
    fit <- fit + getb(evt_coef_name)

    # factor: fire_recent
    fr_coef_name <- paste0("fire_recent", fire_recent_level)
    fit <- fit + getb(fr_coef_name)

    # interaction: fire_frq * fire_recent
    # coefficient name usually "fire_frq:fire_recentLEVEL"
    fr_int_name1 <- paste0("fire_frq:fire_recent", fire_recent_level)
    fr_int_name2 <- paste0("fire_recent", fire_recent_level, ":fire_frq")

    fit <- fit + getb(fr_int_name1) * df$fire_frq
    fit <- fit + getb(fr_int_name2) * df$fire_frq

    df$.fitted <- fit
    df$.group_value <- bio
    df
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL

  qlabs <- c("Cold (Q1)", "Median", "Hot (Q3)")
  out$.group <- factor(out$.group_value,
                       levels = bio1_vals,
                       labels = qlabs[seq_along(bio1_vals)])
  out
}

preds <- manual_preds_Ho(
  model = final_model,
  data = df_copy,
  bio1_vals = as.numeric(quantile(df_copy$bio1, c(0.25, 0.5, 0.75), na.rm = TRUE))
)

f <- formula("Ho ~ glacier + gHM + poly(csi_past, 2, raw = TRUE) + tmean_dif * bio1 + evt + fire_frq * fire_recent")

# final_model <- gls(f, data = mod_sf, method = "REML")
# coef(final_model)["tmean_dif:bio1"] 

final_model <- lm(f, data = mod_sf)

b <- coef(final_model)
V <- vcov(final_model)

bio_seq <- seq(
  min(mod_sf$bio1, na.rm = TRUE),
  max(mod_sf$bio1, na.rm = TRUE),
  length.out = 200
)

slope_df <- data.frame(bio1 = bio_seq)

# slope of tmean_dif at each bio1:
# dY/d(tmean_dif) = beta_tmean_dif + beta_tmean_dif:bio1 * bio1
slope_df$slope <- b["tmean_dif"] + b["tmean_dif:bio1"] * slope_df$bio1

# SE for linear combination a' beta, where a = [1, bio1]
slope_df$se <- sapply(slope_df$bio1, function(z) {
  a <- c(1, z)
  names(a) <- c("tmean_dif", "tmean_dif:bio1")
  sqrt(t(a) %*% V[names(a), names(a)] %*% a)
})

# 95% CI
crit <- qt(0.975, df = df.residual(final_model))

slope_df <- slope_df %>%
  dplyr::mutate(
    conf.low = slope - crit * se,
    conf.high = slope + crit * se
  )

# Get value where slope switches from positive to negative
switch_point <- -b["tmean_dif"] / b["tmean_dif:bio1"]
switch_point

png(here("TEST.png"), width = 4, height = 4, units = "in", res = 300) 
ggplot(slope_df, aes(x = bio1, y = slope)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, col = "gray") +
  geom_vline(xintercept = switch_point, linetype = "dashed", color = "red") +
  annotate("text", x = switch_point, y = max(slope_df$slope), color = "red", vjust = 0, hjust = 1.1, label = paste0(round(switch_point, 2), " °C")) +
  theme_classic() +
  labs(
    x = make_pretty_names("bio1"),
    y = "Effect of climate change on Ho"
  )
dev.off()

# Split bio1 into min, halfway between min and siwthc point, switch point, halfway between switch point and max, and max
# compute breakpoints ONCE
bio_min  <- min(mod_sf$bio1, na.rm = TRUE)
bio_max  <- max(mod_sf$bio1, na.rm = TRUE)

mid_low  <- (bio_min + switch_point) / 2
mid_high <- (switch_point + bio_max) / 2

breaks <- c(bio_min, mid_low, switch_point, mid_high, bio_max)

fmt <- function(x) paste0(format(round(x, 2), nsmall = 2), " °C")
labels <- c(
  paste0(fmt(bio_min), " – ", fmt(mid_low)),
  paste0(fmt(mid_low), " – ", fmt(switch_point)),
  paste0(fmt(switch_point), " – ", fmt(mid_high)),
  paste0(fmt(mid_high), " – ", fmt(bio_max))
)

df_copy2 <- mod_sf %>%
  mutate(
    bio1_bin = cut(
      bio1,
      breaks = breaks,
      include.lowest = TRUE,
      labels = labels
    )
  )

slope_df2 <-
  slope_df %>%
  mutate(
    bio1_bin = cut(
      bio1,
      breaks = breaks,
      include.lowest = TRUE,
      labels = labels
    )
  )

# Get midpoint of each break (or you will have five intervals)
midbreaks <- (breaks[-1] + breaks[-length(breaks)]) / 2
preds <- ggeffects::ggpredict(
  final_model,
  terms = c(
    "tmean_dif [all]",
    paste0("bio1 [", paste(midbreaks, collapse = ", "), "]")
  ),
  data = mod_sf
)
 
png(here("TEST.png"), width = 4, height = 4, units = "in", res = 300)
plot(preds)
dev.off()

#c("#6bb2cc", "#042346", "#f04f32", "#f7a957", "#295580")

cols <- c("#042346","#6bb2cc", "#f7a957", "#e24923")


plt_slope <-
  ggplot(slope_df2, aes(x = bio1, y = slope)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = bio1_bin), alpha = 0.2) +
  geom_line(linewidth = 1.2, aes(col = bio1_bin)) +
  geom_hline(yintercept = 0, col = "gray") +
  geom_vline(xintercept = breaks, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = switch_point, linetype = "dashed", color = "red") +
  annotate("text", x = switch_point, y = max(slope_df$slope), color = "red", vjust = 0, hjust = 1.1, label = paste0(round(switch_point, 2), " °C")) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_classic() +
  labs(
    x = make_pretty_names("bio1"),
    y = "Effect of climate change"
  ) +
  theme(legend.position = "none")

plt_interact <-
  ggplot(preds, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  theme_classic() +
  labs(
    x = make_pretty_names("tmean_dif"),
    y = "Predicted heterozygosity",
  ) +
  theme(legend.position = "none") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)

plt_map <-
  ggplot(df_copy2) +
  geom_sf(data = ca) +
  geom_sf(aes(col = bio1_bin)) +
  scale_color_manual(values = cols) +
  theme_void() +
  labs(col = "Contemporary\ntemperature") +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.75))

plt_raw <-
  ggplot(df_copy2, aes(x = tmean_dif, y = Ho, col = bio1_bin)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, col = "black") +
  theme_minimal() +
  labs(x = make_pretty_names("tmean_dif"), y = "Observed Heterozygosity") +
  theme(legend.title = element_blank()) + 
  facet_wrap(~ bio1_bin, nrow = 1) +
  scale_color_manual(values = cols) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "none", 
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
    axis.line = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 13)
  )   

png(here("TEST.png"), width = 10, height = 6, units = "in", res = 300)
group1 <- plot_grid(plt_interact, plt_map, ncol = 1, labels = c("B", "D"))
group2 <- plot_grid(plt_slope, plt_raw, ncol = 1, labels = c("A", "C"), align = "v")
plot_grid(group2, group1, nrow = 1, rel_widths = c(2, 1))
dev.off()

future <- rast(here("analysis", "gea", "outputs", "scelop_adaptive_env_layers","env_fut_2071-2100_GFDL-ESM4_ssp126_ssp585.tif"))[[2]]
df_copy2$future_bio1 <- extract(future, df_copy2, ID = FALSE, method = "bilinear")[,1]   
df_copy2$future_change <- df_copy2$future_bio1 - df_copy2$bio1

range(df_copy2$future_bio1, na.rm = TRUE)
unique(df_copy2$future_bin)

plt_future <-
  ggplot(df_copy2, aes(x = future_change, y = Ho, col = bio1_bin)) +
  geom_point(alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, col = "black") +
  theme_minimal() +
  labs(x = "Future temperature change", y = "Observed Heterozygosity") +
  theme(legend.title = element_blank()) + 
  facet_wrap(~ bio1_bin, nrow = 1) +
  scale_color_manual(values = cols) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "none", 
    panel.background = element_rect(fill = "white", color = "black", linewidth = 1),
    axis.line = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    plot.margin = margin(5, 5, 5, 13)
  )  
png(here("TEST.png"), width = 10, height = 5, units = "in", res = 300)
plt_future
dev.off()
### ARIDITY
ai <- rast(here("data/env/Global-AI_ET0_v3_annual/ai_v3_yr.tif"))
et <- rast(here("data/env/Global-AI_ET0_v3_annual/et0_v3_yr.tif"))
ai_vals <- terra::extract(ai, mod_sf_unique_scaled[, c("x", "y")])[[2]]
et_vals <- terra::extract(et, mod_sf_unique_scaled[, c("x", "y")])[[2]]

cwd <- rast(here("data/env/cwd_wy2016.nc"))
cwd_annual <- sum(cwd, na.rm = TRUE)

vpd <- rast(list.files(here("data/env/chelsa_vpd"), full.names = TRUE)) 
vpd_vals <- terra::extract(vpd, mod_sf_unique_scaled)[,-1]
stopifnot(ncol(vpd_vals) == 12)
vpd_annual <- rowMeans(vpd_vals, na.rm = TRUE)

cmi <- rast(list.files(here("data/env/chelsa_cmi"), full.names = TRUE))
cmi_vals <- terra::extract(cmi, mod_sf_unique_scaled)[,-1]
stopifnot(ncol(cmi_vals) == 12)
cmi_annual <- rowMeans(cmi_vals, na.rm = TRUE)

cwd_vals <- terra::extract(cwd_annual, mod_sf_unique_scaled[, c("x", "y")])[[2]]

df_copy <- mod_sf_unique_scaled %>% mutate(ai = ai_vals, cwd = cwd_vals, vpd = vpd_annual, cmi = cmi_annual)

cor.test(df_copy$ai, df_copy$bio1) # -0.47

cor.test(df_copy$cwd, df_copy$Ho) # 0.68
cor.test(df_copy$vpd, df_copy$bio1) # 0.67
cor.test(df_copy$cmi, df_copy$bio1) # -0.50

png(here("ai_plot.png"))
ggplot(df_copy) +
  geom_sf(aes(col = ai)) +
  scale_color_viridis_c(option = "plasma", end = 0.1, begin = 0.9, direction = -1) +
  theme_void() +
  labs(title = "Aridity Index")

dev.off()


f <- formula("Ho ~ glacier + gHM + poly(csi_past, 2, raw = TRUE) + tmean_dif * bio1 + evt + fire_frq * fire_recent") 

# Fit model
final_model <- gls(f, correlation = corExp(form = ~ x + y, nugget = FALSE), data = df_copy, method = "ML")

library(ggeffects)

final_model <- lm(Ho ~  glacier + gHM + poly(csi_past, 2, raw = TRUE) + tmean_dif * bio1 + evt + fire_frq * fire_recent, data = df_copy)

preds <- ggpredict(
  final_model,
  terms = c("tmean_dif [all]", "bio1 [quart]"),
  data = df_copy
)

png(here("TEST.png"), width = 4, height = 4, units = "in", res = 300)
plot(preds)
dev.off()

# Results
mod_result <- broom::tidy(final_model)
mod_result
mod_result %>% filter(term == "vpd")
# Significant variables
mod_result %>% filter(p.value < 0.05) %>% arrange((p.value))
sig_vars <- mod_result %>% filter(p.value < 0.05, !term %in% c("(Intercept)", "lambda")) %>% pull(term)
print(sig_vars)
