
buffers <- c(0, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
names(buffers) <- paste0("buffer_", buffers)


buffer_dfs <-
  map(buffers, ~{
    buffer_dist <- .x

    fire_list <- list(pfs_replac, pfs_mixed, pfs_surfac, fri)
    names(fire_list) <- c("pfs_replac", "pfs_mixed", "pfs_surfac", "fri")

    # Make buffered points (use the mapped distance!)
    coords_buffer <- coords_proj %>% st_buffer(dist = buffer_dist)
    if (buffer_dist == 0) coords_buffer <- coords_proj

    # Get vegetation data
    groupveg_vals <- terra::extract(bps, coords_buffer, na.rm = FALSE)
    length(unique(groupveg_vals$ID)) == nrow(coords_proj)

    if (buffer_dist == 0) {
      open_water_ids <- which(groupveg_vals$GROUPVEG == "Open Water")
      open_water_vals <- 
        extract(bps, st_buffer(coords_proj[open_water_ids,], dist = 50)) %>%
        group_by(ID) %>%
        count(GROUPVEG) %>%
        filter(GROUPVEG != "Open Water") %>%
        slice_max(order_by = n, n = 1) %>%
        ungroup() %>%
        select(ID, GROUPVEG)

      groupveg_vals$GROUPVEG[open_water_ids] <- open_water_vals$GROUPVEG
    }

    groupveg_count <-
      groupveg_vals %>%
      filter(GROUPVEG != "Open Water") %>%
      count(ID, GROUPVEG) %>%
      mutate(prop = n / sum(n))

    groupveg_mode_with_ties <-
      groupveg_count %>%
      group_by(ID) %>%
      slice_max(order_by = prop, n = 1, with_ties = TRUE) %>%
      ungroup() %>%
      select(ID, GROUPVEG)

    tied_ids <-
      groupveg_mode_with_ties %>%
      group_by(ID) %>%
      filter(n() > 1) %>%
      pull(ID) %>%
      unique()

    # default if no ties
    groupveg_mode <- groupveg_mode_with_ties

    if (length(tied_ids) > 0) {
      tied_coords <- coords_proj[tied_ids, ]
      tied_groupveg_vals <- terra::extract(bps, tied_coords, na.rm = FALSE)
      tied_groupveg_vals$ID <- tied_ids

      groupveg_mode <-
        groupveg_mode_with_ties %>%
        filter(!(ID %in% tied_ids)) %>%
        bind_rows(tied_groupveg_vals) %>%
        arrange(ID) %>%
        select(ID, GROUPVEG)
    }

    stopifnot(nrow(groupveg_mode) == nrow(coords_proj))

    groupveg_mean <-
      groupveg_count %>%
      select(-n) %>%
      pivot_wider(
        names_from = GROUPVEG,
        values_from = prop,
        values_fill = 0
      ) 
    names(groupveg_mean)[-1] <- paste0("veg_", names(groupveg_mean)[-1])

    veg_wide <-
      groupveg_mean %>%
      ungroup() %>%
      select(-ID)

    pc <- prcomp(veg_wide, scale. = TRUE)

    allowed_groupveg <- unique(terra::extract(bps, coords_proj, na.rm = FALSE)$GROUPVEG)
    allowed_groupveg <- allowed_groupveg[allowed_groupveg != "Open Water"]

    fire_buffer_vals <-
      map(
        fire_list,
        ~ terra::extract(.x, coords_buffer, na.rm = FALSE),
        .progress = TRUE
      )

    fire_vals <-
      map(fire_buffer_vals, ~{
        fire_col <- setdiff(names(.x), "ID")
        stopifnot(length(fire_col) == 1)

        stopifnot(nrow(.x) == nrow(groupveg_vals))
        stopifnot(all(!is.na(groupveg_vals$GROUPVEG)))

        output <-
          .x %>%
          mutate(groupveg = groupveg_vals$GROUPVEG) %>%
          mutate(fire_num = as.numeric(as.character(.data[[fire_col]]))) %>%
          mutate(
            fire_clean = case_when(
              fire_num == -9999 & (groupveg %in% allowed_groupveg)  ~ 0,
              fire_num == -9999 & !(groupveg %in% allowed_groupveg) ~ NA_real_,
              TRUE ~ fire_num
            )
          ) %>%
          group_by(ID) %>%
          summarise(
            fire_mean = median(fire_clean, na.rm = TRUE),
            .groups = "drop"
          )

        names(output) <- c("ID", fire_col)
        output
      }, .progress = TRUE)

    vdep_vals <-
      terra::extract(vdep, coords_buffer) %>%
      rename(vdep = LABEL) %>%
      mutate(
        vdep = as.numeric(as.character(vdep)),
        vdep = ifelse(is.na(vdep), 101, vdep)
      ) %>%
      group_by(ID) %>%
      summarise(vdep = median(vdep, na.rm = TRUE), .groups = "drop")

    frq_vals <-
      terra::extract(frq, coords_buffer) %>%
      group_by(ID) %>%
      summarise(frq = max(conus_1984_2022_FRQ, na.rm = TRUE), .groups = "drop")

    fire_mod_df <-
      reduce(fire_vals, left_join, by = "ID") %>%
      left_join(groupveg_mode, by = "ID") %>%
      rename(groupveg = GROUPVEG) %>%
      left_join(groupveg_mean, by = "ID") %>%
      left_join(vdep_vals, by = "ID") %>%
      left_join(frq_vals, by = "ID") %>%
      select(-ID) %>%
      mutate(SampleID = coords_proj$SampleID) %>%
      rename(
        fri = FRI_ALLFIR,
        pfs_replac = PRC_REPLAC,
        pfs_mixed = PRC_MIXED,
        pfs_surfac = PRC_SURFAC
      ) %>%
      mutate(
        fri_cat = case_when(
          fri == 0 ~ "No Fire Regime",
          fri <= 35 ~ "Frequent Fire (≤35 yrs)",
          fri <= 200 ~ "Infrequent Fire (35–200 yrs)",
          fri > 200 ~ "Very Infrequent Fire (>200 yrs)"
        ),
        pfs_replac_cat = case_when(
          pfs_replac == 0 ~ "No Fire Regime",
          pfs_replac < 33 ~ "Low Replacement (<33%)",
          pfs_replac < 67 ~ "Moderate Replacement (33–66%)",
          pfs_replac >= 67 ~ "High Replacement (>66%)"
        ),
        pfs_sum  = pfs_replac + pfs_mixed + pfs_surfac,
        pfs_low  = (pfs_mixed + pfs_surfac) / pfs_sum,
        pfs_high = pfs_replac / pfs_sum
      )

    cor(
      log(fire_mod_df[, c("pfs_replac","pfs_surfac","pfs_low","pfs_high","fri","frq")] + 1),
      use = "pairwise.complete.obs"
    )

    fire_mod_df <-
      fire_mod_df %>%
      mutate(
        fire_severity = log1p(pfs_replac),
        fire_frequency = log1p(fri),
        fire_severity_cat = factor(
          pfs_replac_cat,
          levels = c("No Fire Regime","Low Replacement (<33%)","Moderate Replacement (33–66%)","High Replacement (>66%)")
        ),
        fire_frequency_cat = factor(
          fri_cat,
          levels = c("No Fire Regime","Frequent Fire (≤35 yrs)","Infrequent Fire (35–200 yrs)","Very Infrequent Fire (>200 yrs)")
        ),
        fire_recent = log1p(frq),
        groupveg = as.character(groupveg)
      )

    cor(fire_mod_df$fire_severity, fire_mod_df$fire_frequency, use = "pairwise.complete.obs")

    fire_pc <- prcomp(
      fire_mod_df[, c("fire_frequency", "fire_severity")],
      center = TRUE, scale. = TRUE
    )$x[, 1]

    fire_mod_df$fire_pc <- fire_pc * -1

    fire_mod_df
  }, .progress = TRUE)


buffer_mod_dfs <- map(buffer_dfs, ~{
  mod_df_sub <- mod_df %>% select(-any_of(names(.x)[-1]), SampleID)
  left_join(.x, mod_df_sub, by = c("SampleID")) %>% 
    group_by(x, y) %>% 
    summarize(
      across(
        c(Ho, bio1, csi_past, tmean_dif, NDVI, glacier, gHM, vdep,
          fire_severity, fire_frequency, fire_recent, fire_pc),
        mean,
        na.rm = TRUE
      ),
      groupveg = first(groupveg),
      evt      = first(evt),
      .groups = "drop"
    ) %>%
    st_as_sf(coords = c("x", "y"), crs = 4326) %>%
    st_transform(crs = 3310) %>%
    bind_cols(st_coordinates(.)) %>%
    ungroup() %>%
    # Scale continuous variables
    mutate(across(
      c(Ho, bio1, csi_past, tmean_dif, NDVI, glacier, gHM, vdep,
        fire_severity, fire_frequency, fire_recent, fire_pc),
      ~ as.numeric(scale(.))
    ))
}, .progress = TRUE)

gls_mods_veg <- map(buffer_mod_dfs, ~{
  f <- formula("Ho ~ bio1 + csi_past + tmean_dif  + glacier + gHM + fire_severity + fire_recent + vdep + groupveg + evt")
  model <- gls(f, correlation = corExp(form = ~ X + Y, nugget = FALSE), data = .x, method = "REML")
  model_result <- model %>% broom::tidy(conf.int = TRUE, conf.level = 0.95) %>% mutate(model = "Vegetation group")
}, .progress = TRUE)

gls_mods_noveg <- map(buffer_mod_dfs, ~{
  f <- formula("Ho ~ bio1 + csi_past + tmean_dif  + glacier + gHM + fire_severity + fire_recent + vdep + evt")
  model <- gls(f, correlation = corExp(form = ~ X + Y, nugget = FALSE), data = .x, method = "REML")
  model_result <- model %>% broom::tidy(conf.int = TRUE, conf.level = 0.95) %>% mutate(model = "No vegetation group")
}, .progress = TRUE)

gls_results <- 
  bind_rows(gls_mods_veg, .id = "buffer") %>%
  bind_rows(bind_rows(gls_mods_noveg, .id = "buffer")) %>%
  mutate(buffer = as.numeric(sub("buffer_", "", buffer))) %>%
  mutate(Significant = ifelse(p.value < 0.05, "Yes", "No"))

pdf(here("analysis/genetic_diversity/plots/fire_buffer_model_results.pdf"), width = 5, height = 6)
ggplot(gls_results %>% filter(term == "fire_severity"), aes(x = buffer, y = estimate))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_line() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 10) +
  geom_point(aes(fill = Significant), pch = 21, cex = 2) +
  theme_classic() +
  facet_wrap(~model, ncol = 1) +
  scale_fill_manual(values = c("Yes" = "black", "No" = "white")) +
  labs(y = "Coefficient of fire severity", x = "Buffer radius (m)")

ggplot(gls_results %>% filter(term == "fire_recent"), aes(x = buffer, y = estimate))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_line() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 10) +
  geom_point(aes(fill = Significant), pch = 21, cex = 2) +
  theme_classic() +
  facet_wrap(~model, ncol = 1) +
  scale_fill_manual(values = c("Yes" = "black", "No" = "white")) +
  labs(y = "Coefficient of fire severity", x = "Buffer radius (m)")
dev.off()


#
gls_model <- gls(Ho ~ bio1 + csi_past + tmean_dif + NDVI + glacier + gHM + fire_severity + fire_recent + vdep + groupveg + evt,
    correlation = corExp(form = ~ X + Y, nugget = TRUE),
    data = buffer_mod_dfs$buffer_100,
    method = "REML"
)

ols_model <- lm(Ho ~ bio1 + csi_past + tmean_dif + NDVI + glacier + gHM + fire_severity + fire_recent + vdep + groupveg + evt,
    data = buffer_mod_dfs$buffer_100
)

anova(AIC(ols_model, gls_model))

df <- buffer_mod_dfs$buffer_100
df$gls_resid <- resid(gls_model)
df$ols_resid <- resid(ols_model)

pdf(here("analysis/genetic_diversity/plots/fire_buffer_model_residuals.pdf"), width = 8, height = 4)
ggplot(df, aes(col = gls_resid)) +
  geom_sf() +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "GLS Model Residuals")
ggplot(df, aes(col = ols_resid)) +
  geom_sf() +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  labs(title = "OLS Model Residuals")
dev.off()


# 

evt2 <- evt
activeCat(evt2) <- "EVT_NAME"

buffer_sizes <- c(50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
names(buffer_sizes) <- paste0("buffer_", buffer_sizes)

shannon_df <- 
  map(buffer_sizes, ~{
    buffer_dist <- .x

    # Make buffered points (use the mapped distance!)
    coords_buffer <- coords_proj %>% st_buffer(dist = buffer_dist)

    # Get counts of each EVT class per ID
    evt_shannon <- 
      terra::extract(evt2, coords_buffer, na.rm = FALSE, ID = TRUE) %>%
      filter(!is.na(EVT_NAME)) %>%
      count(ID, EVT_NAME) %>%
      group_by(ID) %>%
      mutate(p = n / sum(n)) %>%
      summarise(
        evt_shannon = -sum(p * log(p)),
        evt_evenness = evt_shannon / log(n_distinct(EVT_NAME)),
        evt_simple = n_distinct(EVT_NAME),
        .groups = "drop"
      ) %>%
      mutate(evt_evenness = ifelse(is.nan(evt_evenness), 0, evt_evenness)) %>%
      mutate(SampleID = coords_proj$SampleID) %>%
      left_join(read_csv(here("analysis/genetic_diversity/outputs/model_df.csv")), by = "SampleID") %>%
      # Scale Ho, csi_past, tmean_dif, bio1
      mutate(
        Ho = as.numeric(scale(Ho)),
        csi_past = as.numeric(scale(csi_past)),
        tmean_dif = as.numeric(scale(tmean_dif)),
        bio1 = as.numeric(scale(bio1))
      )

      mod1 <- lm(Ho ~ evt_shannon + csi_past + tmean_dif + bio1 + fire_frequency_cat, data = evt_shannon)
      mod2 <- lm(Ho ~ evt_evenness + csi_past + tmean_dif + bio1 + fire_frequency_cat, data = evt_shannon)
      mod3 <- lm(Ho ~ evt_simple + csi_past + tmean_dif + bio1 + fire_frequency_cat, data = evt_shannon)
      
      bind_rows(
        broom::tidy(mod1),
        broom::tidy(mod2),
        broom::tidy(mod3)
      ) %>% filter(term != "(Intercept)") 
  }, .progress = TRUE) 

shannon_results <- 
  bind_rows(shannon_df, .id = "buffer") %>%
  mutate(buffer = as.numeric(sub("buffer_", "", buffer)))

# Check for any significant relationships
shannon_results %>%
  filter(term %in% c("evt_shannon", "evt_evenness", "evt_simple")) %>%
  filter(p.value < 0.05)

shannon_results %>%
  filter(buffer == 1000) %>%
  filter(p.value < 0.05)


# Make buffered points (use the mapped distance!)
coords_buffer <- coords_proj %>% st_buffer(dist = 1000)

# Get counts of each EVT class per ID
evt_shannon <- 
  terra::extract(evt2, coords_buffer, na.rm = FALSE, ID = TRUE) %>%
  filter(!is.na(EVT_NAME)) %>%
  count(ID, EVT_NAME) %>%
  group_by(ID) %>%
  mutate(p = n / sum(n)) %>%
  summarise(
    evt_shannon = -sum(p * log(p)),
    evt_evenness = evt_shannon / log(n_distinct(EVT_NAME)),
    evt_simple = n_distinct(EVT_NAME),
    .groups = "drop"
  ) %>%
  mutate(evt_evenness = ifelse(is.nan(evt_evenness), 0, evt_evenness)) %>%
  mutate(SampleID = coords_proj$SampleID) %>%
  left_join(read_csv(here("analysis/genetic_diversity/outputs/model_df.csv")), by = "SampleID") %>%
  # Scale Ho, csi_past, tmean_dif, bio1
  mutate(
    Ho = as.numeric(scale(Ho)),
    csi_past = as.numeric(scale(csi_past)),
    tmean_dif = as.numeric(scale(tmean_dif)),
    bio1 = as.numeric(scale(bio1))
  )

lm(Ho ~ evt_shannon + csi_past + tmean_dif + bio1 + fire_frequency_cat, data = evt_shannon) %>%
  summary()

# Get partial evt_shannon effect
evt_shannon_resid <-
  lm(evt_shannon ~ csi_past + tmean_dif + bio1 + fire_frequency_cat, data = evt_shannon) %>%
  resid()

evt_shannon_Ho_resid <- 
  lm(Ho ~ csi_past + tmean_dif + bio1 + fire_frequency_cat, data = evt_shannon) %>%
  resid()

evt_shannon <-
  evt_shannon %>%
  mutate(evt_shannon_resid = evt_shannon_resid, 
         Ho_resid = evt_shannon_Ho_resid)

pdf(here("analysis/genetic_diversity/plots/evt_shannon_buffer_effects.pdf"), width = 5, height = 4)
ggplot(evt_shannon, aes(x = evt_shannon, y = Ho)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  labs(title = "Ho vs. EVT Shannon Diversity (1000m buffer)", x = "EVT Shannon Diversity", y = "Ho") +
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top")

ggplot(evt_shannon, aes(x = evt_shannon_resid, y = Ho_resid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  labs(title = "Partial Ho vs. Partial EVT Shannon Diversity (1000m buffer)", x = "Partial EVT Shannon Diversity", y = "Partial Ho") +
  ggpubr::stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top")
dev.off()
