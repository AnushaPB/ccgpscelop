library(tidyverse)
library(cowplot)
library(here)
library(sf)
outpath <- here("analysis", "check_nonsyn", "outputs")
plotpath <- here("analysis", "check_nonsyn", "plots")
source(here("general_functions.R"))

all_het <- read_table(here(outpath, "58-Sceloporus.het"), col_names = TRUE) %>% mutate(Ho =  1 - (`O(HOM)` / `N(NM)`))
nonsyn_het <- read_table(here(outpath, "all_nonsynonymous.het"), col_names = TRUE) %>% mutate(Ho_nonsyn = 1 - (`O(HOM)` / `N(NM)`))
syn_het <- read_table(here(outpath, "all_synonymous.het"), col_names = TRUE) %>% mutate(Ho_syn = 1 - (`O(HOM)` / `N(NM)`))

model_df <- read_csv(here("analysis", "genetic_diversity", "outputs", "model_df.csv"))
pop_df <- get_pops()

het <- 
  all_het %>%
  select(IID, Ho) %>%
  left_join(select(nonsyn_het, IID, Ho_nonsyn)) %>%
  left_join(select(syn_het, IID, Ho_syn)) %>%
  rename(SampleID = IID) %>%
  left_join(get_coords(sf = TRUE)) %>%
  left_join(select(model_df, -Ho)) %>%
  left_join(pop_df) %>%
  mutate(
    Ho_nonsyn_resid = resid(lm(Ho_nonsyn ~ Ho)),
    Ho_syn_resid = resid(lm(Ho_syn ~ Ho)),
    Ho_nonsyn_resid_syn = resid(lm(Ho_nonsyn ~ Ho_syn))
  ) %>%
  st_as_sf()

simple_plot <- function(x, y, data) {
  ggplot(data, aes_string(x = x, y = y)) +
    geom_point(size = 0.8) +
    theme_classic() +
    geom_smooth(method = "lm") +
    ggpubr::stat_cor() +
    labs(x = make_pretty_names(x), y = make_pretty_names(y)) +
    theme(
      axis.title = element_text(size = 14)
    )
}
pdf(here(plotpath, "het_check.pdf"), width=4, height=4)
simple_plot("Ho", "Ho_nonsyn", het)
simple_plot("tmean_dif", "Ho_nonsyn", het)
simple_plot("tmean_dif", "Ho_syn", het) 
simple_plot("tmean_dif", "Ho_nonsyn_resid", het) + ylab("Residuals(non-synonymous ~ genome-wide Ho)")
simple_plot("tmean_dif", "Ho_syn_resid", het)  + ylab("Residuals(synonymous ~ genome-wide Ho)")
simple_plot("tmean_dif", "Ho_nonsyn_resid_syn", het) + ylab("Residuals(non-synonymous ~ synonymous Ho)")
dev.off()


het_long <- het %>%
  select(SampleID, tmean_dif, Ho_nonsyn_resid, Ho_syn_resid) %>%
  pivot_longer(cols = c(Ho_nonsyn_resid, Ho_syn_resid), names_to = "type", values_to = "residuals") %>%
  mutate(
    type = case_when(
      type == "Ho_nonsyn_resid" ~ "Non-synonymous",
      type == "Ho_syn_resid" ~ "Synonymous"
    )
  )


color_values <- c("Non-synonymous" = "red2", "Synonymous" = "blue2", "All" = "#5f5f5f")

pdf(here(plotpath, "het_check_synonymous.pdf"), width=6, height=4.5)
ggplot(het_long, aes(x = tmean_dif, y = residuals, color = type)) +
  theme_classic() +
  geom_point(size = 0.8) +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor() +
  labs(x = "Recent temperature change", y = "Residual Ho", color = "Variant set") +
  scale_color_manual(values = color_values) +
  theme(
    axis.title = element_text(size = 14)
  )
dev.off()

pdf(here(plotpath, "het_check_clean.pdf"), width=3*5, height=5)
plot_grid(
  simple_plot("tmean_dif", "Ho_nonsyn", het),
  simple_plot("tmean_dif", "Ho_syn", het),
  simple_plot("tmean_dif", "Ho_nonsyn_resid_syn", het) + ylab("Residuals(non-synonymous ~ synonymous Ho)"),
  nrow = 1
)
dev.off()

southern_group <- c(6, 8)
eastern_group <- c(2, 1, 5)
western_group <- c(3, 4, 7)

het <- het %>%
  mutate(
    group = case_when(
      cluster %in% southern_group ~ "Southern",
      cluster %in% eastern_group ~ "Eastern",
      cluster %in% western_group ~ "Western",
      TRUE ~ "Other"
    )
  )


pdf(here(plotpath, "het_check_by_pop.pdf"), width=6, height=6)

ggplot(het, aes(x = tmean_dif, y = Ho_nonsyn_resid)) +
  geom_point(size = 1.5, aes(color = cluster)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = make_pretty_names("tmean_dif")) +
  ggpubr::stat_cor()

ggplot(het) + 
  geom_sf(aes(col = cluster), size = 3) 
 
ggplot(het, aes(x = tmean_dif, y = Ho_nonsyn_resid, color = cluster)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = make_pretty_names("tmean_dif")) +
  ggpubr::stat_cor(label.y = 0.015) +
  facet_wrap(~cluster) 
  

ggplot(het, aes(x = tmean_dif, y = Ho_nonsyn_resid)) +
  geom_point(size = 1.5, aes(color = group)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = make_pretty_names("tmean_dif")) +
  ggpubr::stat_cor()

ggplot(het) + 
  geom_sf(aes(col = group), size = 3) 
 
ggplot(het, aes(x = tmean_dif, y = Ho_nonsyn_resid, color = group)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = make_pretty_names("tmean_dif")) +
  ggpubr::stat_cor() + 
  facet_wrap(~group)

ggplot(filter(het, group != "Southern"), aes(x = tmean_dif, y = Ho_nonsyn_resid)) +
  geom_point(size = 1.5, aes(color = group)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = make_pretty_names("tmean_dif")) +
  ggpubr::stat_cor()

dev.off()


ca <- get_ca()

cluster_colors <- viridis::turbo(8)
names(cluster_colors) <- 1:8

plt1 <- ggplot(het) + 
  geom_sf(data = ca) +
  geom_sf(aes(col = cluster), size = 3) + 
  theme_void() +
  scale_color_manual(values = cluster_colors) 

plt2 <- ggplot(het, aes(x = tmean_dif, y = Ho_nonsyn_resid)) +
  geom_point(size = 1.5, aes(color = cluster)) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = make_pretty_names("tmean_dif"), y = "Residual(non-synonymous ~ genome-wide Ho)") +
  ggpubr::stat_cor() +
  theme_classic() +
  scale_color_manual(values = cluster_colors)

plt3 <- ggplot(filter(het, group != "Southern"), aes(x = tmean_dif, y = Ho_nonsyn_resid)) +
  geom_point(size = 1.5, aes(color = cluster)) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = make_pretty_names("tmean_dif"), y = "Residual(non-synonymous ~ genome-wide Ho)", title = "Southern clusters excluded") +
  ggpubr::stat_cor() + 
  theme_classic() +
  scale_color_manual(values = cluster_colors)

plt4 <- ggplot(filter(het, group == "Southern"), aes(x = tmean_dif, y = Ho_nonsyn_resid)) +
  geom_point(size = 1.5, aes(color = cluster)) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = make_pretty_names("tmean_dif"), y = "Residual(non-synonymous ~ genome-wide Ho)", title = "Southern clusters only") +
  ggpubr::stat_cor() + 
  theme_classic() +
  scale_color_manual(values = cluster_colors)


png(here(plotpath, "het_skeptical.png"), width=12*300, height=11*300, res = 300)
plot_grid(plt1, plt2, plt3, plt4, nrow = 2)
dev.off()
