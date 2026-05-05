library(tidyverse)
library(sf)
library(here)
source("general_functions.R")
coords <- read_csv(here("data", "WGS_METADATA_DB_04-29-2025.csv"))
old_coords <- get_coords(sf = TRUE)
specimens <- read_csv(here("specimens", "tracker.csv")) %>% st_as_sf(coords = c("x", "y"), crs = 4326) %>% drop_na(updated_note)
scelop_coords <- 
  coords %>% 
  filter(`ccgp-project-id` == "58-Sceloporus") %>%
  dplyr::select("long", "lat", `*organism`, `*sample_name`, `minicore_seq_id`) %>% 
  st_as_sf(
  coords = c("long", "lat"),
  crs = 4326
) %>%
rename(
  MinicoreID = `minicore_seq_id`,
  SampleID = `*sample_name`,
  Organism = `*organism`
) 

new_coords <- 
  scelop_coords %>%
  filter(!SampleID %in% old_coords$SampleID & !MinicoreID %in% old_coords$SampleID) 

new_specimens <-
  new_coords %>%
  filter(grepl("MVZ|CAS", SampleID))

new_specimens_not_in_old <-
  new_specimens %>%
  filter(!SampleID %in% specimens$SampleID & !MinicoreID %in% specimens$SampleID)
new_specimens_not_in_old

new_coords <- scelop_coords %>% filter(!SampleID %in% old_coords$SampleID & !MinicoreID %in% old_coords$SampleID)
nrow(new_coords)

ca <- get_ca()


plt1 <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = old_coords, aes(color = "old sequenced"), cex = 1) +
  scale_color_manual(values = c("old sequenced" = "black", "new sequenced" = "red", "old specimens" = "black", "new specimens" = "red")) +
  theme_void() +
  theme(legend.position = "top", legend.title = element_blank()) 

plt2 <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = new_coords, aes(color = "new sequenced"), cex = 1, pch = 1) +
  scale_color_manual(values = c("old sequenced" = "black", "new sequenced" = "red", "old specimens" = "black", "new specimens" = "red")) +
  theme_void() +
  theme(legend.position = "top", legend.title = element_blank()) 


plt3 <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = new_specimens, aes(color = "new specimens"), cex = 2, pch  = 4) +
  scale_color_manual(values = c("old sequenced" = "black", "new sequenced" = "red", "old specimens" = "black", "new specimens" = "red")) +
  theme_void() +
  theme(legend.position = "top", legend.title = element_blank()) 



plt4 <-
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = specimens, aes(color = "old specimens"), cex = 2, pch = 4) +
  scale_color_manual(values = c("old sequenced" = "black", "new sequenced" = "red", "old specimens" = "black", "new specimens" = "red")) +
  theme_void() +
  theme(legend.position = "top", legend.title = element_blank()) 


plt5 <- 
  ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = old_coords, aes(color = "old sequenced"), cex = 2) +
  geom_sf(data = new_coords, aes(color = "new sequenced"), cex = 2, pch = 1) +
  geom_sf(data = new_specimens, aes(color = "new specimens"), cex = 4, pch  = 4) +
  geom_sf(data = specimens, aes(color = "old specimens"), cex = 4, pch = 4) +
  scale_color_manual(values = c("old sequenced" = "black", "new sequenced" = "red", "old specimens" = "black", "new specimens" = "red")) +
  theme_void()


library(cowplot)
group1 <- plot_grid(plt1, plt4, plt2, plt3, ncol = 2)
group2 <- plot_grid(plt5, ncol = 1)
plot_grid(group1, group2, nrow = 1)


# Coordinate to remove
bad_coord <- scelop_coords %>% filter(grepl("1382", MinicoreID)) 

ggplot() +
  geom_sf(data = ca) +
  geom_sf(data = old_coords, cex = 1) +
  geom_sf(data = bad_coord, aes(color = "bad coord"), cex = 4, pch = 4) +
  scale_color_manual(values = c("bad coord" = "red")) +
  theme_void() +
  theme(legend.position = "top", legend.title = element_blank())

bad_coord %>% select(SampleID, MinicoreID) %>% st_drop_geometry()


# Sex check
sex_ids <- c(
  "Scelocci_IW1426",
  "Scelocci_IW2916",
  "Scelocci_IW3203",
  "Scelocci_Array_33",
  "Scelocci_NS19-002-S18",
  "Scelocci_MLY75",
  "Scelocci_JOS99",
  "Scelocci_MNW15-040-S84",
  "Scelocci_MNW16-025-S28",
  "Scelocci_IW2782",
  "Scelocci_CAS224151",
  "Scelocci_CAS251463",
  "Scelocci_IW2893",
  "Scelocci_LACM188840",
  "Scelocci_IW2779",
  "Scelocci_CAS227430",
  "Scelocci_MVZ_257313",
  "Scelocci_CAS252891",
  "Scelocci_IW3281",
  "Scelocci_IW3292",
  "Scelocci_IW3210",
  "Scelocci_CAS220918",
  "Scelocci_CAS227001",
  "Scelocci_HBS34959",
  "Scelocci_CAS212620",
  "Scelocci_LACM188004",
  "Scelocci_CAS225254",
  "Scelocci_IW3278",
  "Scelocci_CAS224099",
  "Scelocci_CAS252970",
  "Scelocci_WAP1301",
  "Scelocci_IW3213",
  "Scelocci_IW3284",
  "Scelocci_IW3291",
  "Scelocci_IW1501",
  "Scelocci_CAS212756",
  "Scelocci_MVZ_243334",
  "Scelocci_MVZ_232848",
  "Scelocci_IW3247",
  "Scelocci_IW2789",
  "Scelocci_CAS241780",
  "Scelocci_IW1537",
  "Scelocci_CAS262691",
  "Scelocci_HBS135925_or_135926",
  "Scelocci_HBS_136912",
  "Scelocci_MVZ_233482",
  "IW3488",
  "IW3491",
  "IW3223",
  "IW2907",
  "IW2909",
  "IW3027"
)

scelop_coords <- 
  coords %>% 
  filter(`ccgp-project-id` == "58-Sceloporus") %>%
  st_as_sf(
  coords = c("long", "lat"),
  crs = 4326
) %>%
select(
  MinicoreID = `minicore_seq_id`,
  SampleID = `*sample_name`,
  Organism = `*organism`,
  sex, `*sex`
) 

nona_sex <- scelop_coords %>% filter(sex != "NaN" | `*sex` != "NaN") %>% st_drop_geometry()

nona_sex %>% filter(SampleID %in% sex_ids) 
nona_sex %>% filter(!SampleID %in% sex_ids) 

write_csv(scelop_coords, "sceloporus_metadata.csv")
