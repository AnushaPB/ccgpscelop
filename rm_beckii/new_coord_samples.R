library(tidyverse)
library(sf)
library(here)
source("general_functions.R")
coords <- read_csv(here("rm_beckii", "WGS_METADATA_DB_04-29-2025.csv"))
scelop_coords <- 
  coords %>% 
  filter(`ccgp-project-id` == "58-Sceloporus") %>%
  dplyr::select("long", "lat", `*organism`, `*sample_name`, `minicore_seq_id`) %>% 
  st_as_sf(
  coords = c("long", "lat"),
  crs = 4326
) %>%
rename(
  SampleID = `*sample_name`,
  Organism = `*organism`,
  MinicoreID = `minicore_seq_id`
)

# Plot beckii vs occidentalis
ggplot(scelop_coords) +
  geom_sf(aes(col = Organism)) +
  theme_void()
 
scelop_coords <-
  scelop_coords %>%
  mutate(
    corrected_species = case_when(
      Organism == "Sceloporus becki" ~ "beckii",
      # This was labelled as Sceloporus occidentalis in the database, but it is actually Sceloporus
      # beckii
      SampleID == "MW01-3_-14" ~ "beckii",
      TRUE ~ "occidentalis"
    )
  )

pdf(here("rm_beckii", "scelop_coords.pdf"), width = 20, height = 20)
ggplot(scelop_coords) +
  #geom_sf(data = ca) +
  geom_sf(aes(col = corrected_species)) +
  #geom_sf(data = old_coords, color = "black") +
  ggrepel::geom_text_repel(data = scelop_coords, aes(label = SampleID, geometry = geometry), size = 3, max.overlaps = 100, stat = "sf_coordinates") +
  theme_void()
dev.off()

beckii <- 
  scelop_coords %>% 
  filter(corrected_species == "beckii") %>%
  st_drop_geometry() %>% 
  pull(MinicoreID)

writeLines(beckii, here("rm_beckii", "beckii_sampleids.txt"))
