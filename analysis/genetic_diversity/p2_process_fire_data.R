
library(here)
library(tidyverse)
library(sf)
library(terra)
source(here("general_functions.R"))

outpath <- here("analysis", "genetic_diversity", "outputs")
plotpath <- here("analysis", "genetic_diversity", "plots")

# California boundary -----------------------------------------------------------------------------------
ca <- get_ca()
ca_proj <- ca %>% st_transform(3310)

# Coordinates -------------------------------------------------------------------------------------------
coords <- get_coords(sf = TRUE)
coords_proj <- coords %>% st_transform(crs = st_crs(ca_proj))

# Vegetation type (EVT) ---------------------------------------------------------------------------------
ca_proj_buffer <-  ca_proj %>% st_transform(5070) %>% st_buffer(dist = 10000) # buffer to ensure full coverage when cropping 

evt_og <- rast(here("data", "env", "LF2016_EVT_200_CONUS", "Tif", "LC16_EVT_200.tif")) 
evt <- crop(evt_og, ca_proj_buffer)
# EVT PHYS is EVT physiogomy 
# It is a good balance between detail and generalization that still captures relevant vegetation types for sceloporus habitat
activeCat(evt) <- "EVT_PHYS"
names(evt) <- "evt"

# Extract EVT values in a 100 m buffer around each point
coords_buffer <- coords_proj %>% st_buffer(dist = 100)

# Get counts of each EVT class per ID
evt_vals <- terra::extract(evt, coords_buffer, na.rm = FALSE, ID = TRUE)

evt_counts <-
  evt_vals  %>% 
  # Remove open water and roads
  filter(!evt %in% c("Open Water", "Developed-Roads")) %>%
  count(ID, evt, name = "n")     # counts per ID x class

# This keeps all top classes per ID (i.e., ties are retained)
evt_mod_df_with_ties <-
  evt_counts %>%
  group_by(ID) %>%
  slice_max(order_by = n, n = 1, with_ties = TRUE) %>%
  ungroup()

# IDs that have ties for the top class
ties <-
  evt_mod_df_with_ties %>%
  add_count(ID, name = "n_top") %>%  # how many "winners" for that ID
  filter(n_top > 1) %>%
  arrange(ID, desc(n), evt)

# For ties, take the class of the original coordinates
tie_ids <- unique(ties$ID)
tie_coords <- coords_proj[tie_ids, ]
tie_evt_vals <- 
  terra::extract(evt, tie_coords, na.rm = FALSE) %>%
  mutate(ID = tie_ids)

# Replace tied classes with the value from the original coordinates
evt_mod_df <-
  evt_mod_df_with_ties %>%
  filter(!(ID %in% tie_ids)) %>%  # remove tied IDs
  bind_rows(tie_evt_vals %>% dplyr::select(ID, evt)) %>%  # add back in tied IDs with original coord values
  arrange(ID) %>%
  mutate(SampleID = coords_proj$SampleID) %>%
  dplyr::select(-ID) 

# Get counts in each category
evt_mod_df %>% 
  count(evt) %>%
  arrange(desc(n))
#    evt                            n
#    <fct>                      <int>
#  1 Shrubland                     83
#  2 Developed                     41
#  3 Conifer                       38
#  4 Exotic Herbaceous             26
#  5 Hardwood                      16
#  6 Riparian                      14
#  7 Sparsely Vegetated             9
#  8 Developed-Medium Intensity     7
#  9 Agricultural                   7
# 10 Exotic Tree-Shrub              6
# 11 Developed-Low Intensity        5
# 12 Grassland                      3
# 13 Developed-Roads                2
# Based on these counts, I think it is safe to collapse all developed classes into one "Developed" category

# Replace any evt classes that are developed with just "Developed"
# Doing this collapsing after the counting to avoid losing information about the other vegetation types present in the buffer
evt_mod_df <-
  evt_mod_df %>%
  mutate(evt = as.character(evt)) %>%
  mutate(evt = ifelse(grepl("Developed", evt), "Developed", evt)) %>%
  dplyr::select(-n)

# Count categories
cat_count <- evt_mod_df %>% count(evt) %>% arrange(desc(n))
cat_count

# Order levels by count
evt_mod_df$evt <- factor(evt_mod_df$evt, levels = cat_count$evt)

# Check NAs
evt_mod_df %>%
  summarise(across(everything(), ~sum(is.na(.))))

# Historical fire regime --------------------------------------------------------------------------------------

ca_proj_buffer <-  ca_proj %>% st_transform(5070) %>% st_buffer(dist = 10000) # buffer to ensure full coverage when cropping 

# Historical fire regime data
# Fire Return Interval (FRI)
fri_og <- rast(here("data", "env", "LF2016_FRI_200_CONUS", "Tif", "LC16_FRI_200.tif")) 
fri <- crop(fri_og, ca_proj_buffer)
activeCat(fri) <- "FRI_ALLFIR"

# Percent Fire Severity (PFS) 
pfs_og <- rast(here("data", "env", "LF2016_PFS_200_CONUS", "Tif", "LC16_PFS_200.tif"))
pfs <- crop(pfs_og, ca_proj_buffer)

# Vegetation group (BPS)
bps_og <- rast(here("data", "env", "LF2016_BPS_200_CONUS", "Tif", "LC16_BPS_200.tif"))
bps <- crop(bps_og, ca_proj_buffer)
activeCat(bps) <- "GROUPVEG"

# Longest fire free interval from 1984-2022
lffi_og <- rast(here("data", "env", "conus_fire_history_metrics_1984_2020", "conus_1984_2020_LFFI.tif"))
lffi <- crop(lffi_og, ca_proj_buffer)
names(lffi) <- "lffi"

# Recent fire frequency (FRQ) from 1984-2022
frq_og <- rast(here("data", "env", "conus_fire_history_metrics_1984_2020", "conus_1984_2020_FRQ.tif"))
frq <- crop(frq_og, ca_proj_buffer)
names(frq) <- "frq"

# Only using percent replacement for fire severity (high severity fire) - other PFS categories are low/moderate and together those percentages (low + moderate + high) sum to 100, so high severity captures the relevant variation (high = 100 - low + moderate)
pfs_replac <- pfs
activeCat(pfs_replac) <- "PRC_REPLAC"

# Make buffered points
coords_buffer <- coords_proj %>% st_buffer(dist = 100)

# Get vegetation data
groupveg_vals <- terra::extract(bps, coords_buffer, na.rm = FALSE)

# Count each vegetation type within buffer
groupveg_count <- 
  groupveg_vals %>%
  filter(GROUPVEG != "Open Water") %>%  # remove open water points
  count(ID, GROUPVEG) %>%               # count points per vegetation type
  group_by(ID) 

# Get mode vegetation type per ID
groupveg_mode_with_ties <- 
  groupveg_count %>%
  group_by(ID) %>%
  slice_max(order_by = n, n = 1, with_ties = TRUE) %>%
  ungroup() %>%
  select(ID, GROUPVEG)

# Pull out tied values and replace with the coordinate value
tied_ids <- 
  groupveg_mode_with_ties %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  pull(ID) %>%
  unique()

if (length(tied_ids) > 0) {
  tied_coords <- coords_proj[tied_ids, ]
  tied_groupveg_vals <- terra::extract(bps, tied_coords, na.rm = FALSE)
  tied_groupveg_vals$ID <- tied_ids # ensure ID column is correct
  # Replace in groupveg_mode
  groupveg_mode <- 
    groupveg_mode_with_ties %>%
    filter(!(ID %in% tied_ids)) %>%
    bind_rows(tied_groupveg_vals) %>%
    arrange(ID)
}

# Confirm that number of rows matches number of coords
stopifnot(nrow(groupveg_mode) == nrow(coords_proj))

# Get unique vegetation values for each coordinate to get a vector of allowed GROUPVEG types
allowed_groupveg <- unique(terra::extract(bps, coords_proj, na.rm = FALSE)$GROUPVEG)
# Remove open water 
allowed_groupveg <- allowed_groupveg[allowed_groupveg != "Open Water"]
# Anything that is -9999 and NOT in allowed_groupveg will be coded as NA (since it is not sceloporus habitat), any other value of -9999 will get a 0 (no fire regime)
print(allowed_groupveg)

# Create fire regime groups based on FRI and PFS based on an adaptation of bins from LANDFIRE
# NEW FRG bins (Too fine scale): https://www.landfire.gov/sites/default/files/DataDictionary/LF2016/LF16_FRGADD.pdf
# Adapted based on OLD FRG bins: https://www.landfire.gov/sites/default/files/DataDictionary/LF2014/LF14_FRGADD.pdf

# Create list of historical fire variables
pfs_replac_vals <- terra::extract(pfs_replac, coords_buffer, na.rm = FALSE, ID = TRUE)
fri_vals <- terra::extract(fri, coords_buffer, na.rm = FALSE, ID = TRUE)

frg_mode_with_ties <-
  data.frame(
    ID = fri_vals$ID,
    groupveg = groupveg_vals$GROUPVEG,
    pfs_replac = pfs_replac_vals$PRC_REPLAC,
    fri =  fri_vals$FRI_ALLFIR
  ) %>%
  # Replace -9999 with NAs
  mutate(replac_clean = ifelse(pfs_replac == -9999, NA, pfs_replac), fri_clean = ifelse(fri == -9999, NA, fri)) %>%
  mutate(
    frg = 
      case_when(
        # If vegetation group is not in allowed list, code as NA
        !(groupveg %in% allowed_groupveg) ~ NA_character_,
        # If FRI is NA and vegetation type is in allowed list, code as No Fire Regime
        is.na(fri_clean) & (groupveg %in% allowed_groupveg) ~ "No Fire Regime",
        # New FRG bins
        fri_clean <= 35 & replac_clean >= 66.7 ~ "High Frequency + High Severity",
        fri_clean <= 35 & replac_clean < 66.7 ~ "High Frequency + Low/Moderate Severity",
        fri_clean > 35 & fri_clean < 200 & replac_clean >= 66.7 ~ "Intermediate Frequency + High Severity",
        fri_clean > 35 & fri_clean < 200 & replac_clean < 66.7 ~ "Intermediate Frequency + Low/Moderate Severity",
        fri_clean >= 200 ~ "Infrequent Fire Regime",
        TRUE ~ NA_character_
      )
  ) %>%
  # Remove NA values before counting
  drop_na(frg) %>%
  # Count each fire regime group per ID
  count(ID, frg, name = "n") %>%
  group_by(ID) %>%
  # Get mode fire regime group per ID (keeping ties)
  slice_max(order_by = n, n = 1, with_ties = TRUE) %>%
  ungroup() 

# Identify IDs that have ties
frg_ties <-
  frg_mode_with_ties %>%
  add_count(ID, name = "n_top") %>%   # number of "winners" for that ID
  filter(n_top > 1) %>%
  arrange(ID, desc(n), frg)

tie_ids <- unique(frg_ties$ID)

# For tied IDs, use the original (unbuffered) coordinate value as the tiebreaker
if (length(tie_ids) > 0) {
  tie_coords <- coords_proj[tie_ids, ]

  pfs_replac_tie_vals <- terra::extract(pfs_replac, tie_coords, na.rm = FALSE)
  fri_tie_vals <- terra::extract(fri, tie_coords, na.rm = FALSE)

  tie_frg_vals <-
    data.frame(
      ID = tie_ids,
      pfs_replac = pfs_replac_tie_vals$PRC_REPLAC,
      fri = fri_tie_vals$FRI_ALLFIR
    ) %>%
    mutate(replac_clean = ifelse(pfs_replac == -9999, NA, pfs_replac), fri_clean = ifelse(fri == -9999, NA, fri)) %>%
    # No vegetation check here because these are all in allowed vegetation types since they were selected from the main coords
    mutate(
      frg = 
        case_when(
          is.na(fri_clean) ~ "No Fire Regime",
          fri_clean <= 35 & replac_clean >= 66.7 ~ "High Frequency + High Severity",
          fri_clean <= 35 & replac_clean < 66.7 ~ "High Frequency + Low/Moderate Severity",
          fri_clean > 35 & fri_clean < 200 & replac_clean >= 66.7 ~ "Intermediate Frequency + High Severity",
          fri_clean > 35 & fri_clean < 200 & replac_clean < 66.7 ~ "Intermediate Frequency + Low/Moderate Severity",
          fri_clean >= 200 ~ "Infrequent Fire Regime",
          TRUE ~ NA_character_
        ),
      ID = tie_ids
    ) %>%
    select(ID, frg)

  # Replace tied IDs in the buffer-mode result with point-based category
  frg_mode_resolved <-
    frg_mode_with_ties %>%
    filter(!(ID %in% tie_ids)) %>%          # drop tied IDs
    bind_rows(tie_frg_vals) %>%             # add back resolved values
    arrange(ID)
} else {
  frg_mode_resolved <- frg_mode_with_ties %>% arrange(ID)
}

stopifnot(nrow(frg_mode_resolved) == nrow(coords_proj))  # confirm number of rows matches number of coords

# Create dataframe with SampleID and fire regime group
frg_df <- 
  frg_mode_resolved %>%
  mutate(SampleID = coords_proj$SampleID) %>%
  select(SampleID, frg) %>%
  # Set factor levels
  mutate(frg = factor(frg, levels = c(
    "No Fire Regime",
    "Infrequent Fire Regime",
    "Intermediate Frequency + Low/Moderate Severity",
    "Intermediate Frequency + High Severity",
    "High Frequency + Low/Moderate Severity",
    "High Frequency + High Severity"
  )))

# Get recent fire data
# (Not same scale as evt or vegetation group, so just take median over buffer)
recent_fire_vals <-
  terra::extract(c(lffi, frq), coords_buffer, na.rm = FALSE, ID = TRUE) %>%
  group_by(ID) %>%
  mutate(lffi_nona = ifelse(is.na(lffi), 36, lffi)) %>%
  # Calculate LFFI in the buffer, ONLY using the pixels that experienced fire (i.e., not NA values)
  # Why: we want to know the fire free interval where fire did occur, not the average over the whole buffer which would be biased towards unburned areas
  # If you do it the other way (i.e., turn NAs into 39 first) you get biased high values
  # Either way, if you remove all the points whre lffi is NA, there is a significant relationship between lffi and Ho, so I am comfortable with this approach
  mutate(frq = ifelse(is.na(frq), 0, frq)) %>% # Treat NA as 0
  summarise(
    lffi = median(lffi, na.rm = TRUE),
    lffi_nona = median(lffi_nona, na.rm = TRUE),
    frq_median = median(frq, na.rm = TRUE),
    burned = any(frq > 0, na.rm = TRUE),
    burned_prop = mean(frq > 0, na.rm = TRUE),
    burned_sum = sum(frq > 0, na.rm = TRUE),
    frq_max = max(frq, na.rm = TRUE)
  ) %>%
  # Define recent fire based on whether the majority of the pixels in the buffer burned at least one time (>50% burned proportion)
  # Don't use frq_max for bining because that would code any point that had a single pixel burned as "recent fire" and that is too sensitive (I have tested both ways and bin(frq_max) is not significant)
  mutate(fire_recent = factor(ifelse(burned_prop > 0.5, 1, 0), levels = c(0,1), labels = c("Unburned", "Burned"))) %>%
  mutate(SampleID = coords_proj$SampleID) 

# Combine into a single dataframe
fire_mod_df <- 
  # Combine data frames
  groupveg_mode %>% 
  rename(groupveg = GROUPVEG) %>%    
  left_join(recent_fire_vals) %>%
  select(-ID) %>%
  # Add frg
  left_join(frg_df)  %>%
  # Create simple frequency column 
  mutate(fire_frq = case_when(
    grepl("High Frequency", frg) ~ "High Frequency",
    grepl("Intermediate Frequency", frg) ~ "Intermediate Frequency",
    frg == "Infrequent Fire Regime" ~ "No Fire/Low Frequency",
    frg == "No Fire Regime" ~ "No Fire/Low Frequency",
    TRUE ~ NA_character_
  )) %>%
  # Create factors 
  mutate(
    fire_frg = factor(frg, levels = rev(c("High Frequency + High Severity", "High Frequency + Low/Moderate Severity", "Intermediate Frequency + High Severity", "Intermediate Frequency + Low/Moderate Severity", "Infrequent Fire Regime", "No Fire Regime"))),
    fire_frq = factor(fire_frq, levels = rev(c("High Frequency", "Intermediate Frequency", "No Fire/Low Frequency"))),
    groupveg = as.character(groupveg),
    fire_recent = fire_recent
  )  %>%
  # Add vegetation type
  left_join(evt_mod_df %>% select(SampleID, evt), by = "SampleID") %>%
  # Put SampleID as first column
  select(SampleID, everything())

# Get counts of fire groups for methods
fire_mod_df %>%
  group_by(fire_frg) %>%
  count()
#   fire_frg                                           n
#   <fct>                                          <int>
# 1 No Fire Regime                                    14
# 2 Infrequent Fire Regime                            12
# 3 Intermediate Frequency + Low/Moderate Severity    10
# 4 Intermediate Frequency + High Severity            95
# 5 High Frequency + Low/Moderate Severity           115
# 6 High Frequency + High Severity                    11

fire_mod_df %>%
  group_by(fire_frq) %>%
  count()
#   fire_frq                   n
#   <fct>                  <int>
# 1 No Fire/Low Frequency     26
# 2 Intermediate Frequency   105
# 3 High Frequency           126

# Visualize EVT/Fire variables
source(here("analysis", "genetic_diversity", "functions_genetic_diversity.R"))
test_df <- 
  fire_mod_df %>%
  left_join(evt_mod_df, by = "SampleID") %>%
  left_join(get_het()) %>%
  mutate(fire_frq = factor(fire_frq, levels = c("No Fire/Low Frequency", "Intermediate Frequency", "High Frequency"), labels = c("No Fire/Low", "Intermediate", "High"))) %>%
  mutate(fire_recent = ifelse(burned_prop > 0.5, "Burned", "Unburned")) %>%
  mutate(fire_recent = factor(fire_recent, levels = c("Unburned", "Burned")))

pltA <-
  ggplot(test_df, aes(x = fire_frq, y = Ho)) +
  geom_boxplot(aes(fill = fire_recent)) +
  # jitter points grouped by recent fire
  geom_jitter(aes(group = fire_recent, size = frq_median, col = burned_prop), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), pch = 16, alpha = 0.8) +
  theme_classic() +
  scale_color_viridis_c(option = "plasma") +
  labs(
    fill = "Contemporary\nFire",
    size = "Median\nFire Count",
    color = "Proportion\nBurned",
    x = make_pretty_names("fire_frq"),
    y = "Observed Heterozygosity"
  ) +
  # Change order of legends 
  guides(
    fill = guide_legend(order = 1),
    size = guide_legend(order = 2),
    color = guide_colorbar(order = 3)
  ) +
  scale_fill_manual(values = c("white", "gray")) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    text = element_text(size = 18)
  )

pltB <-
  ggplot(test_df, aes(x = fire_frq, y = Ho)) +
  geom_boxplot(aes(fill = fire_recent)) +
  # jitter points grouped by recent fire
  geom_jitter(aes(group = fire_recent, color = evt), position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75), pch = 16, alpha = 1, size = 2) +
  theme_classic() +
  labs(
    fill = "Contemporary\nFire",
    color = "Vegetation\nType",
    x = make_pretty_names("fire_frq"),
    y = "Observed Heterozygosity"
  ) +
  scale_color_viridis_d(option = "turbo") +
  # Change order of legends 
  guides(
    fill = guide_legend(order = 1),
    size = guide_legend(order = 2),
    color = guide_legend(order = 3)
  ) +
  scale_fill_manual(values = c("white", "gray")) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    text = element_text(size = 18)
  )
png(here(plotpath,"recent_fire_plots.png"), width = 16, height = 5.5, units = "in", res = 300)
cowplot::plot_grid(pltA, pltB, nrow = 1, labels = c("A", "B"), align = "hv", label_size = 20)
dev.off()


# Pull out variables for modeling
fire_mod_df %>%
  select(
    SampleID, 
    # Historical fire frequency
    fire_frq, 
    # Contemporary fire
    fire_recent, 
    # Contemporary vegetation
    evt
  ) %>%
  write_csv(here(outpath, "landfire_data.csv"))

