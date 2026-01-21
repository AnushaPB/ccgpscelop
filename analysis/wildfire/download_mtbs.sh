# TEMPORARY LINK: https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/download-tool/Fire_data_bundles_AExAOW3fglZELZrjm3bZ.zip generated from https://www.mtbs.gov/direct-download
wget -c -O Fire_data_bundles_AExAOW3fglZELZrjm3bZ.zip 'https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/download-tool/Fire_data_bundles_AExAOW3fglZELZrjm3bZ.zip'
unzip Fire_data_bundles_AExAOW3fglZELZrjm3bZ.zip -d mtbs_fire
rm Fire_data_bundles_AExAOW3fglZELZrjm3bZ.zip

# Unzip nested files
for z in mtbs_fire/composite_data/MTBS_BSmosaics/*/mtbs_CA_*.zip; do
  unzip -o "$z" -d "$(dirname "$z")"
done

# Move all files from year subfolders into MTBS_Fire/
find mtbs_fire/composite_data/MTBS_BSmosaics -type f -exec mv -t mtbs_fire {} +

# Remove composite_data folder
rm -rf mtbs_fire/composite_data



# Background/No Data - (Value/code = 0) Areas outside of the burned area boundary and excluded from the MTBS mapping process.
# Increased Greenness - (Value/code = 5) Areas that burned but display more vegetation cover, density, and/or productivity (vigor), usually within one growing season after fire. This is a fire-caused effect from release of nutrients into soil, and/or reduced competition for nutrients, light and water (much like a thinning effect). These areas are usually herbaceous or low shrub communities that undergo little change in species composition after fire.
# Unburned to Low - (Value/code = 1) Areas that are either unburned, or when visible fire effects occupy a small proportion of the site, on the order of less than 5 percent. If more of the site is burned, then effects are limited to a few biophysical components. The class may also include areas that recover very quickly after fire, such as grasslands or light surface burns under dense, non-impacted forest canopies.
# Low - (Value/code = 2) Areas where more than a small proportion of the site burned. Collectively, all strata are slightly altered from the pre-fire state. Duff, woody debris and newly exposed mineral soil typically exhibit some change. Low vegetation (<1 meter) and shrubs or trees (1-5 meters) may show significant aboveground scorch, char or consumption, and vegetation density or cover may be greatly altered. These prefire plants are generally still viable and recover quickly (within a year or two), with little change in species composition. An exception is western conifers, where sapling-sized trees may exhibit 50 percent or more mortality. Intermediate and large overstory trees may exhibit up to 25 percent mortality evidenced by crown char or scorch. Where charring does not kill tree crowns, as is common in the southeast, higher percentages of black char may occur. Char height from ground flames is typically less than 3 meters.
# Moderate - (Value/code = 3) The moderate class is difficult if not impossible to briefly describe. Indicators may be fairly consistent across biophysical strata and will exhibit traits between the low and high severity classes. On the other hand, numerous potential combinations of distinct low and high indicators may occur to yield a moderate classification overall within the minimum mapping unit. Conditions are transitional in magnitude and/or uniformity between the low and high characteristics described.
# High - (Value/code = 4) This class is characterized by fairly consistent effects across a site. In forested ecosystems, litter is totally consumed; duff is typically nearly entirely consumed. Medium and heavy woody debris are at least partially consumed and at least deeply charred with mostly ash and charcoal remaining. Overstory trees typically exhibit greater than 75 percent mortality. Biomass consumption and above-ground changes in carbon balances are significant. Crown char is frequently 100 percent from torching fire, and significant branch loss is evident at the highest crown levels. Where crown torching did not occur, char height from ground flames often exceeds 4 meters. Overstory tree effects are generally long lasting. New tree establishment may occur 1-3 years post-fire, but forest development often takes many decades. Herbaceous plants and shrubs are almost completely charred or consumed above ground, often with notable branch loss on taller shrubs, which may be reduced to small stubs. Resprouting from perennial plants, except grasses, is strongly reduced, as most individuals lose viability with a significant reduction in cover.
# Non-Processing Mask - (Value/Code = 6). Areas within the burned area boundary representing missing data due to image sensor problems (Landsat 7 scan line omissions) or atmospheric/terrain interference (clouds, smoke, shadow, snow). No attmepts are made to fill in missing data areas through interpolation or other methods, though MTBS data users may consider this for a given analysis objective.
