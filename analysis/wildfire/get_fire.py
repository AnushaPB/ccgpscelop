import ee

# Authenticate and initialize the Earth Engine API
#ee.Authenticate()
ee.Initialize(project='ee-anushabishop')

# Define California region of interest (approximate bounding box) --------------
# (lon_min, lat_min, lon_max, lat_max)
roi = ee.Geometry.Rectangle([-125.0, 32.0, -113.0, 42.5])


# MTBS Burn Severity -----------------------------------------------------------
# Dataset: USFS/GTAC/MTBS/annual_burn_severity_mosaics/v1
# - Some images have band 'Severity', others have 'Burn_Severity'
# - There are images for AK and CONUS; we only want CONUS for CA

# Define the date range for MTBS annual mosaics
start_date = '1984-01-01'
end_date   = '2024-12-31'

# Load MTBS annual burn severity data (don't select bands yet)
mtbs_raw = ee.ImageCollection('USFS/GTAC/MTBS/annual_burn_severity_mosaics/v1') \
    .filterDate(start_date, end_date) \
    .filter(ee.Filter.stringContains('system:index', 'CONUS'))  # keep only CONUS

# Function to:
#   1) pick the correct severity band (Severity or Burn_Severity)
#   2) rename it to 'Severity'
#   3) clip to ROI
def harmonize_and_clip(image):
    band_names = image.bandNames()
    has_burn_sev = band_names.contains('Burn_Severity')

    severity_img = ee.Image(
        ee.Algorithms.If(
            has_burn_sev,
            image.select('Burn_Severity').rename('Severity'),
            image.select('Severity')
        )
    )

    # Clip to ROI and keep original time
    return severity_img.clip(roi).set('system:time_start', image.get('system:time_start'))

# Apply harmonization + clipping
mtbs_clipped = mtbs_raw.map(harmonize_and_clip)

# Sort by time to ensure chronological band order
mtbs_sorted = mtbs_clipped.sort('system:time_start')

# Build a multiband stack manually: one band per year (Severity_YYYY) ----------

# Convert the collection to a list so we can iterate
n_images = mtbs_sorted.size()
image_list = mtbs_sorted.toList(n_images)

# Initialize the stack with the first image
first = ee.Image(image_list.get(0))
first_year_str = ee.Date(first.get('system:time_start')).format('YYYY')
first_band_name = ee.String('Severity_').cat(first_year_str)
stack0 = first.select('Severity').rename(first_band_name)

# Function for iterate(): add each subsequent year as a new band
def add_band(index, prev_image):
    prev_image = ee.Image(prev_image)
    img = ee.Image(image_list.get(index))
    year_str = ee.Date(img.get('system:time_start')).format('YYYY')
    band_name = ee.String('Severity_').cat(year_str)
    # Each img has one band: 'Severity'
    band = img.select('Severity').rename(band_name)
    return prev_image.addBands(band)

# Use iterate to build the full stack
indices = ee.List.sequence(1, n_images.subtract(1))
severity_stack = ee.Image(indices.iterate(add_band, stack0))

# Export the multiband burn severity stack to Google Drive ---------------------

task = ee.batch.Export.image.toDrive(
    image=severity_stack,
    description='california_mtbs_burn_severity_stack_1984_2024',
    folder='EarthEngineExports',  # or None for root
    fileNamePrefix='california_mtbs_burn_severity_stack_1984_2024',
    region=roi,
    crs='EPSG:4326',
    scale=30,          # Landsat-resolution
    maxPixels=1e13
)

task.start()
print('Export started: MTBS burn severity stack for California (1984–2024)')