import ee

# Authenticate and initialize the Earth Engine API
ee.Authenticate()
ee.Initialize()

# Load California shape file
states_dataset = ee.FeatureCollection('TIGER/2018/States')
california = states_dataset.filter(ee.Filter.eq('NAME', 'California'))
roi = california.geometry()  # Use the geometry of the feature collection

# NDVI -----------------------------------------------------------------------------------------

# Define the date range
start_date = '2000-01-01'
end_date = '2020-12-31'

# Load MODIS NDVI data
modis_ndvi = ee.ImageCollection('MODIS/061/MOD13Q1') \
    .select('NDVI') \
    .filterDate(start_date, end_date)

# Function to clip each image in the collection to the ROI
def clip_image(image):
    return image.clip(roi)

# Apply the clipping function to each image in the collection
clipped_ndvi = modis_ndvi.map(clip_image)

# Calculate the mean NDVI over the specified period
mean_ndvi = clipped_ndvi.mean().multiply(0.0001)  # Scale factor for MODIS NDVI

# Export the mean NDVI layer to Google Drive
task = ee.batch.Export.image.toDrive(
    image=mean_ndvi,
    description='california_ndvi_mean_2000_2020',
    folder='EarthEngineExports',
    fileNamePrefix='california_ndvi_mean_2000_2020',
    region=roi,
    crs='EPSG:4326',
    maxPixels=1e13
)

# Start the export task
task.start()