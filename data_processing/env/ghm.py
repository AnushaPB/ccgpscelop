import ee

# Authenticate and initialize the Earth Engine API
ee.Authenticate()
ee.Initialize()

# Load California shape file
states_dataset = ee.FeatureCollection('TIGER/2018/States')
california = states_dataset.filter(ee.Filter.eq('NAME', 'California'))
roi = california.geometry()  # Use the geometry of the feature collection

# GLOBAL HUMAN MODIFICATION --------------------------------------------------------------------
# Load the Global Human Modification dataset
gHM = ee.ImageCollection('CSP/HM/GlobalHumanModification')

# Function to clip each image in the collection to the ROI
def clip_image(image):
    return image.clip(roi)

# Apply the clipping function to each image in the collection
clipped_gHM = gHM.map(clip_image)

# Create a mosaic of the clipped images to get a single image
mosaic_clipped_gHM = clipped_gHM.mosaic()
mosaic_gHM = gHM.mosaic()

# Export the 'Human modification' layer to Google Drive
task = ee.batch.Export.image.toDrive(
    image=mosaic_clipped_gHM,
    description='california_global_human_modification_2016_1km',
    folder='EarthEngineExports',
    fileNamePrefix='california_global_human_modification_2016_1km',
    region=roi,
    scale=1000,
    crs='EPSG:4326',
    maxPixels=1e13
)

# Start the tasks
task.start()
