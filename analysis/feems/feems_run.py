# source activate feems_env

# base
import numpy as np
import pkg_resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink

# viz
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz
import pandas as pd

# grid
import geopandas as gpd
from shapely.geometry import Polygon

# set data path
data_path = pkg_resources.resource_filename("feems", "data/")

# change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"

# get working direcotyr in python
import os
os.getcwd()

# read the genotype data and mean impute missing data
scelop_data_path = "../../58-Sceloporus"
plink_path = "{}/58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_r60".format(scelop_data_path)
(bim, fam, G) = read_plink(plink_path)
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)
print("n_samples={}, n_snps={}".format(genotypes.shape[0], genotypes.shape[1]))

# setup graph
# Read the coordinates from the CSV file
coords_df = pd.read_csv("../../data/58-Sceloporus.coords.txt", delimiter="\t", header=None)
# Filter coords_df by SampleID using ids from the plink file
sample_ids = fam.iid.tolist()
filtered_coords_df = coords_df[coords_df[0].isin(sample_ids)]
# Extract the x and y columns from the filtered coords_df
coord = filtered_coords_df.iloc[:, [1, 2]]
# convert to numpy array
coord = coord.values

# setup graph
coordd = np.loadtxt("{}/wolvesadmix.coord".format(data_path))  # sample coordinates
outerd = np.loadtxt("{}/wolvesadmix.outer".format(data_path)) 

# Calculate the outer coordinates by getting the minimum and maximum values from coord
min_x = np.min(coord[:, 0])
max_x = np.max(coord[:, 0])
min_y = np.min(coord[:, 1])
max_y = np.max(coord[:, 1])

# Add a small buffer to the outer coordinates
buffer = 0.01  # Adjust the buffer size as needed
outer = np.array([[min_x - buffer, min_y - buffer],
          [min_x - buffer, max_y + buffer],
          [max_x + buffer, max_y + buffer],
          [max_x + buffer, min_y - buffer]])

# Continue with the rest of the code
grid_pathd = "{}/grid_100.shp".format(data_path)  # path to discrete global grid
grid_path = "triangle_res8.shp"

# Load the shapefile
#gdf = gpd.read_file(grid_path)

# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                       ggrid=grid_path,
                       translated=True, 
                       buffer=0,
                       outer=outer)

# set up the spatial graph
sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)

# plot coords/graph
projection = ccrs.EquidistantConic(central_longitude=-108.842926, central_latitude=66.037547)
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
        edge_alpha=1, edge_zorder=100, sample_pt_size=10, 
        obs_node_size=7.5, sample_pt_color="black", 
        cbar_font_size=10)
v.draw_map()
v.draw_samples()
v.draw_edges(use_weights=False)
v.draw_obs_nodes(use_ids=False)

# fit
sp_graph.fit(lamb = 20.0)


fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)  
v = Viz(ax, sp_graph, projection=projection, edge_width=.5, 
  edge_alpha=1, edge_zorder=100, sample_pt_size=20, 
  obs_node_size=7.5, sample_pt_color="black", 
  cbar_font_size=10)
v.draw_map()
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False) 
v.draw_edge_colorbar()

# Save the figure
plt.savefig('feems.png')

# If you want to show the plot as well, uncomment the next line
# plt.show()

# EXPORT

# Export edges
weights = sp_graph.w
edges = sp_graph.edges

# Create a LineString for each edge
from shapely.geometry import LineString
geometry = [LineString([sp_graph.node_pos[edge[0]], sp_graph.node_pos[edge[1]]]) for edge in edges]

# Create a DataFrame with the weights
df = pd.DataFrame(weights, columns=['weight'])

# Create a GeoDataFrame with the edges and weights
gdf = gpd.GeoDataFrame(df, geometry=geometry)

# Write out shp file
gdf.to_file("feems_edges.shp", layer = "weight")

# Extract sample size for each node
from feems.spatial_graph import query_node_attributes
permuted_idx = query_node_attributes(sp_graph, "permuted_idx") 
obs_perm_ids = permuted_idx[: sp_graph.n_observed_nodes]
node_positions = grid[obs_perm_ids, :]
node_sizes = sp_graph.n_samples_per_obs_node_permuted

# Create a Point for each node
from shapely.geometry import Point
geometry = [Point(pos) for pos in node_positions]

# Create a DataFrame with the node sizes
df = pd.DataFrame({'size': node_sizes})

# Create a GeoDataFrame with the nodes and sizes
gdf = gpd.GeoDataFrame(df, geometry=geometry)

# Write out shapefile
gdf.to_file("feems_nodes.shp", layer = "weight")

