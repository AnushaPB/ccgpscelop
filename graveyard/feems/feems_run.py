# first set-up feems with: feems_setup.sh
# conda activate feems_e
# cd analysis/feems
# IMPORTANT: do not name this file feems.py because then it won't recognize the package
# base
import numpy as np
import pkg_resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink

# viz
import matplotlib.pyplot as plt
# causes core dump:
import cartopy.crs as ccrs

# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz
import pandas as pd

# grid
#import geopandas as gpd
from shapely.geometry import Polygon

# change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"

# get working directory in python
import os
os.getcwd()

# read the genotype data and mean impute missing data
scelop_data_path = "../../data/processed_data"
plink_path = "{}/58-Sceloporus_maf05_minDP5_maxDP50_rmsamp60_mm80_rmsamp20_r60".format(scelop_data_path)
(bim, fam, G) = read_plink(plink_path)
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)
print("n_samples={}, n_snps={}".format(genotypes.shape[0], genotypes.shape[1]))

# setup graph
# Read the coordinates from the CSV file
coords_df = pd.read_csv("../../data/raw_data/58-Sceloporus.coords.txt", delimiter="\t", header=None)
# Filter coords_df by SampleID using ids from the plink file
sample_ids = fam.iid.tolist()
filtered_coords_df = coords_df[coords_df[0].isin(sample_ids)]
# Extract the x and y columns from the filtered coords_df
coord = filtered_coords_df.iloc[:, [1, 2]]
# convert to numpy array
coord = coord.values

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
grid_path = "outputs/triangle_res8.shp"

# Load the shapefile (if you want to plot it)
#gdf = gpd.read_file(grid_path)

# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                       ggrid=grid_path,
                       translated=True, 
                       buffer=0,
                       outer=outer)

# set up the spatial graph
sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)

# EXPORT --------------------------------

# Code from this issue: https://github.com/NovembreLab/feems/issues/34
feems_nodes = sp_graph.nodes
pd.DataFrame(feems_nodes).to_csv('outputs/feems_nodes.csv', header=False, index=False)

feems_node_pos = np.vstack((sp_graph.node_pos.T,[sp_graph.nodes[n]['n_samples'] for n in range(len(sp_graph.nodes))])).T
pd.DataFrame(feems_node_pos).to_csv('outputs/feems_node_pos.csv', header=False, index=False)

feems_edges = sp_graph.edges
pd.DataFrame(feems_edges).to_csv('outputs/feems_edges.csv', header=False, index=False)

feems_w = sp_graph.w
pd.DataFrame(feems_w).to_csv('outputs/feems_w.csv', header=False, index=False)

# LOT ---

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


