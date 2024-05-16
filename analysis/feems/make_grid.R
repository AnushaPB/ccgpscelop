#source activate feems_e
Sys.setenv(PROJ_LIB = "/home/wanglab/miniconda3/envs/feems_e/share/proj")
Sys.setenv(PKG_CONFIG_PATH = "/home/wanglab/miniconda3/envs/feems_e/share/proj")

library("here")
library("dggridR")
library("sf")
library("ggplot2")
library("terra")
source(here("general_functions.R"))

#Generate a dggs specifying an intercell spacing of ~25 miles
dggs <- dgconstruct(res = 8, projection = "ISEA", aperture = 4, topology = "TRIANGLE")

#Read in ca border
ca <- get_ca()

#Write out shp file
shp_path <- here("analysis", "feems", "ca.shp")
sf::st_write(ca, shp_path, append = FALSE)

#Get a grid covering ca
grid <- dgshptogrid(dggs, shp_path)

#Plot borders and the associated grid
ggplot() +
    geom_sf(data=ca, fill=NA, color="black")   +
    geom_sf(data=grid, fill=alpha("blue", 0.4), color=alpha("white", 0.4))

# write out shp file
st_crs(grid) <- NA
sf::st_write(grid, here("analysis", "feems", "outputs", "triangle_res8.shp"), append = FALSE)
