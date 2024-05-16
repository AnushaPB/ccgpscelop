# Run in analysis/feems directory
# Set up feems environment
bash feems_run.sh

# Create grid for feems
Rscript make_grid.R

# Visualize results 
Rscript -e "rmarkdown::render(here::here("analysis", "feems", "feems_results.Rmd"))"