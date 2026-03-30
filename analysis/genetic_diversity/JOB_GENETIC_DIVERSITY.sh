# Calculate genetic diversity statistics
bash p1_genetic_diversity.sh

# Pre-process fire data for modelling
Rscript p2_process_fire_data.R

# Render the R Markdown document for genetic diversity results
Rscript -e "rmarkdown::render('p3_notebook_genetic_diversity_model.Rmd')"

# Render the R Markdown document for non-synonymous genetic diversity
Rscript -e "rmarkdown::render('p4_notebook_nonsyn_diversity.Rmd')"

# Render the R Markdown document for maps of genome-wide and non-synonymous diversity
Rscript -e "rmarkdown::render('p5_notebook_diversity_maps.Rmd')"