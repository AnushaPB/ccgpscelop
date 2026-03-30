# Download the required RDA files
bash p0_download_rda.sh

# Convert into a bed file that will be used to extract genes
Rscript p1_process_rda.R

# Calculate and process correlations between SNPs
bash p2_calc_correlations.sh
Rscript p3_process_correlations.R

# Extract genes using the bed file and the genome annotation
bash p4_get_genes.sh

# Intersect the extracted genes with the original RDA output
Rscript p5_intersect_genes.R

# Intersect SNPs with consequences (e.g., missense, synonymous) from the genome annotation
Rscript p6_csq_intersect.R

# Calculate statistics and make basic plots
bash p7_calc_stats.sh
Rscript p8_summary_stats_and_plots.R

# Calculate windowed pi
bash p9_window_pi.sh

# Render the R Markdown document for genes of interest
Rscript -e "rmarkdown::render('p10_notebook_genes_of_interest.Rmd')"

# Render the R Markdown document for windowed pi
Rscript -e "rmarkdown::render('p11_notebook_window_pi.Rmd')"