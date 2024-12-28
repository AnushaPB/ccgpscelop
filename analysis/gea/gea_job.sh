# Download the required RDA files
bash p0_download_rda.sh

# Convert into a bed file that will be used to extract genes
Rscript p1_process_rda.R

# Extract genes using the bed file and the genome annotation
bash p2_get_genes.sh

# Intersect the extracted genes with the original RDA output
Rscript p3_intersect_genes.R

# Calculate stats from the intersected genes
bash p4_calc_gea_het.sh

# Render the R Markdown document for neutral adaptive test
Rscript -e "rmarkdown::render('gea_neutral_adaptive_test.Rmd')"

# Render the R Markdown document for gene ontology analysis
Rscript -e "rmarkdown::render('go_analysis.Rmd')"

# Render the R Markdown document for windowed nucleotide diversity (pi) analysis
Rscript -e "rmarkdown::render('window_pi.Rmd')"