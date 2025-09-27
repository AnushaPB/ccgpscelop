# Run in analysis/adaptive

# Generate vcf file with RDA outliers found within genes only and generate Plink distances used for running GDM
bash p1_process_snps.sh

# Run GDM using RDA outliers and calculate offset
# Args are `path` and `prefix`
Rscript p2_gdm_gea.R "./outputs" "58-Sceloporus_bio1ndvi_gea"

# Render the R Markdown document for GDM-based adaptive index results
Rscript -e "rmarkdown::render('adaptive_index.Rmd')"

# Render the R Markdown document for genomic offset
Rscript -e "rmarkdown::render('genomic_offset.Rmd')"

# Render the R Markdown document for building adaptive plots
Rscript -e "rmarkdown::render('adaptive_plots.Rmd')"
