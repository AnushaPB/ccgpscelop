# Outstanding things to discuss:
#       - resamp vs just crop for future PCs?
#       - NA values in env_pres PCA -- IMPUTE?
#       - maps looking wonky?
#       - multiple time points - global offset?
#       - troubleshoot calculating AI from pc-corrected RDA results

# Generate env data files, including future layers; generate fig comparing present and future PCs
Rscript p0_env_data.R

# Generate vcf file with RDA outliers found within genes only; alternatively, can just use genes.raw dosage matrix (skip to p2)
bash p1_process_snps.sh

# Run new RDA using RDA outliers found within genes
# Args are `gea_method` and `nPC`
Rscript p2_rda_geagenes.R "bio1ndvi"

# Render the R Markdown document to visualize present and future environmental layers
Rscript -e "rmarkdown::render('env_viz.Rmd')"

# Render the R Markdown document for RDA results
Rscript -e "rmarkdown::render('RDA_genes.Rmd')"

# Render the R Markdown document for RDA adaptive index results within genes
Rscript -e "rmarkdown::render('adaptive_index.Rmd')"

# Render the R Markdown document for genomic offset calculations
Rscript -e "rmarkdown::render('genomic_offset.Rmd')"

