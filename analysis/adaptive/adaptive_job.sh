# Run in analysis/adaptive

# Outstanding things to discuss:
#       - resamp vs just crop for future PCs?
#       - NA values in env_pres PCA -- IMPUTE?
#       - maps looking wonky?
#       - multiple time points - global offset?
#       - troubleshoot calculating AI from pc-corrected RDA results

# Generate env data files, including future layers; generate fig comparing present and future PCs
# Running the following is only necessary if environmental PCs are to be used
Rscript p0_env_data.R

# Generate vcf file with RDA outliers found within genes only; 
# alternatively, can just use genes.raw dosage matrix (skip to p2)
# unless Plink distances are needed for running GDM, which is also done here:
bash p1_process_snps.sh

# Run new RDA using RDA outliers found within genes
# Args are `gea_method`, `nPC`, `path`, and `prefix`
Rscript p2_rda_geagenes.R "bio1ndvi" "./outputs" "58-Sceloporus_bio1ndvi_gea"
# Rscript p2_rda_geagenes.R "bio1ndvi" "./outputs" "58-Sceloporus_bio1ndvi_gea_genes_nonsyn"

# Run GDM using RDA outliers and calculate offset
Rscript p3_gdm_gea.R "bio1ndvi" "./outputs" "58-Sceloporus_bio1ndvi_gea"

# Render the R Markdown document to visualize present and future environmental layers
Rscript -e "rmarkdown::render('env_viz.Rmd')"

# Render the R Markdown document for RDA results
Rscript -e "rmarkdown::render('RDA_genes.Rmd')"

# Render the R Markdown document for RDA adaptive index results within genes
Rscript -e "rmarkdown::render('adaptive_index.Rmd')"

# Render the R Markdown document for genomic offset calculations
Rscript -e "rmarkdown::render('genomic_offset.Rmd')"

