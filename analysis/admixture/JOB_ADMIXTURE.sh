# Run admixture
bash p1_run_admixture.sh

# Run PCA
bash p2_run_pca.sh

# Render the R Markdown document for admixture results
Rscript -e "rmarkdown::render('p3_notebook_admixture.Rmd')"