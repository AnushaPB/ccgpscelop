# Calculate genetic distances
bash p1_gendist.sh

# Render the R Markdown document for the IBD/IBE analysis
Rscript -e "rmarkdown::render('p2_notebook_ibdibe.Rmd')"