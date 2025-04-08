# Run in analysis/gea/
# RDA outputs:
mkdir -p outputs
rsync -avz hgdownload.soe.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA_no_pca/58-Sceloporus_RDA_outliers_full_Zscores.csv ../analysis/gea/outputs
rsync -avz hgdownload.soe.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA_no_pca/58-Sceloporus_RDA_cortest_full.csv  ../analysis/gea/outputs