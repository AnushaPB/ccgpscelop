PLINK=../../data/ccgp_data/58-Sceloporus_complete_coords_annotated

plink2 --bfile $PLINK \
  --allow-extra-chr \
  --extract range outputs/melanization_genes.bed \
  --recode A  \
  --out outputs/melanization_genes

# RUN PCA on plink file
plink2 --pca 3 --bfile $PLINK --out outputs/pca --allow-extra-chr

