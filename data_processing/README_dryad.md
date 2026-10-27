# Data from `The genetic consequences of climate change and wildfire depend on evolutionary context`

## Dataset description

This dataset accompanies the manuscript "The genetic consequences of climate change and wildfire depend on evolutionary context."

The dataset contains whole-genome genetic data from 257 western fence lizards (*Sceloporus occidentalis*) sampled across California as part of the California Conservation Genomics Project (CCGP). It includes genetic data, sampling coordinates, callable genomic regions, sequencing quality statistics, and the reference genome and annotation used in the study.

These files provide the genomic data necessary to reproduce the main analyses presented in the manuscript, using the analysis scripts available in the `analysis/` folder of the corresponding GitHub repository AnushaPB/ccgpscelop (https://github.com/anushapb/ccgpscelop/tree/main/analysis), which is also archived on Zenodo.

## File descriptions

| File | Description |
|---|---|
| `58-Sceloporus_complete_coords_annotated_final_samples.vcf.gz` | Filtered VCF containing genome-wide variants and genotypes for the 257 individuals retained in the study. |
| `58-Sceloporus_complete_coords_annotated_chr.{bed,bim,fam}` | PLINK dataset used for genetic diversity and genotype-environment association analyses. |
| `58-Sceloporus_annotated_pruned_0.6_chr.{bed,bim,fam}` | LD-pruned PLINK dataset used for PCA and genetic distance analyses. |
| `58-Sceloporus_pruned_0.6_thinned_10kb_chr.{bed,bim,fam}` | LD-pruned and thinned PLINK dataset used for ADMIXTURE analyses. |
| `58-Sceloporus.coords.txt` | Sampling coordinates. Column 1: sample ID; Column 2: longitude (decimal degrees); Column 3: latitude (decimal degrees). |
| `58-Sceloporus_callable_sites.bed` | Callable genomic regions. Columns: chromosome/scaffold, start position, end position (base pairs). |
| `58-Sceloporus.idepth` | Individual sequencing depth (VCFtools output). Columns: INDV (sample ID), N_SITES (number of sites), MEAN_DEPTH (mean sequencing depth). |
| `58-Sceloporus.imiss` | Individual genotype missingness (VCFtools output). Columns: INDV (sample ID), N_DATA (sites evaluated), N_GENOTYPES_FILTERED (filtered genotypes), N_MISS (missing genotypes), F_MISS (proportion missing). |
| `rSceOcc1.20250428.p_ctg_chr_names.fna` | Revised *S. occidentalis* reference genome used in the study. |
| `rSceOcc1.20250428.p_ctg_chr_names.fna.fai` | FASTA index for the revised reference genome. |
| `complete.genomic.gff` | Structural annotation of the revised reference genome, including genes and associated genomic features. |