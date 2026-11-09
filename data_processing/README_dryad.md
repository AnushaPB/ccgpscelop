# Data from 'Eco-evolutionary processes mediate the effects of climate change and wildfire on genetic diversity'

## Dataset description

This dataset accompanies the manuscript "Eco-evolutionary processes mediate the effects of climate change and wildfire on genetic diversity" and contains whole-genome data from 257 western fence lizards (*Sceloporus occidentalis*) sampled across California as part of the California Conservation Genomics Project (CCGP).

It includes genetic data, sampling coordinates, callable genomic regions, and the reference genome and annotation used in the study.

These files provide the genomic data used in the main analyses presented in the manuscript. For the code used for these analyses see the corresponding GitHub repository [AnushaPB/ccgpscelop](https://github.com/anushapb/ccgpscelop/tree/main/analysis) which is archived on Zenodo ([![DOI](https://zenodo.org/badge/696150420.svg)](https://doi.org/10.5281/zenodo.24091513)).

## File descriptions

| File | Description |
|---|---|
| `58-Sceloporus_complete_coords_annotated_final_samples.vcf.gz` | Full VCF containing quality-filtered genotypes for the 257 individuals retained in the study. See the manuscript for details on bioinformatic processing and filtering. |
| `58-Sceloporus_complete_coords_annotated_chr.bed`, `58-Sceloporus_complete_coords_annotated_chr.bim`, `58-Sceloporus_complete_coords_annotated_chr.fam` | PLINK dataset containing genotypes for the 10 autosomes. |
| `58-Sceloporus_annotated_pruned_0.6_chr.bed`, `58-Sceloporus_annotated_pruned_0.6_chr.bim`, `58-Sceloporus_annotated_pruned_0.6_chr.fam` | LD-pruned PLINK dataset for the 10 autosomes. |
| `58-Sceloporus_pruned_0.6_thinned_10kb_chr.bed`, `58-Sceloporus_pruned_0.6_thinned_10kb_chr.bim`, `58-Sceloporus_pruned_0.6_thinned_10kb_chr.fam` | LD-pruned PLINK dataset for the 10 autosomes, further thinned to a minimum spacing of 10 kb between variants. |
| `58-Sceloporus.idepth` | Individual-level sequencing depth statistics. |
| `58-Sceloporus.imiss` | Individual-level genotype missingness statistics. |
| `58-Sceloporus.coords.txt` | Sampling coordinates. Column 1: sample ID; Column 2: longitude (decimal degrees); Column 3: latitude (decimal degrees). |
| `58-Sceloporus_callable_sites.bed` | Callable genomic regions. Columns: chromosome/scaffold, start position, end position (base pairs). |
| `rSceOcc1.20250428.p_ctg_chr_names.fna` | *S. occidentalis* reference genome used in the study. |
| `rSceOcc1.20250428.p_ctg_chr_names.fna.fai` | FASTA index for the reference genome. |
| `complete.genomic.gff` | Structural annotation of the reference genome, including genes and associated genomic features. |