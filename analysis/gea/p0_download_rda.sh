# Run in analysis/gea/
# RDA outputs:
mkdir -p outputs
#rsync -avz --include '*/' --include '*rdadapt.csv*' --exclude '*' hgdownload.soe.ucsc.edu::ccgp/algatr/58-Sceloporus/subsets/RDA/ outputs/

rsync --list-only hgdownload.soe.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA/

rsync -avz hgdownload.soe.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA/58-Sceloporus_RDA_peakRAM.csv outputs/
rsync -avz hgdownload.soe.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA/58-Sceloporus_RDA_outliers_full_Zscores.csv outputs/
rsync -avz hgdownload.soe.ucsc.edu::ccgp/algatr/58-Sceloporus/RDA/58-Sceloporus_RDA_outliers_full_rdadapt.csv outputs/