# Create a new conda environment
conda create -n=feems_env
conda activate feems_env
conda install pandas pandas-plink geopandas
conda install numpy==1.22.3
conda install -c bioconda feems -c conda-forge

# Export the environment
conda env export --name feems_env > feems_env.yml
#conda env create --name feems_env --file feems_env.yml

#https://github.com/NovembreLab/feems/issues/15
conda create --name feems3 --file feems3.txt
pip install git+https://github.com/NovembreLab/feems
conda install -c conda-forge proj

conda create --name feems4 --file erik_feems.txt
conda activate feems4
conda install -c conda-forge proj gdal
export PROJ_LIB=$(conda info --base)/envs/feems4/share/proj


conda create --name feems6 --file erik_feems.txt
pip install git+https://github.com/NovembreLab/feems
conda install -c conda-forge proj gdal