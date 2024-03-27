# create a new conda environment
conda create -n=feems_e python=3.8.3 
conda activate feems_e

# install dependencies
conda install numpy==1.22.3 scipy==1.5.0 scikit-learn==0.23.1
conda install matplotlib==3.2.2 pyproj==2.6.1.post1 networkx==2.4.0 
conda install shapely==1.7.1 
conda install fiona
conda install pytest==5.4.3 pep8==1.7.1 flake8==3.8.3
conda install click==7.1.2 setuptools pandas-plink
conda install msprime==1.0.0 statsmodels==0.12.2 PyYAML==5.4.1
conda install xlrd==2.0.1 
conda install openpyxl==3.0.7
conda install suitesparse=5.7.2
conda install scikit-sparse=0.4.4 
conda install cartopy=0.18.0

# install feems
pip install git+https://github.com/NovembreLab/feems