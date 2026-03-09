# Create a new conda environment
conda create -n=env_feems python=3.8.3 
conda activate env_feems

# Install dependenciesconda install -y numpy==1.22.3 scipy==1.5.0 scikit-learn==0.23.1
conda install -y matplotlib==3.2.2 pyproj==2.6.1.post1 networkx==2.4.0 
conda install -y shapely==1.7.1 
conda install -y fiona
conda install -y pytest==5.4.3 pep8==1.7.1 flake8==3.8.3
conda install -y click==7.1.2 setuptools pandas-plink
conda install -y msprime==1.0.0 statsmodels==0.12.2 PyYAML==5.4.1
conda install -y xlrd==2.0.1 
conda install -y openpyxl==3.0.7
conda install -y suitesparse=5.7.2
conda install -y scikit-sparse=0.4.4 
conda install -y cartopy=0.18.0

# Install feems
pip install git+https://github.com/NovembreLab/feems
