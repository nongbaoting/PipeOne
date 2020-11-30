conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

## lncRNA
conda env create  --name pipeOne_lncRNA --file ymls/pipeOne_lncRNA.yml
conda env create  --name pipeOne_py2 --file ymls/pipeOne_py2.yml
conda env create  --name pipeOne_py3 --file ymls/pipeOne_py3.yml

## circRNA
conda env create  --name pipeOne_CIRIquant --file ymls/pipeOne_CIRIquant.yml

conda env create  --name pipeOne_apa --file ymls/pipeOne_apa.yml
## install qapa
cd qapa-1.3.0
conda activate pipeOne_apa
python setup.py install
conda deactivate
cd ..

## retro
conda env create  --name pipeOne_telescope --file ymls/pipeOne_telescope.yml

## fusion
conda create -y -n pipeOne_fusion python=3 STAR=2.7.3a arriba fire

## editing
conda env create  --name pipeOne_RnaEditing --file ymls/pipeOne_RnaEditing.yml
###  install sprint
git clone https://github.com/jumphone/SPRINT.git
cd SPRINT
conda activate pipeOne_RnaEditing
python setup.py install
cd ..

## install AS 
conda env create  --name pipeOne_AS --file ymls/pipeOne_AS.yml
conda activate pipeOne_AS
pip install spladder
conda activate base

## gatk
conda env create  --name pipeOne_gatk3.8 --file ymls/pipeOne_gatk3.8.yml
conda activate pipeOne_gatk3.8
gatk3-register GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar
conda activate base

## machine learning
conda env create  --name pipeOne_ml --file ymls/pipeOne_ml.yml
R -e 'install.packages(c("survival", "survminer", "ggplot2", "tidyverse", "data.table"))'


## packages we have download and include in the current directory
# git clone https://github.com/jumphone/SPRINT.git

# wget http://www.rnabinding.com/CPPred/CPPred/CPPred.tar.gz
# tar -xvf CPPred.tar.gz

## install CPAT
# wget  https://sourceforge.net/projects/rna-cpat/files/v1.2.4/CPAT-1.2.4.tar.gz
# tar -xvf CPAT-1.2.4.tar.gz

# wget https://github.com/morrislab/qapa/archive/v1.3.0.tar.gz
# tar -xvf v1.3.0.tar.gz

# wget https://console.cloud.google.com/storage/browser/_details/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
# tar -xvf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

#  rm -rf  *gz *bz2

# R