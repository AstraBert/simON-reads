#!/bin/bash

echo "Installation process started, checking if Anaconda/Miniconda and Mamba are installed..."


command -v conda >/dev/null 2>&1 || {
   echo >&2 "The installation pipeline requires Anaconda/Miniconda but it is not installed. Please check here: https://anaconda.org/ for more details. Aborting."
   exit 1
}


command -v mamba >/dev/null 2>&1 || {
   echo >&2 "The installation pipeline requires Mamba but it is not installed. Please check here: https://github.com/conda-forge/miniforge#mambaforge for more details. Aborting."
   exit 2
}

echo "Anaconda/Miniconda and Mamba are correctly installed, proceeding with the installation of necessary softwares: it may take a while..."

sleep 1

tmp=$(dirname $0)
basedir=$(realpath $tmp)

ConPath=$(which conda)
tmp=${ConPath#* }
Conda=${tmp%%/bin/co*}

# CREATE A NEW FOLDER WHERE TO STORE THE ENVIRONMENT
mkdir -p $basedir/environments

# CREATE A MAMBA/CONDA ENVIRONMENT
mamba create \
  -p $basedir/environments/simON-reads  \
  -y \
  -c conda-forge \
  -c bioconda \
  python=3.10

source ${Conda}/etc/profile.d/conda.sh
conda activate \
  $basedir/environments/simON-reads

## INSTALL necessary dependencies
mamba install \
   -c conda-forge \
   -y \
   matplotlib 

mamba install \
   -c conda-forge \
   -y \
   biopython 

mamba install \
   -c conda-forge \
   -y \
   dos2unix

## DEACTIVATE
conda deactivate

# make exectuable and change the interpreter to the python installation in the environment
cp -n ${basedir}/simON_reads_module.py ${basedir}/environments/simON-reads/bin/simON_reads_module.py
echo \#\!${basedir}/environments/simON-reads/bin/python3.10 > ${basedir}/environments/simON-reads/bin/simON_reads.py
awk 'NR>1' ${basedir}/simON_reads.py >>${basedir}/environments/simON-reads/bin/simON_reads.py
dos2unix ${basedir}/environments/simON-reads/bin/simON_reads.py
chmod 755 ${basedir}/environments/simON-reads/bin/simON_reads.py

echo "Installation process is now completed"
