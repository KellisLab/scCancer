#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p kellis
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ruiwenfu@mit.edu
#SBATCH --job-name=LIana
module load miniconda3

source activate r4

Rscript /net/bmc-lab5/data/kellis/users/ruiwenfu/scRNA/metastatic_all/cpdb/LIana.R
echo 'LIana done'
