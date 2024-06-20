#!/bin/bash
#SBATCH -N 1
#SBATCH -c 16
#SBATCH -p kellis
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ruiwenfu@mit.edu
#SBATCH --job-name=SCENIC
module load miniconda3

for cellType in 'cDC1' 'cDC2' 'cDC3'	; do

  source activate r4
  Rscript /net/bmc-lab5/data/kellis/users/ruiwenfu/SCENIC/bashScripts/scenic_scRNA_cellType_part1Coexpression.r $cellType
  echo 'scenic part1 done'
  conda deactivate

  source activate /net/bmc-lab5/data/kellis/users/ruiwenfu/conda_related/envs/python3.9
  python /net/bmc-lab5/data/kellis/group/Takeda_final/Takeda_single_cell/GRNBoost2.py '/net/bmc-lab5/data/kellis/users/ruiwenfu/SCENIC/'$cellType'/'
  echo 'GRNBoost done'
  conda deactivate

  source activate r4
  Rscript /net/bmc-lab5/data/kellis/users/ruiwenfu/SCENIC/bashScripts/scenic_scRNA_cellType_part2scoreGRN.r $cellType
  echo 'scenic part2 done'
  conda deactivate
done
