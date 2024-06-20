#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p kellis
#SBATCH --array=0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ruiwenfu@mit.edu
#SBATCH --job-name=velocyto
module load samtools
module load miniconda3
source activate /net/bmc-lab5/data/kellis/users/ruiwenfu/conda_related/envs/cellrank

SAMPLE_LIST=($(<"scCancer_10x_path.txt"))
FOLDER=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}
echo $FOLDER
SAMPLE=$(basename $FOLDER)
echo $SAMPLE
OUT='/net/bmc-lab5/data/kellis/group/Fu_Doris/velocyto/'
GTF='/net/bmc-lab5/data/kellis/group/Fu_Doris/velocyto/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
MASK='/net/bmc-lab5/data/kellis/group/Fu_Doris/velocyto/hg38_rmsk.gtf'

##### make directory and copy 10x files
OUT=$OUT$SAMPLE/

echo $OUT
# mkdir -p $OUT

cp $FOLDER'/outs/possorted_genome_bam.bam' $OUT
cp $FOLDER'/outs/possorted_genome_bam.bam.bai' $OUT
cp $FOLDER'/outs/filtered_feature_bc_matrix/barcodes.tsv.gz' $OUT

BARCODE=$OUT'barcodes.tsv.gz'
BAM=$OUT'possorted_genome_bam.bam'

echo 'finished copying'
# run velocyto
velocyto run -b $BARCODE  -o $OUT -m $MASK $BAM $GTF

echo 'done'
