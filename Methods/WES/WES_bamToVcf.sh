#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p kellis
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ruiwenfu@mit.edu
#SBATCH --job-name=updateVCF


module load bcftools

BAM=/net/bmc-lab5/data/kellis/group/031621-WES/RP-1381_1000N_v1_Exome_OnPrem.hg38.bam

SAMPLE=$(basename $BAM)
SAMPLE=${SAMPLE%%.*}

VCF=/net/bmc-lab5/data/kellis/group/031621-WES/prospective_hg38_variant/"$SAMPLE".vcf
SELECTED=/net/bmc-lab5/data/kellis/group/031621-WES/prospective_hg38_variant/selected/"$SAMPLE"Selected.vcf
REHEADED=/net/bmc-lab5/data/kellis/group/031621-WES/prospective_hg38_variant/selected/"$SAMPLE"SelectedReheaded.vcf


bcftools mpileup -Ob -I -f /net/bmc-lab5/data/kellis/group/011121-RNAseq/ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta $BAM -R selectChr.txt\
| bcftools call -mv -Ov -o $VCF

grep -w '^#\|^#CHROM\|^chr[1-9]\|chr[1-2][0-9]\|chr[X]' $VCF > $SELECTED
bcftools reheader -f broad_hg38_v0_selected.fasta.fai $SELECTED -o $REHEADED


echo 'done'
