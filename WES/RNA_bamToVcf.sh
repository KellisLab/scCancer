#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p kellis
#SBATCH --array=0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ruiwenfu@mit.edu
#SBATCH --job-name=rna_VCF

module load bcftools
module load gatk
module load samtools

SAMPLE_LIST=($(<"pathRetrospctiveRNAseq.bam.txt"))
BAM=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}

SAMPLE=$(basename $BAM)
SAMPLE=${SAMPLE%%.*}

VCF=/net/bmc-lab5/data/kellis/group/011121-RNAseq/retrospective_hg38_variant/"$SAMPLE".vcf
SELECTED=/net/bmc-lab5/data/kellis/group/011121-RNAseq/retrospective_hg38_variant/selected/"$SAMPLE"RenamedSelected.vcf
REHEADED=/net/bmc-lab5/data/kellis/group/011121-RNAseq/retrospective_hg38_variant/selected/"$SAMPLE"RenamedSelectedReheaded.vcf
REORDERED=/net/bmc-lab5/data/kellis/group/011121-RNAseq/retrospective_hg38_variant/selected/"$SAMPLE"RenamedSelectedReheadedReordered.vcf

samtools index $BAM

bcftools mpileup -Ob -I -f /net/bmc-lab5/data/kellis/group/011121-RNAseq/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa $BAM \
| bcftools call -mv -Ov -o $VCF

bcftools annotate --rename-chrs chr_altNames.txt $VCF | grep -w '^#\|^#CHROM\|^chr[1-9]\|chr[1-2][0-9]\|chr[X]' > $SELECTED
bcftools reheader -f broad_hg38_v0_selected.fasta.fai $SELECTED -o $REHEADED

gatk UpdateVCFSequenceDictionary \
    -V $REHEADED \
    --source-dictionary broad_hg38_v0_selected.dict \
    --output $REORDERED \
    --replace=true

echo 'done'
