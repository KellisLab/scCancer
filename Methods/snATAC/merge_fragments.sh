#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p kellis
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ruiwenfu@mit.edu
#SBATCH --job-name=merge_frag

module load samtools/1.10
cd /net/bmc-lab5/data/kellis/users/ruiwenfu/scATAC/metastatic_all

# while read F  ; do
#         echo $F
#         ID=$F'_'
#         gzip -dc '/net/bmc-lab5/data/kellis/group/Takeda_final/Takeda_single_cell/atac_cellranger_2.0.0/'$F'/outs/fragments.tsv.gz' | \
#         awk -v x="$ID" 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,x$4,$5}' - > 'fragments/'$F'_fragments.tsv'
#
# done <ATAC_ids_striped.txt
cd fragments

for file in *_fragments.tsv
do
  sed '/^#/d' $file > 'noHeader_'$file
done


# sort all files together
sort -k1,1V -k2,2n noHeader_* > noHeader_fragments.tsv

# block gzip compress the merged file
bgzip -@ 8 noHeader_fragments.tsv # -@ 4 uses 4 threads

# index the bgzipped file
tabix -p bed noHeader_fragments.tsv.gz
