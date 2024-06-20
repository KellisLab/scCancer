
module load miniconda3

source activate /net/bmc-lab5/data/kellis/users/ruiwenfu/conda_related/envs/r-monocle

Rscript monocle3.R /home/ruiwenfu/data-lab5/scRNA/metastatic_all/takeda_39metastaticSamples_mt10_SCTmerge_CD4T_noTreg020622.rds \
/home/ruiwenfu/data-lab5/scRNA/metastatic_all/monocle3/takeda_39metastaticSamples_mt10_SCTmerge_CD4T_noTreg020622_ 0
echo 'monocle3 done'
