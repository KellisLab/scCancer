####################
# read in argument
####################

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("cell type name must be supplied", call.=FALSE)
} else {
  sample_name = do.call(paste, c(as.list(args), sep = "_"))
}

print(paste("running part 1 on", sample_name))

####################
# load library
####################

library('Seurat')
library('dplyr')
library("tidyverse")
library(SCENIC)

####################
# extract data from seurat object
####################

subtyped <-readRDS("/net/bmc-lab5/data/kellis/users/ruiwenfu/scRNA/metastatic_all/takeda_39metastaticSamples_mt10_SCT_RPCAintegrated_cDC021022.rds")

sample <- subset(subtyped, subset = (Ident == gsub('_', ' ', sample_name)))


exprMat <- GetAssayData(object = sample, assay = "RNA", slot = "counts")
cellInfo <- sample@meta.data[c("nCount_RNA", "nFeature_RNA", "percent_mt", "tissue","uid")]

print(c('exprMat dim is', dim(exprMat)))

####################
# iniitialize scenic
####################

setwd("/net/bmc-lab5/data/kellis/users/ruiwenfu/SCENIC")
dir.create(sample_name)
setwd(sample_name)

org <- "hgnc"
myDatasetTitle <- paste("SCENIC PRE/ON ICI", sample_name) # choose a name for your analysis
dbDir <- "/net/bmc-lab5/data/kellis/users/ruiwenfu/cisTarget_databases" # RcisTarget databases location
hg38_dbs <- c('500bp' = 'hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather', '10kb' = 'hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')

scenicOptions <- initializeScenic(org=org, datasetTitle=myDatasetTitle, dbs = hg38_dbs,
                                dbDir = dbDir, nCores=16)
scenicOptions@settings$db_mcVersion <- 'v9'


dir.create("int")

saveRDS(cellInfo, file="int/cellInfo.Rds")
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

print('iniitialize scenic done')

####################
# filter genes
####################

exprMat <- as.matrix(exprMat)

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.1*ncol(exprMat),
                           minSamples=ncol(exprMat)*.1)

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

####################
# co-exression network and export for arboreto
####################

runCorrelation(exprMat_filtered, scenicOptions)

exprMat_filtered_log <- log2(exprMat_filtered+1)
exportsForArboreto(exprMat_filtered_log, scenicOptions)

print('export for arboreto done')
