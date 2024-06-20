####################
# read in argument
####################

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("cell type name must be supplied", call.=FALSE)
} else {
  sample_name = do.call(paste, c(as.list(args), sep = " "))
}

print(paste("running part 2 on", sample_name))

####################
# load library
####################

library('dplyr')
library("tidyverse")
library(SCENIC)
library(pheatmap)
library('ggplot2')
library(AUCell)
####################
# load scenicOptions
####################

setwd("/net/bmc-lab5/data/kellis/users/ruiwenfu/SCENIC")
setwd(sample_name)

scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 16
scenicOptions@settings$seed <- 17

####################
# build and score GRN
####################

GRNBoost_output <- read.delim("int/1.1_network.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(GRNBoost_output) <- c("TF", "Target", "weight")
saveRDS(GRNBoost_output, file = "int/1.4_GENIE3_linkList.Rds")

exprMat_filtered_log <- t(read.table("int/1.1_exprMatrix_filtered_t.txt", header = T, stringsAsFactors = F, sep = "\t"))
cellInfo <- readRDS('int/cellInfo.Rds')
colnames(exprMat_filtered_log) <- rownames(cellInfo)


runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered_log)

print("runSCENIC1 2 3 done")

####################
# explore output
####################

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonAUC <- getAUC(regulonAUC)

n <- min(nrow(regulonAUC), 30)
selected <- order(apply(regulonAUC, 1, var), decreasing=TRUE)[1:n]
cellInfo.order <- cellInfo %>% arrange(match(state, c('PRE', 'ON')))

pheatmap <- pheatmap(regulonAUC[selected,rownames(cellInfo.order)],
     cluster_rows = T,
     cluster_cols = F,
     show_rownames = T,
     show_colnames = F,
     annotation = cellInfo.order,
     border_color = NA,
     fontsize = 8,
     fontsize_row = 8,
     height =8,
     width = 16,
     gaps_col = nrow(filter(cellInfo.order, state == 'PRE')),
    )

ggsave(filename=paste0(getwd(), '/output/postStep3_topRegulonActivity_heatmap.pdf'), plot=pheatmap,
       width = 18, height = 6)
print('heatmap done')
