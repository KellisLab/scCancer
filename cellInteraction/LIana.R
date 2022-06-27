library('Seurat')
library('dplyr')
library('tidyr')
library('ggplot2')
library('Matrix')

## merge objects
obj1  <- readRDS("/net/bmc-lab5/data/kellis/users/ruiwenfu/scRNA/metastatic_all/takeda_39metastaticSamples_mt10_SCT_RPCAintegrated_cDC021022.rds")
#obj2  <- readRDS("/net/bmc-lab5/data/kellis/users/ruiwenfu/scRNA/metastatic_all/takeda_39metastaticSamples_mt10_SCTmerge_CD4T011422.rds")
obj2 <- readRDS( "/net/bmc-lab5/data/kellis/users/ruiwenfu/scRNA/metastatic_all/takeda_39metastaticSamples_mt10_SCTmerge_CD8Tcleaned020922.rds")


DefaultAssay(object = obj1) <- "RNA"
obj1 <- DietSeurat(obj1, assays = 'RNA')

DefaultAssay(object = obj2) <- "RNA"
obj2 <- DietSeurat(obj2, assays = 'RNA')
meta <- obj2@meta.data %>%
    mutate(Ident = paste0('CD8T_', obj2[[]]$Ident))
obj2@meta.data  <- meta

obj <- merge(obj1, obj2)
Idents(object = obj) <- "Ident"

obj.list <- SplitObject(obj, split.by = "uid")

## run LIana

require(liana)
require(data.table)

cpdb_results <- vector(mode = "list", length = 0)
for (id in names(obj.list)){
    if(ncol(obj.list[[id]]) < 50 | length(unique(obj.list[[id]]$Ident)) <2) next
    res <- liana_wrap(obj.list[[id]],
                        method = 'cellphonedb',
                        resource = c('CellPhoneDB'),
                        permutation.params = list(nperms=10000,
                                                  parallelize=TRUE,
                                                  workers=16))
    res$sample <- id
    cpdb_results[[id]] <- res
}
cpdb <- rbindlist(cpdb_results)
write.csv(cpdb, '/net/bmc-lab5/data/kellis/users/ruiwenfu/scRNA/metastatic_all/cpdb/CD8T_cDC_cpdbResults042622.csv')
