library(Seurat)
library(SeuratWrappers)
library(monocle3)

####################
# read in argument
####################

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("2 arguments must be supplied", call.=FALSE)
} else {
  objPATH <- args[1]
  figPATH <- args[2]
  root_cell <- args[3]
}

print(objPATH)
print(figPATH)
print(cdsPATH)
print(root_cell)


#############

obj<- readRDS(objPATH)
cds <- as.cell_data_set(obj)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

### plot tree
p <- plot_cells(cds, label_groups_by_cluster = TRUE, label_leaves = FALSE, label_branch_points = FALSE)
pdf(paste0(figPATH,'_monocle3_tree.pdf'))
print(p)
dev.off()

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin = root_cell){
  cell_ids <- which(colData(cds)[, "Ident"] == time_bin)

  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

### plot pseudotime
p <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
    label_branch_points = FALSE)
pdf(paste0(figPATH,'_monocle3_pseudotime.pdf'))
print(p)
dev.off()


obj <- AddMetaData(
  object = obj,
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "monocle3_pseudotime"
)


saveRDS(obj, objPATH)
