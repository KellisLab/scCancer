library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(patchwork)

############# load annotation

#### extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
### change to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"


############ load directory
path = '/net/bmc-lab5/data/kellis/group/Takeda_final/Takeda_single_cell/atac_cellranger_2.0.0/'
sample_id <- c('D21-194106','D19-9141','D19-11988','D19-11985',
               'D19-11984','D19-11983','D19-11981','D19-11980','D19-11973','D19-11972')

############# loop through all samples
for (sample in sample_id){
    print(paste("loading", sample))
    ############# load ATAC data

    counts <- Read10X_h5(filename = paste0(path, sample,'/outs/raw_peak_bc_matrix.h5'))
    metadata <- read.csv( file = paste0(path, sample,'/outs/singlecell.csv'),
                  header = TRUE, row.names = 1)
    chrom_assay <- CreateChromatinAssay(
              counts = counts,
              sep = c(":", "-"),
              genome = 'hg38',
              fragments = paste0(path, sample,'/outs/fragments.tsv.gz'),
              min.cells = 10,
              min.features = 200 )
    ATAC <- CreateSeuratObject(
              counts = chrom_assay,
              assay = "peaks",
              meta.data = metadata,
            project = sample)

    ### add the gene information to the object
    Annotation(ATAC) <- annotations
    print("adding annotation")
    ############# QC

    ### compute nucleosome signal score per cell
    ATAC <- NucleosomeSignal(object = ATAC)

    # compute TSS enrichment score per cell
    ATAC <- TSSEnrichment(object = ATAC, fast = FALSE)

    # add blacklist ratio and fraction of reads in peaks
    ATAC $pct_reads_in_peaks <- ATAC $peak_region_fragments / ATAC $passed_filters * 100
    ATAC $blacklist_ratio <- ATAC $blacklist_region_fragments / ATAC $peak_region_fragments

    #filter cells
    ATAC  <- subset(
      x = ATAC,
      subset = peak_region_fragments > 5000 &
        peak_region_fragments < 20000 &
        pct_reads_in_peaks > 25 &
        blacklist_ratio < 0.05 &
        nucleosome_signal < 4 &
        TSS.enrichment > 2
    )
    print("QC done")

    ############# create gene activity matrix
    gene.activities <- GeneActivity(ATAC)

    # add the gene activity matrix to the Seurat object as a new assay and normalize it
    ATAC[['gene.activities']] <- CreateAssayObject(counts = gene.activities)
    ATAC <- NormalizeData(
      object = ATAC,
      assay = 'gene.activities',
      normalization.method = 'LogNormalize',
      scale.factor = median(ATAC$nCount_gene.activities)
    )
    print("gene activity done")
    saveRDS(ATAC, file = paste("/home/ruiwenfu/data-lab5/scATAC/filtered/takedaICI_ATAC_filtered_" ,sample , "_102821.rds"))
}
