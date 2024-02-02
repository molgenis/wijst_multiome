# reinstall to make sure we have correct versions
install.packages('Seurat')
remotes::install_github("stuart-lab/signac", ref = 'develop')
# I need these for the annotations apparently
BiocManager::install("biovizBase")
BiocManager::install("EnsDb.Hsapiens.v86")
# I need these for the annotations apparently
BiocManager::install("biovizBase")
BiocManager::install("EnsDb.Hsapiens.v86")
install.packages('irlba')

# load libraries
library(Seurat)
library(Signac)
# and need to load this
library(EnsDb.Hsapiens.v86)

# get annotations from ensemble database
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# we use UCSC gencode
seqlevelsStyle(annotations) <- "UCSC"
# and this was aligned on build 38
genome(annotations) <- "hg38"

# load the lane I rsynced
lane <- readRDS('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_230302_lane6multimodal_azi_mapped.rds')
# set to the ATAC assay
DefaultAssay(lane) <- 'peaks'

# set the annotations to the object now
Annotation(lane) <- annotations
# hack the path to be vaxtron compatible
lane$peaks@fragments[[1]]@path <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/processed/joint/alignment/b38/230302_lane6/outs/atac_fragments.tsv.gz'
# compute nucleosome signal score per cell
lane <- NucleosomeSignal(object = lane)
# compute TSS enrichment score per cell
lane <- TSSEnrichment(object = lane, fast = FALSE)
# show TSS
lane$high.tss <- ifelse(lane$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(lane, group.by = 'high.tss') + NoLegend()
# RunTFIDF will fail unless you reinstall this
install.packages('irlba')
# do normalization, already done
#lane <- RunTFIDF(lane)
#lane <- FindTopFeatures(lane, min.cutoff = 'q0')
#lane <- RunSVD(lane)
# calculate gene activity
gene.activities <- GeneActivity(lane)
# set the cell type as the identity
Idents(lane) <- 'predicted.mo_10x_cell_type'
# calculate the peaks
da_peaks <- FindAllMarkers(
  object = lane,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
# we only want the ones that are more open for a specific cell type
da_peaks_open <- da_peaks[da_peaks[['p_val_adj']] < 0.05 & da_peaks[['avg_log2FC']] > 0, ]
# now get the closest genes
da_closest_open_features <- ClosestFeature(lane, regions = "NEED TO KNOW WHAT THE COLUMN IS CALLED")

# better stack trace
options(error = function() {
  calls <- sys.calls()
  if (length(calls) >= 2L) {
    sink(stderr())
    on.exit(sink(NULL))
    cat("Backtrace:\n")
    calls <- rev(calls[-length(calls)])
    for (i in seq_along(calls)) {
      cat(i, ": ", deparse(calls[[i]], nlines = 1L), "\n", sep = "")
    }
  }
  if (!interactive()) {
    q(status = 1)
  }
})
