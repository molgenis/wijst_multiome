#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_cpeaks_to_signac.R
# Function: perform normalization
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(Signac)
# this database is needed for annotations
library(EnsDb.Hsapiens.v86)
# for plots
library(ggplot2)


####################
# Functions        #
####################


read_cpeaks <- function(cpeaks_loc, cellranger_loc, lane, annotations) {
  # read the features
  features <- read.table(paste(cpeaks_loc, '/', lane, '/features.tsv.gz', sep = ''), sep = '\t', header = F)
  features[['peak']] <- paste(features[[1]], ':',  features[[2]], '-', features[[3]], sep = '')
  write.table(features, gzfile(paste(cpeaks_loc, '/', lane, '/named_features.tsv.gz', sep = '')), sep = '\t', row.names = F, col.names = F, quote = F)
  # read the actual counts
  peaks <- ReadMtx(mtx = paste(cpeaks_loc, '/', lane, '/atac_matrix.mtx.gz', sep = ''),
                   cells = paste(cpeaks_loc, '/', lane, '/barcodes.tsv.gz', sep = ''),
                   features = paste(cpeaks_loc, '/', lane, '/named_features.tsv.gz', sep = ''),
                   feature.column = 4)
  
  # create some metadata, for now, we'll first just store the lane here
  barcodes_short <- gsub('(-\\d+)', '', colnames(peaks))
  barcodes_lane <- paste(barcodes_short, rep(lane, times = length(barcodes_short)), sep = '_')
  # create metadata
  metadata <- data.frame(barcode=barcodes_short,
                         barcode_1=colnames(peaks),
                         barcode_lane=barcodes_lane)
  # save the lane the samples are from
  metadata$lane <- lane
  metadata$batch <- lane
  #rownames(metadata) <- metadata$barcode_lane
  # change to a barcode unique across lanes
  #colnames(peaks) <- metadata$barcode_lane
  # get the fragments
  fragments_loc <- paste(cellranger_loc, '/', lane, '/outs/atac_fragments.tsv.gz', sep = '')
  # create the chromatin assay
  chrom_assay <- CreateChromatinAssay(
    counts = peaks,
    sep = c(":", "-"),
    fragments = fragments_loc,
    min.cells = 10,
    min.features = 200
  )
  # create object
  seurat_object <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata,
    project = 'wijst_multiome'
  )
  # set the annotations to the object now
  Annotation(seurat_object) <- annotations
  # rename to lane and barcode
  seurat_object <- RenameCells(seurat_object, new.names = seurat_object@meta.data[['barcode_lane']])
  # return result
  return(seurat_object)
}


merge_cpeaks_objects <- function(cpeaks_loc, cellranger_loc, lanes, annotations, cpeak_prepend='mo_cpeaks_', cpeak_append='.rds', frag_append='outs/atac_fragments.tsv.gz') {
  # we will store all the metadata per lane
  lane_metadata <- list()
  # we also have to create a merged count data matrix, as the Seurat implementation is slow
  merged_atac_matrix <- NULL
  # we also have to store all the fragments
  lane_fragments <- list()
  # we'll read each lane
  for (lane in lanes) {
    # get the fragments location
    fragments_loc <- paste(cellranger_loc, '/', lane, '/', frag_append, sep = '')
    # get the seurat object
    cpeaks_object_loc <- paste(cpeaks_loc, '/', lane, '/', cpeak_prepend, lane, cpeak_append, sep = '')
    cpeaks_object <- readRDS(cpeaks_object_loc)
    # extract the metadata
    cpeaks_metadata <- cpeaks_object@meta.data
    # add that to the list
    lane_metadata[[lane]] <- cpeaks_metadata
    # get fragments info
    #fragments_barcodes <- lane_metadata[['barcode_1']]
    # give those vector names
    #names(fragments_barcodes) <- lane_metadata[['barcode_lane']]
    # create a fragments object
    #fragments <- CreateFragmentObject(
    #  path = fragments_loc,
    #  cells = fragments_barcodes, 
    #  validate.fragments = FALSE
    #)
    # put that into a list
    #lane_fragments[[lane]] <- fragments
    lane_fragments[[lane]] <- Fragments(cpeaks_object)[[1]]
    # extract the count matrix
    counts <- cpeaks_object@assays$peaks@counts
    print(paste('merging', lane))
    # check if we already have count data
    if (is.null(merged_atac_matrix)) {
      merged_atac_matrix <- counts
    }
    else {
      # get what is only present in the new count matrix
      only_mo_2_locs <- setdiff(rownames(counts), rownames(merged_atac_matrix))
      # get what is only present in the exisiting matrix
      only_mo_1_locs <- setdiff(rownames(merged_atac_matrix), rownames(counts))
      # add zero counts for the locations only present in each of the modalities
      merged_atac_matrix <- rbind(merged_atac_matrix, SparseEmptyMatrix(nrow = length(only_mo_2_locs), ncol = ncol(merged_atac_matrix), rownames = only_mo_2_locs, colnames = colnames(merged_atac_matrix)))
      counts <- rbind(counts, SparseEmptyMatrix(nrow = length(only_mo_1_locs), ncol = ncol(counts), rownames = only_mo_1_locs, colnames = colnames(counts)))
      # sort both of them
      merged_atac_matrix <- merged_atac_matrix[order(rownames(merged_atac_matrix)), ]
      counts <- counts[order(rownames(counts)), ]
      # now merge them
      merged_atac_matrix <- cbind(merged_atac_matrix, counts)
    }
  }
  # now merge all the metadata
  metadata_all <- do.call('rbind', lane_metadata)
  # create the chromatin assay
  chrom_assay <- CreateChromatinAssay(
    counts = merged_atac_matrix,
    sep = c(":", "-"),
    fragments = as.vector(unlist(lane_fragments)),
    min.cells = 10,
    min.features = 200
  )
  # create object
  seurat_object <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata_all,
    project = 'wijst_multiome'
  )
  # set the annotations to the object now
  Annotation(seurat_object) <- annotations
  return(seurat_object)
}


get_arc_metadata <- function(cellranger_loc, lanes, metadata_append='outs/per_barcode_metrics.csv') {
  # store per lane first
  metadata_per_lane <- list()
  # check each lane
  for (lane in lanes) {
    # paste the path together
    metadata_loc <- paste(cellranger_loc, '/', lane, '/', metadata_append, sep = '')
    # read the table
    metadata_table <- read.table(metadata_loc, sep = ',', header = T)
    # add lane as column
    metadata_table[['lane']] <- lane
    # create some metadata, for now, we'll first just store the lane here
    barcodes_short <- gsub('(-\\d+)', '', metadata_table[['barcode']])
    barcodes_lane <- paste(barcodes_short, rep(lane, times = length(barcodes_short)), sep = '_')
    # add short barcode and lane+barcode
    metadata_table[['barcode_1']] <- metadata_table[['barcode']]
    metadata_table[['barcode_lane']] <- barcodes_lane
    metadata_table[['barcode']] <- barcodes_short
    # put in list
    metadata_per_lane[[lane]] <- metadata_table
  }
  # merge all together
  metadata_all <- do.call('rbind', metadata_per_lane)
  rownames(metadata_all) <- metadata_all[['barcode_lane']]
  return(metadata_all)
}

ref10xmo_predictions_to_lower_res_mapping <- function() {
  high_to_low <- list()
  high_to_low[['CD4 Naive']] <- 'CD4T'
  high_to_low[['CD4 TCM']] <- 'CD4T'
  high_to_low[['CD8 Naive']] <- 'CD8T'
  high_to_low[['CD16 Mono']] <- 'monocyte'
  high_to_low[['NK']] <- 'NK'
  high_to_low[['Treg']] <- 'T_other'
  high_to_low[['CD14 Mono']] <- 'monocyte'
  high_to_low[['cDC']] <- 'DC'
  high_to_low[['CD8 TEM_1']] <- 'CD8T'
  high_to_low[['Intermediate B']] <- 'B'
  high_to_low[['Naive B']] <- 'B'
  high_to_low[['Plasma']] <- 'plasmablast'
  high_to_low[['CD4 TEM']] <- 'CD4T'
  high_to_low[['MAIT']] <- 'T_other'
  high_to_low[['Memory B']] <- 'B'
  high_to_low[['gdT']] <- 'T_other'
  high_to_low[['pDC']] <- 'DC'
  high_to_low[['CD8 TEM_2']] <- 'CD8T'
  high_to_low[['HSPC']] <- 'hemapoietic_stem'
  return(high_to_low)
}


get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UT"]] <- "lightgrey"
  color_coding[["24hCA"]] <- "forestgreen"
  color_coding[["24hCa"]] <- "forestgreen"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  color_coding[["CD4+ T"]] <- "#153057"
  color_coding[["CD8+ T"]] <- "#009DDB"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  color_coding[["T_other"]] <- "#FF63B6"
  color_coding[["hemapoietic_stem"]] <- "#8B8000"
  color_coding[["hemapoietic stem"]] <- "#8B8000"
  return(color_coding)
}


label_dict <- function(){
  label_dict <- list()
  # condition combinations
  label_dict[["UT"]] <- "UT"
  label_dict[["24hCa"]] <- "24hCA"
  label_dict[["24hCA"]] <- "24hCA"
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
  label_dict[["CD4T"]] <- "CD4+ T"
  label_dict[["CD8T"]] <- "CD8+ T"
  label_dict[["monocyte"]] <- "monocyte"
  label_dict[["NK"]] <- "NK"
  label_dict[["B"]] <- "B"
  label_dict[["DC"]] <- "DC"
  label_dict[["HSPC"]] <- "HSPC"
  label_dict[["plasmablast"]] <- "plasmablast"
  label_dict[["platelet"]] <- "platelet"
  label_dict[["T_other"]] <- "other T"
  label_dict[["hemapoietic_stem"]] <- "hemapoietic stem"
  # minor cell types
  label_dict[["CD4_TCM"]] <- "CD4 TCM"
  label_dict[["Treg"]] <- "T regulatory"
  label_dict[["CD4_Naive"]] <- "CD4 naive"
  label_dict[["CD4_CTL"]] <- "CD4 CTL"
  label_dict[["CD8_TEM"]] <- "CD8 TEM"
  label_dict[["cMono"]] <- "cMono"
  label_dict[["CD8_TCM"]] <- "CD8 TCM"
  label_dict[["ncMono"]] <- "ncMono"
  label_dict[["cDC2"]] <- "cDC2"
  label_dict[["B_intermediate"]] <- "B intermediate"
  label_dict[["NKdim"]] <- "NK dim"
  label_dict[["pDC"]] <- "pDC"
  label_dict[["ASDC"]] <- "ASDC"
  label_dict[["CD8_Naive"]] <- "CD8 naive"
  label_dict[["MAIT"]] <- "MAIT"
  label_dict[["CD8_Proliferating"]] <- "CD8 proliferating"
  label_dict[["CD4_TEM"]] <- "CD4 TEM"
  label_dict[["B_memory"]] <- "B memory"
  label_dict[["NKbright"]] <- "NK bright"
  label_dict[["B_naive"]] <- "B naive"
  label_dict[["gdT"]] <- "gamma delta T"
  label_dict[["CD4_Proliferating"]] <- "CD4 proliferating"
  label_dict[["NK_Proliferating"]] <- "NK proliferating"
  label_dict[["cDC1"]] <- "cDC1"
  label_dict[["ILC"]] <- "ILC"
  label_dict[["dnT"]] <- "double negative T"
  return(label_dict)
}



####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')
set.seed(7777)

# get annotations from ensemble database
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# we use UCSC gencode
seqlevelsStyle(annotations) <- "UCSC"
# and this was aligned on build 38
genome(annotations) <- "hg38"


####################
# Main Code        #
####################

# these are the lanes
lanes <- c('230105_lane1', '230105_lane2', '230105_lane3', '230105_lane4',
           '230105_lane5', '230105_lane6', '230105_lane7', '230105_lane8'
           # '230112_lane1', '230112_lane2', '230112_lane3', '230112_lane4',
           # '230112_lane5', '230112_lane6', '230112_lane7', '230112_lane8',
           # '230120_lane1', '230120_lane2', '230120_lane3', '230120_lane4',
           # '230120_lane5', '230120_lane6', '230120_lane7', '230120_lane8',
           # '230127_lane1', '230127_lane2', '230127_lane3', '230127_lane4',
           # '230127_lane5', '230127_lane6', '230127_lane7', '230127_lane8',
           # '230202_lane1', '230202_lane2', '230202_lane3', '230202_lane4',
           # '230202_lane5', '230202_lane6', '230202_lane7', '230202_lane8',
           # '230209_lane1', '230209_lane2', '230209_lane3', '230209_lane4',
           # '230209_lane5', '230209_lane6', '230209_lane7', '230209_lane8',
           # '230216_lane1', '230216_lane2', '230216_lane3', '230216_lane4',
           # '230216_lane5', '230216_lane6', '230216_lane7', '230216_lane8',
           # '230223_lane1', '230223_lane2', '230223_lane3', '230223_lane4',
           # '230223_lane5', '230223_lane6', '230223_lane7', '230223_lane8',
           # '230302_lane1', '230302_lane2', '230302_lane3', '230302_lane4',
           # '230302_lane5', '230302_lane6', '230302_lane7', '230302_lane8',
           # '230316_lane1', '230316_lane2', '230316_lane3', '230316_lane4',
           # '230316_lane5', '230316_lane6', '230316_lane7', '230316_lane8'
)

# location of the cpeaks objects
cpeaks_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/cpeaks_peak_calling/output/default/'
# get cellranger loc
cellranger_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/processed/joint/alignment/b38/'
# location of the gene annotations
gtf_loc <- '/groups/umcg-franke-scrna/tmp03/external_datasets/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'

# check each lane
for (lane in lanes) {
  # read the object
  object_lane <- read_cpeaks(cpeaks_loc, cellranger_loc, lane = lane, annotations = annotations)
  # save the result
  cpeaks_object_loc <- paste(cpeaks_loc, '/', lane, '/', 'mo_cpeaks_', lane, '.rds', sep = '')
  saveRDS(object_lane, cpeaks_object_loc)
  rm(object_lane)
}

# add all together
mo_all <- merge_cpeaks_objects(cpeaks_loc, cellranger_loc, lanes, annotations)
# reremoved cells not having enough counts
mo_all <- mo_all[, mo_all$nFeature_peaks > 200]
# compute nucleosome signal score per cell
mo_all <- NucleosomeSignal(object = mo_all)
# compute TSS enrichment score per cell
mo_all <- TSSEnrichment(object = mo_all, fast = FALSE)

# get the metdatadata
arc_metadata <- get_arc_metadata(cellranger_loc, lanes)
mo_all <- AddMetaData(mo_all, arc_metadata['atac_peak_region_fragments'], 'peak_region_fragments')
mo_all <- AddMetaData(mo_all, arc_metadata['atac_fragments'], 'atac_fragments')

# now the blacklist region
mo_all$blacklist_fraction <- FractionCountsInRegion(
  object = mo_all, 
  assay = 'peaks',
  regions = blacklist_hg38
)

# add blacklist ratio and fraction of reads in peaks
#mo_all$pct_reads_in_peaks <- mo_all$peak_region_fragments / mo_all$passed_filters * 100
mo_all$pct_reads_in_peaks <- mo_all$peak_region_fragments / mo_all$atac_fragments * 100
mo_all$blacklist_ratio <- mo_all$blacklist_fraction / mo_all$peak_region_fragments
mo_all$nucleosome_group <- ifelse(mo_all$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
# save intermediate result
saveRDS(mo_all, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/cpeaks_peak_calling/output/default/mo_cpeaks_unfiltered_230105_lane1-8.rds')
# remove outliers
ncol(mo_all)
# [1] 96919
mo_all <- subset(x = mo_all, subset = nCount_peaks > 3000)
ncol(mo_all)
# [1] 88481
mo_all <- subset(x = mo_all, subset = nCount_peaks < 30000)
ncol(mo_all)
# [1] 81367
mo_all <- subset(x = mo_all, subset = pct_reads_in_peaks > 15)
ncol(mo_all)
# [1] 81367
mo_all <- subset(x = mo_all, subset = blacklist_ratio < 0.05)
ncol(mo_all)
# [1] 81367
mo_all <- subset(x = mo_all, subset = nucleosome_signal < 4)
ncol(mo_all)
# [1] 81367
mo_all <- subset(x = mo_all, subset = TSS.enrichment > 3)
ncol(mo_all)
# [1] 81265

# read the RNA level metadata
rna_metadata <- read.table('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz', header = T, sep = '\t')
# add a final cell type, by taking the prediction and adding the imputed where it was empty
rna_metadata[['cell_type_final']] <- rna_metadata[['predicted.mo_10x_cell_type']]
rna_metadata[is.na(rna_metadata[['cell_type_final']]), 'cell_type_final'] <- rna_metadata[is.na(rna_metadata[['cell_type_final']]), 'celltype_imputed']
# set rownames to be barcode and lane
rownames(rna_metadata) <- rna_metadata[['barcode_lane']]
# add more metadata
mo_all <- AddMetaData(mo_all, rna_metadata['cell_type_final'], 'cell_type_final')
mo_all <- AddMetaData(mo_all, rna_metadata['sample_final'], 'sample_final')
mo_all <- AddMetaData(mo_all, rna_metadata['final_condition'], 'final_condition')
mo_all <- AddMetaData(mo_all, rna_metadata['soup_status'], 'soup_status')
# remove doublets
mo_all <- mo_all[, !is.na(mo_all@meta.data[['soup_status']]) & mo_all@meta.data[['soup_status']] == 'singlet']
ncol(mo_all)
# [1] 67596

# do normalization
mo_all <- RunTFIDF(mo_all)
mo_all <- FindTopFeatures(mo_all, min.cutoff = 'q0')
mo_all <- RunSVD(mo_all)
# clustering and UMAP
mo_all <- RunUMAP(object = mo_all, reduction = 'lsi', dims = 1:30)
mo_all <- FindNeighbors(object = mo_all, reduction = 'lsi', dims = 1:30)
mo_all <- FindClusters(object = mo_all, verbose = FALSE, algorithm = 3)
DimPlot(object = mo_all, label = TRUE) + NoLegend()

# add lower cell type classification
mo_all@meta.data[['cell_type_final_lowerres']] <- as.vector(unlist(ref10xmo_predictions_to_lower_res_mapping()[mo_all@meta.data[['cell_type_final']]]))

# add the gene activity matrix to the Seurat object as a new assay and normalize it
gene_activities <- GeneActivity(mo_all)
mo_all[['activity']] <- CreateAssayObject(counts = gene_activities)
mo_all <- NormalizeData(
  object = mo_all,
  assay = 'activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(mo_all$nCount_activity)
)
saveRDS(mo_all, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/cpeaks_peak_calling/output/default/mo_cpeaks_230105_lane1-8.rds')
