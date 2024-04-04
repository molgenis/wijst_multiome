#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_cpeaks_to_celltypes.R
# Function: 
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
library(cowplot)
# for upset plot of peak sharing
library(UpSetR)


####################
# Functions        #
####################

process_signac_object <- function(signac_object, arc_metadata, npeaks_cutoff=200, fragments_loc=NULL, fragment_prepend='', fragment_append='_rounded_fragments.tsv.gz') {
  # update the fragment loc if requested
  if (!is.null(fragments_loc)) {
    # check each lane
    for (lane in names(signac_object@assays$peaks@fragments)) {
      # and update the location
      signac_object@assays$peaks@fragments[[lane]]@path <- paste(fragments_loc, '/', fragment_prepend, lane, fragment_append, sep = '')
    }
  }
  # add relevant metadata
  signac_object <- AddMetaData(signac_object, arc_metadata['atac_peak_region_fragments'], 'peak_region_fragments')
  signac_object <- AddMetaData(signac_object, arc_metadata['atac_fragments'], 'atac_fragments')
  print(paste('cells before nfeatures in npeaks', ncol(signac_object)))
  # reremoved cells not having enough counts
  signac_object <- signac_object[, signac_object$nFeature_peaks > npeaks_cutoff]
  print(paste('cells after nfeatures in npeaks', ncol(signac_object)))
  # compute nucleosome signal score per cell
  signac_object <- NucleosomeSignal(object = signac_object)
  # compute TSS enrichment score per cell
  signac_object <- TSSEnrichment(object = signac_object, fast = FALSE)
  # now the blacklist region
  signac_object$blacklist_fraction <- FractionCountsInRegion(
    object = signac_object, 
    assay = 'peaks',
    regions = blacklist_hg38
  )
  signac_object$pct_reads_in_peaks <- signac_object$peak_region_fragments / signac_object$atac_fragments * 100
  signac_object$blacklist_ratio <- signac_object$blacklist_fraction / signac_object$peak_region_fragments
  signac_object$nucleosome_group <- ifelse(signac_object$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  return(signac_object)
}

qc_signac_object <- function(signac_object, min_nCount_peaks=3000, max_nCount_peaks=30000, min_pct_reads_in_peaks=50, max_blacklist_ratio=0.05, max_nucleosome_signal=4, min_TSS.enrichment=3) {
  # remove outliers
  print(paste('ncells pre QC:', as.character(ncol(signac_object))))
  signac_object <- subset(x = signac_object, subset = nCount_peaks > min_nCount_peaks)
  print(paste('ncells post-min_nCount_peaks QC:', as.character(ncol(signac_object))))
  signac_object <- subset(x = signac_object, subset = nCount_peaks < max_nCount_peaks)
  print(paste('ncells post-max_nCount_peaks:', as.character(ncol(signac_object))))
  signac_object <- subset(x = signac_object, subset = pct_reads_in_peaks >= min_pct_reads_in_peaks)
  print(paste('ncells post-min_pct_reads_in_peaks:', as.character(ncol(signac_object))))
  signac_object <- subset(x = signac_object, subset = blacklist_ratio < max_blacklist_ratio)
  print(paste('ncells post-max_blacklist_ratio:', as.character(ncol(signac_object))))
  signac_object <- subset(x = signac_object, subset = nucleosome_signal < max_nucleosome_signal)
  print(paste('ncells post-max_nucleosome_signal:', as.character(ncol(signac_object))))
  signac_object <- subset(x = signac_object, subset = TSS.enrichment > min_TSS.enrichment)
  print(paste('ncells post-min_TSS.enrichment:', as.character(ncol(signac_object))))
  return(signac_object)
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
  color_coding[["T-other"]] <- "#FF63B6"
  color_coding[["hemapoietic_stem"]] <- "#8B8000"
  color_coding[["hemapoietic stem"]] <- "#8B8000"
  color_coding[["hemapoietic-stem"]] <- "#8B8000"
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
  label_dict[["T-other"]] <- "other T"
  label_dict[["hemapoietic_stem"]] <- "hemapoietic stem"
  label_dict[["hemapoietic-stem"]] <- "hemapoietic stem"
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


peaks_to_bed <- function(exp_per_group, output_prepend='mo_peaks_', output_append='.bed', min_log10_peak=0) {
  # check each of the columns
  for (identity in colnames(exp_per_group)) {
    # subset to that identity
    peaks_ident <- exp_per_group[, c(identity), drop = F]
    # get the locations
    positions <- data.frame(do.call('rbind', strsplit(rownames(peaks_ident), '-')))
    # set the colnames properly to bed format
    colnames(positions) <- c('#chrom', 'start', 'end')
    # add the rownames themselves as the name
    positions[['name']] <- rownames(peaks_ident)
    # and the counts as score
    positions[['score']] <- log10(peaks_ident[, c(identity)])
    # remove empty entries
    positions <- positions[positions[['score']] >= min_log10_peak, ]
    # make the start and stop numeric
    positions[['start']] <- as.numeric(positions[['start']])
    positions[['end']] <- as.numeric(positions[['end']])
    # and sort
    positions <- positions[order(positions[['#chrom']], positions[['start']], positions[['end']]), ]
    # create the output location
    output_loc <- paste(output_prepend, identity, output_append, sep = '')
    # write the result
    write.table(positions, output_loc, row.names = F, col.names = T, quote = F, sep = '\t')
  }
}

get_npeaks_per_threshold <- function(exp_per_group, log10_thresholds=c(0,1,2,3,4)) {
  # we'll do each identity, and store the result in per ident in a list
  numbers_per_identity <- list()
  # check each of the columns
  for (identity in colnames(exp_per_group)) {
    # subset to that identity
    peaks_ident <- exp_per_group[, c(identity), drop = F]
    # get the locations
    positions <- data.frame(do.call('rbind', strsplit(rownames(peaks_ident), '-')))
    # set the colnames properly to bed format
    colnames(positions) <- c('#chrom', 'start', 'end')
    # add the rownames themselves as the name
    positions[['name']] <- rownames(peaks_ident)
    # and the counts as score
    positions[['score']] <- log10(peaks_ident[, c(identity)])
    # init dataframe
    numbers = data.frame(identity = rep(identity, times = length(log10_thresholds)), 
                         log10_thresh=log10_thresholds, 
                         npeaks=rep(NA, times = length(log10_thresholds)))
    # check each threshold and fill in that dataframe
    for (thresh_i in 1 : length(log10_thresholds)) {
      # get that threshold
      threshold <- log10_thresholds[thresh_i]
      # check now many peaks are left
      npeaks <- nrow(positions[positions[['score']] >= threshold, ])
      # put that into the dataframe
      numbers[thresh_i, 'npeaks'] <- npeaks
    }
    # put the result in the list we created
    numbers_per_identity[[identity]] <- numbers
  }
  # now merge all of the idents together in one dataframe
  numbers_all <- do.call('rbind', numbers_per_identity)
  return(numbers_all)
}


get_peaks_per_identity <- function(exp_per_group, threshold=0) {
  # we'll do each identity, and store the result in per ident in a list
  peaks_per_identity <- list()
  # check each of the columns
  for (identity in colnames(exp_per_group)) {
    # subset to that identity
    peaks_ident <- exp_per_group[, c(identity), drop = F]
    # get the locations
    positions <- data.frame(do.call('rbind', strsplit(rownames(peaks_ident), '-')))
    # set the colnames properly to bed format
    colnames(positions) <- c('#chrom', 'start', 'end')
    # add the rownames themselves as the name
    positions[['name']] <- rownames(peaks_ident)
    # and the counts as score
    positions[['score']] <- log10(peaks_ident[, c(identity)])
    # check now many peaks are left
    peaks <- positions[positions[['score']] >= threshold, 'name']
    # put the result in the list we created
    peaks_per_identity[[identity]] <- peaks
  }
  return(peaks_per_identity)
}


plot_peak_sharing_per_celltype <- function(peaks_per_ct, use_label_dict=T, use_color_dict=T){
  # rename cell types if requested
  if(use_label_dict){
    names(peaks_per_ct) <- label_dict()[names(peaks_per_ct)]
  }
  # we'll need this list for colouring
  queries <- list()
  sets.bar.color <- 'black'
  # do colouring
  if(use_color_dict){
    # create df to store the number of each set, so we know how to order
    nrs_df <- NULL
    # add the colors for the cell types
    i <- 1
    for(cell_type in (names(peaks_per_ct))){
      # get all peaks except this cell type
      peaks_except_c <- do.call('c', peaks_per_ct[setdiff(names(peaks_per_ct), cell_type)])
      # then this cell type
      peaks_celltype<- peaks_per_ct[[cell_type]]
      # check if this singlet exists
      if (length(setdiff(peaks_celltype, peaks_except_c)) > 0) {
        # add for the singles in the intersection sizes
        ct_list <- list(
          query = intersects,
          params = list(cell_type),
          color = get_color_coding_dict()[[cell_type]],
          active = T)
        queries[[i]] <- ct_list
        i <- i + 1
      }
      else {
        print('blegh')
      }
      
      # add for the DF to order the set sizes
      numbers_row <- data.frame(ct=c(cell_type), nr=c(length(peaks_per_ct[[cell_type]])), stringsAsFactors = F)
      if(is.null(nrs_df)){
        nrs_df <- numbers_row
      }
      else{
        nrs_df <- rbind(nrs_df, numbers_row)
      }
    }
    # get the order of the sets
    ordered_cts <- nrs_df[order(nrs_df$nr, decreasing = T), 'ct']
    # add the colors for the sets
    sets.bar.color <- unlist(get_color_coding_dict()[ordered_cts])
  }
  else {
    queries = NULL
  }
  #upset(fromList(peaks_per_ct), order.by = 'freq', nsets = length(peaks_per_ct), queries = queries, sets.bar.color=sets.bar.color, nintersects = length(peaks_per_ct)*length(peaks_per_ct))
  upset(fromList(peaks_per_ct), order.by = 'freq', nsets = length(peaks_per_ct), sets.bar.color=sets.bar.color, nintersects = length(peaks_per_ct)*length(peaks_per_ct))
  #return(peaks_per_ct)
}


merge_signac_objects <- function(signac_object_vector) {
  signac_metadata_list <- list()
  signac_fragments_lists <- list()
  merged_atac_matrix <- NULL
  # check each object
  for (i in 1:length(signac_object_vector)) {
    # extract the metadata
    signac_metadata_list[[i]] <- signac_object_vector[[i]]@meta.data
    signac_metadata_list[[i]] <- Fragments(signac_object_vector[[i]])
    # extract the count data
    counts <- signac_object_vector[[i]]@assays$peaks@counts
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
  signac_metadata <- do.call('rbind', signac_metadata_list)
  # now merge all the fragments
  signac_fragments <- do.call('c', signac_fragments_lists)
  # create the chromatin assay
  chrom_assay <- CreateChromatinAssay(
    counts = merged_atac_matrix,
    sep = c(":", "-"),
    fragments = signac_fragments,
    min.cells = 10,
    min.features = 200
  )
  # create object
  seurat_object <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = signac_metadata,
    project = 'wijst_multiome'
  )
  # set the annotations to the object now
  Annotation(seurat_object) <- annotations
  return(seurat_object)
}


merge_signac_per_celltypes <- function(signac_object_vector, cell_type_column='cell_type_final') {
  # save an object per cell type
  merged_objects_celltypes <- list()
  # extract the cell types
  cell_types_to_do <- unique(signac_object_vector[[1]]@meta.data[[cell_type_column]])
  # remove NA
  cell_types_to_do <- cell_types_to_do[!is.na(cell_types_to_do)]
  # check each cell type
  for (cell_type in cell_types_to_do) {
    # check each object first
    objects_cell_type <- list()
    # do each object
    for (i in 1:length(signac_object_vector)) {
      # subset to cell type
      signac_object_celltype <- signac_object_vector[[i]][, signac_object_vector[[i]]@meta.data[[cell_type_column]] == cell_type]
      # add to the vector
      objects_cell_type[[i]] <- signac_object_celltype
    }
    # now merge for all of the celltype
    merged_object_celltype <- merge_signac_objects(objects_cell_type)
    # add result to list
    merged_objects_celltypes[[cell_type]] <- merged_object_celltype
  }
  return(merged_objects_celltypes)
}

make_celltypes_safe <- function(cell_types){
  # get a safe file name
  cell_type_safes <- gsub(' |/', '_', cell_types)
  cell_type_safes <- gsub('-', '_negative', cell_type_safes)
  cell_type_safes <- gsub('\\+', '_positive', cell_type_safes)
  cell_type_safes <- gsub('\\)', '', cell_type_safes)
  cell_type_safes <- gsub('\\(', '', cell_type_safes)
  return(cell_type_safes)
}


signac_dimreduc_and_cluster <- function(signac_object) {
  # do normalization
  signac_object <- RunTFIDF(signac_object)
  signac_object <- FindTopFeatures(signac_object, min.cutoff = 'q0')
  signac_object <- RunSVD(signac_object)
  # clustering and UMAP
  signac_object <- RunUMAP(object = signac_object, reduction = 'lsi', dims = 1:30)
  signac_object <- FindNeighbors(object = signac_object, reduction = 'lsi', dims = 1:30)
  signac_object <- FindClusters(object = signac_object, verbose = FALSE, algorithm = 3)
  DimPlot(object = signac_object, label = TRUE) + NoLegend()
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  gene_activities <- GeneActivity(signac_object)
  signac_object[['activity']] <- CreateAssayObject(counts = gene_activities)
  signac_object <- NormalizeData(
    object = signac_object,
    assay = 'activity',
    normalization.method = 'LogNormalize',
    scale.factor = median(signac_object$nCount_activity)
  )
  return(signac_object)
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

# location of the metadata
arc_metadata_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/arc_metadata.tsv.gz'

# location of the fragments
fragments_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/rounded_fragments/'

# these are the objects
mo_object_1_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_unfiltered1_16.rds'
mo_object_1_filtered_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered1_16.rds'
mo_object_1_clustered_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_clustered1_16.rds'
mo_object_2_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_unfiltered17_32.rds'
mo_object_2_filtered_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered17_32.rds'
mo_object_2_clustered_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_clustered17_32.rds'
mo_object_3_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_unfiltered33_48.rds'
mo_object_3_filtered_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered33_48.rds'
mo_object_3_clustered_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_clustered33_48.rds'
mo_object_4_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_unfiltered49_64.rds'
mo_object_4_filtered_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered49_64.rds'
mo_object_4_clustered_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_clustered49_64.rds'


# load metadata
arc_metadata <- read.table(arc_metadata_loc, header = T, sep = '\t')
rownames(arc_metadata) <- arc_metadata[['barcode_lane']]

# read the RNA level metadata
rna_metadata <- read.table('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz', header = T, sep = '\t')
# add a final cell type, by taking the prediction and adding the imputed where it was empty
rna_metadata[['cell_type_final']] <- rna_metadata[['predicted.mo_10x_cell_type']]
rna_metadata[is.na(rna_metadata[['cell_type_final']]), 'cell_type_final'] <- rna_metadata[is.na(rna_metadata[['cell_type_final']]), 'celltype_imputed']
# set rownames to be barcode and lane
rownames(rna_metadata) <- rna_metadata[['barcode_lane']]

# read the first object
mo_object_2 <- readRDS(mo_object_2_loc)
# process
mo_object_2 <- process_signac_object(mo_object_2, arc_metadata, fragments_loc = fragments_loc)
mo_object_2 <- qc_signac_object(mo_object_2)
#[1] "ncells pre QC: 151251"
#[1] "ncells post-min_nCount_peaks QC: 132108"
#[1] "ncells post-max_nCount_peaks: 130057"
#[1] "ncells post-min_pct_reads_in_peaks: 124448"
#[1] "ncells post-max_blacklist_ratio: 124448"
#[1] "ncells post-max_nucleosome_signal: 124448"
#[1] "ncells post-min_TSS.enrichment: 124448"
saveRDS(mo_object_2, mo_object_2_filtered_loc)

# read the first object
mo_object_1 <- readRDS(mo_object_1_loc)
# process
mo_object_1 <- process_signac_object(mo_object_1, arc_metadata, fragments_loc = fragments_loc)
mo_object_1 <- qc_signac_object(mo_object_1)
#[1] "ncells pre QC: 184692"
#[1] "ncells post-min_nCount_peaks QC: 168673"
#[1] "ncells post-max_nCount_peaks: 159033"
#[1] "ncells post-min_pct_reads_in_peaks: 142113"
#[1] "ncells post-max_blacklist_ratio: 142113"
#[1] "ncells post-max_nucleosome_signal: 142113"
#[1] "ncells post-min_TSS.enrichment: 142113"
saveRDS(mo_object_1, mo_object_1_filtered_loc)

# read the first object
mo_object_3 <- readRDS(mo_object_3_loc)
# process
mo_object_3 <- process_signac_object(mo_object_3, arc_metadata, fragments_loc = fragments_loc)
mo_object_3 <- qc_signac_object(mo_object_3)h
#[1] "ncells pre QC: 143544"
#[1] "ncells post-min_nCount_peaks QC: 116650"
#[1] "ncells post-max_nCount_peaks: 115266"
#[1] "ncells post-min_pct_reads_in_peaks: 98627"
#[1] "ncells post-max_blacklist_ratio: 98627"
#[1] "ncells post-max_nucleosome_signal: 98627"
#[1] "ncells post-min_TSS.enrichment: 98627"
saveRDS(mo_object_3, mo_object_3_filtered_loc)

# read the first object
mo_object_4 <- readRDS(mo_object_4_loc)
# process
mo_object_4 <- process_signac_object(mo_object_4, arc_metadata, fragments_loc = fragments_loc)
mo_object_4 <- qc_signac_object(mo_object_4)
#[1] "ncells pre QC: 158193"
#[1] "ncells post-min_nCount_peaks QC: 136293"
#[1] "ncells post-max_nCount_peaks: 135285"
#[1] "ncells post-min_pct_reads_in_peaks: 121580"
#[1] "ncells post-max_blacklist_ratio: 121580"
#[1] "ncells post-max_nucleosome_signal: 121580"
#[1] "ncells post-min_TSS.enrichment: 121580"
saveRDS(mo_object_4, mo_object_4_filtered_loc)

# add more metadata
mo_object_1 <- AddMetaData(mo_object_1, rna_metadata['cell_type_final'], 'cell_type_final')
mo_object_1 <- AddMetaData(mo_object_1, rna_metadata['sample_final'], 'sample_final')
mo_object_1 <- AddMetaData(mo_object_1, rna_metadata['final_condition'], 'final_condition')
mo_object_1 <- AddMetaData(mo_object_1, rna_metadata['soup_status'], 'soup_status')
mo_object_2 <- AddMetaData(mo_object_2, rna_metadata['cell_type_final'], 'cell_type_final')
mo_object_2 <- AddMetaData(mo_object_2, rna_metadata['sample_final'], 'sample_final')
mo_object_2 <- AddMetaData(mo_object_2, rna_metadata['final_condition'], 'final_condition')
mo_object_2 <- AddMetaData(mo_object_2, rna_metadata['soup_status'], 'soup_status')
mo_object_3 <- AddMetaData(mo_object_3, rna_metadata['cell_type_final'], 'cell_type_final')
mo_object_3 <- AddMetaData(mo_object_3, rna_metadata['sample_final'], 'sample_final')
mo_object_3 <- AddMetaData(mo_object_3, rna_metadata['final_condition'], 'final_condition')
mo_object_3 <- AddMetaData(mo_object_3, rna_metadata['soup_status'], 'soup_status')
mo_object_4 <- AddMetaData(mo_object_4, rna_metadata['cell_type_final'], 'cell_type_final')
mo_object_4 <- AddMetaData(mo_object_4, rna_metadata['sample_final'], 'sample_final')
mo_object_4 <- AddMetaData(mo_object_4, rna_metadata['final_condition'], 'final_condition')
mo_object_4 <- AddMetaData(mo_object_4, rna_metadata['soup_status'], 'soup_status')

# merge per celltype
mo_all_per_celltype <- merge_signac_per_celltypes(c(mo_object_1, mo_object_2, mo_object_3, mo_object_4))
# make the names posix safe
names(mo_all_per_celltype) <- make_celltypes_safe(names(mo_all_per_celltype))
# save result
saveRDS(mo_all_per_celltype, '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_percelltype_1_64.rds')

# we should also check the cell type annotation. We can not do this for the full data, so let's do it for each separately
# remove doublets
mo_object_1 <- mo_object_1[, !is.na(mo_object_1@meta.data[['soup_status']]) & mo_object_1@meta.data[['soup_status']] == 'singlet']
ncol(mo_object_1)
# [1] 129983
# now do the dimreduc and clustering
mo_object_1 <- signac_dimreduc_and_cluster(mo_object_1)
# add lower classification celltypes
mo_object_1@meta.data[['cell_type_final_lowerres']] <- as.vector(unlist(ref10xmo_predictions_to_lower_res_mapping()[mo_object_1@meta.data[['cell_type_final']]]))
# plot
plot_grid(
  DimPlot(mo_object_1, group.by = 'seurat_clusters') + theme(legend.position = "none"),
  DimPlot(mo_object_1, group.by = 'cell_type_final_lowerres') + scale_color_manual(values = get_color_coding_dict()),
  nrow = 2
)
saveRDS(mo_object_1, mo_object_1_clustered_loc)
# remove doublets
mo_object_2 <- mo_object_2[, !is.na(mo_object_2@meta.data[['soup_status']]) & mo_object_2@meta.data[['soup_status']] == 'singlet']
ncol(mo_object_2)
# [1] 116330
mo_object_2 <- signac_dimreduc_and_cluster(mo_object_2)
mo_object_2@meta.data[['cell_type_final_lowerres']] <- as.vector(unlist(ref10xmo_predictions_to_lower_res_mapping()[mo_object_2@meta.data[['cell_type_final']]]))
plot_grid(
  DimPlot(mo_object_2, group.by = 'seurat_clusters') + theme(legend.position = "none"),
  DimPlot(mo_object_2, group.by = 'cell_type_final_lowerres') + scale_color_manual(values = get_color_coding_dict()),
  nrow = 2
)
saveRDS(mo_object_2, mo_object_2_clustered_loc)
