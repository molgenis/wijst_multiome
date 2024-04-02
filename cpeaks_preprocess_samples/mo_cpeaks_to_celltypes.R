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
  ncol(signac_object)
  print('ncells pre QC:', str(ncol(signac_object)))
  signac_object <- subset(x = signac_object, subset = nCount_peaks > min_nCount_peaks)
  print('ncells post-min_nCount_peaks QC:', str(ncol(signac_object)))
  signac_object <- subset(x = signac_object, subset = nCount_peaks < max_nCount_peaks)
  print('ncells post-max_nCount_peaks:', str(ncol(signac_object)))
  signac_object <- subset(x = signac_object, subset = pct_reads_in_peaks >= min_pct_reads_in_peaks)
  print('ncells post-min_pct_reads_in_peaks:', str(ncol(signac_object)))
  signac_object <- subset(x = signac_object, subset = blacklist_ratio < max_blacklist_ratio)
  print('ncells post-max_blacklist_ratio:', str(ncol(signac_object)))
  signac_object <- subset(x = signac_object, subset = nucleosome_signal < max_nucleosome_signal)
  print('ncells post-max_nucleosome_signal:', str(ncol(signac_object)))
  signac_object <- subset(x = signac_object, subset = TSS.enrichment > min_TSS.enrichment)
  print('ncells post-min_TSS.enrichment:', str(ncol(signac_object)))
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
mo_object_2 <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_unfiltered17_32.rds'

# load metadata
arc_metadata <- read.table(arc_metadata_loc, header = T, sep = '\t')
rownames(arc_metadata) <- arc_metadata[['barcode_lane']]

# read the first object
mo_object_2 <- readRDS(mo_object_2)
# process
mo_object_2 <- process_signac_object(mo_object_2, arc_metadata, fragments_loc = fragments_loc)
