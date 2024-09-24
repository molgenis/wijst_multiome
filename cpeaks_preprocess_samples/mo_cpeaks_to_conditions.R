#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_cpeaks_to_conditions.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

# for data
library(Seurat)
library(Signac)
# for plots
library(ggplot2)
library(cowplot)
# for upset plot of peak sharing
library(UpSetR)


####################
# Functions        #
####################


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
  
  # now for condition specific
  color_coding[["24hCA_Bulk"]] <- "black"
  color_coding[["24hCA Bulk"]] <- "black"
  color_coding[["24hCA_CD4T"]] <- "#153057"
  color_coding[["24hCA_CD8T"]] <- "#009DDB"
  color_coding[["24hCA monocyte"]] <- "#EDBA1B"
  color_coding[["24hCA_monocyte"]] <- "#EDBA1B"
  color_coding[["24hCA NK"]] <- "#E64B50"
  color_coding[["24hCA_NK"]] <- "#E64B50"
  color_coding[["24hCA B"]] <- "#71BC4B"
  color_coding[["24hCA_B"]] <- "#71BC4B"
  color_coding[["24hCA DC"]] <- "#965EC8"
  color_coding[["24hCA_DC"]] <- "#965EC8"
  color_coding[["24hCA CD4+ T"]] <- "#153057"
  color_coding[["24hCA CD8+ T"]] <- "#009DDB"
  color_coding[["24hCA HSPC"]] <- "#009E94"
  color_coding[["24hCA_HSPC"]] <- "#009E94"
  color_coding[["24hCA platelet"]] <- "#9E1C00"
  color_coding[["24hCA_platelet"]] <- "#9E1C00"
  color_coding[["24hCA plasmablast"]] <- "#DB8E00"
  color_coding[["24hCA_plasmablast"]] <- "#DB8E00"
  color_coding[["24hCA other T"]] <- "#FF63B6"
  color_coding[["24hCA_T_other"]] <- "#FF63B6"
  color_coding[["24hCA T-other"]] <- "#FF63B6"
  color_coding[["24hCA_T-other"]] <- "#FF63B6"
  color_coding[["24hCA_hemapoietic_stem"]] <- "#8B8000"
  color_coding[["24hCA hemapoietic stem"]] <- "#8B8000"
  color_coding[["24hCA hemapoietic-stem"]] <- "#8B8000"
  color_coding[["24hCA_hemapoietic-stem"]] <- "#8B8000"
  pct_whitening=40
  color_coding[["UT_Bulk"]] <- colorRampPalette(c(color_coding[["24hCA_Bulk"]], "white"))(100)[pct_whitening]
  color_coding[["UT Bulk"]] <- colorRampPalette(c(color_coding[["24hCA Bulk"]], "white"))(100)[pct_whitening]
  color_coding[["UT_CD4T"]] <- colorRampPalette(c(color_coding[["24hCA_CD4T"]], "white"))(100)[pct_whitening]
  color_coding[["UT_CD4T"]] <- colorRampPalette(c(color_coding[["24hCA_CD4T"]], "white"))(100)[pct_whitening]
  color_coding[["UT monocyte"]] <- colorRampPalette(c(color_coding[["24hCA monocyte"]], "white"))(100)[pct_whitening]
  color_coding[["UT_monocyte"]] <- colorRampPalette(c(color_coding[["24hCA_monocyte"]], "white"))(100)[pct_whitening]
  color_coding[["UT NK"]] <- colorRampPalette(c(color_coding[["24hCA NK"]], "white"))(100)[pct_whitening]
  color_coding[["UT_NK"]] <- colorRampPalette(c(color_coding[["24hCA_NK"]], "white"))(100)[pct_whitening]
  color_coding[["UT B"]] <- colorRampPalette(c(color_coding[["24hCA B"]], "white"))(100)[pct_whitening]
  color_coding[["UT_B"]] <- colorRampPalette(c(color_coding[["24hCA_B"]], "white"))(100)[pct_whitening]
  color_coding[["UT DC"]] <- colorRampPalette(c(color_coding[["24hCA DC"]], "white"))(100)[pct_whitening]
  color_coding[["UT_DC"]] <- colorRampPalette(c(color_coding[["24hCA_DC"]], "white"))(100)[pct_whitening]
  color_coding[["UT CD4+ T"]] <- colorRampPalette(c(color_coding[["24hCA CD4+ T"]], "white"))(100)[pct_whitening]
  color_coding[["UT CD8+ T"]] <- colorRampPalette(c(color_coding[["24hCA CD8+ T"]], "white"))(100)[pct_whitening]
  color_coding[["UT HSPC"]] <- colorRampPalette(c(color_coding[["24hCA HSPC"]], "white"))(100)[pct_whitening]
  color_coding[["UT_HSPC"]] <- colorRampPalette(c(color_coding[["24hCA_HSPC"]], "white"))(100)[pct_whitening]
  color_coding[["UT platelet"]] <- colorRampPalette(c(color_coding[["24hCA platelet"]], "white"))(100)[pct_whitening]
  color_coding[["UT_platelet"]] <- colorRampPalette(c(color_coding[["24hCA_platelet"]], "white"))(100)[pct_whitening]
  color_coding[["UT plasmablast"]] <- colorRampPalette(c(color_coding[["24hCA plasmablast"]], "white"))(100)[pct_whitening]
  color_coding[["UT_plasmablast"]] <- colorRampPalette(c(color_coding[["24hCA_plasmablast"]], "white"))(100)[pct_whitening]
  color_coding[["UT other T"]] <- colorRampPalette(c(color_coding[["24hCA other T"]], "white"))(100)[pct_whitening]
  color_coding[["UT_T_other"]] <- colorRampPalette(c(color_coding[["24hCA_T_other"]], "white"))(100)[pct_whitening]
  color_coding[["UT T-other"]] <- colorRampPalette(c(color_coding[["24hCA T-other"]], "white"))(100)[pct_whitening]
  color_coding[["UT_T-other"]] <- colorRampPalette(c(color_coding[["24hCA_T-other"]], "white"))(100)[pct_whitening]
  color_coding[["UT_hemapoietic_stem"]] <- colorRampPalette(c(color_coding[["24hCA_hemapoietic_stem"]], "white"))(100)[pct_whitening]
  color_coding[["UT hemapoietic stem"]] <- colorRampPalette(c(color_coding[["24hCA hemapoietic stem"]], "white"))(100)[pct_whitening]
  color_coding[["UT hemapoietic-stem"]] <- colorRampPalette(c(color_coding[["24hCA hemapoietic-stem"]], "white"))(100)[pct_whitening]
  color_coding[["UT_hemapoietic-stem"]] <- colorRampPalette(c(color_coding[["24hCA_hemapoietic-stem"]], "white"))(100)[pct_whitening]
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
  # with the conditions
  label_dict[["UT_Bulk"]] <- "UT bulk-like"
  label_dict[["24hCA_Bulk"]] <- "24hCA bulk-like"
  label_dict[["UT_CD4T"]] <- "UT CD4+ T"
  label_dict[["24hCA_CD4T"]] <- "24hCA CD4+ T"
  label_dict[["UT_CD8T"]] <- "UT CD8+ T"
  label_dict[["24hCA_CD8T"]] <- "24hCA CD8+ T"
  label_dict[["UT_monocyte"]] <- "UT monocyte"
  label_dict[["24hCA_monocyte"]] <- "24hCA monocyte"
  label_dict[["UT_NK"]] <- "UT NK"
  label_dict[["24hCA_NK"]] <- "24hCA NK"
  label_dict[["UT_B"]] <- "UT B"
  label_dict[["24hCA_B"]] <- "24hCA B"
  label_dict[["UT_DC"]] <- "UT DC"
  label_dict[["24hCA_DC"]] <- "24hCA DC"
  label_dict[["UT_HSPC"]] <- "UT HSPC"
  label_dict[["24hCA HSPC"]] <- "24hCA HSPC"
  label_dict[["UT_plasmablast"]] <- "UT plasmablast"
  label_dict[["24hCA_plasmablast"]] <- "24hCA plasmablast"
  label_dict[["UT_platelet"]] <- "UT platelet"
  label_dict[["24hCA_platelet"]] <- "24hCA platelet"
  label_dict[["UT_T_other"]] <- "UT other T"
  label_dict[["24hCA_T_other"]] <- "24hCA other T"
  label_dict[["UT_T-other"]] <- "UT other T"
  label_dict[["24hCA_T-other"]] <- "24hCA other T"
  label_dict[["UT_hemapoietic_stem"]] <- "UT hemapoietic stem"
  label_dict[["24hCA_hemapoietic_stem"]] <- "24hCA hemapoietic stem"
  label_dict[["UT_hemapoietic-stem"]] <- "UT hemapoietic stem"
  label_dict[["24hCA_hemapoietic-stem"]] <- "24hCA hemapoietic stem"
  return(label_dict)
}


read_barcode_and_lane <- function(seurat_object) {
  # do the split first
  seurat_object_rowsnames_split <- strsplit(colnames(seurat_object), split = '_')
  # now do a list apply
  df_per_barcode <- lapply(seurat_object_rowsnames_split, FUN = function(x) {
    data.frame(lane = paste(x[2], x[3], sep = '_'), barcode = x[1])
  })
  # merge all of them together
  extra_metadata <- do.call('rbind', df_per_barcode)
  # set the rownames to be the original ones
  rownames(extra_metadata) <- colnames(seurat_object)
  # now add the extra data we have
  seurat_object <- AddMetaData(seurat_object, extra_metadata)
  return(seurat_object)
}


#' add the inflammation assignments  to the Seurat object
#' 
#' @param seurat_object The Seurat object to add the inflammation status to
#' @param sample_sheet The sample sheet containing lanes, participants and inflammation statuses
#' @param seurat_lane_column The column in the Seurat metadata denoting the 10x lane
#' @param sheet_lane_column The column in the sample sheet denoting the 10x lane
#' @param seurat_participant_column The column in the Seurat metadata denoting the participant assignment
#' @param sheet_participants_column The column in the sample sheet containing the participants per lane
#' @param seurat_inflammation_column The column in the Seurat metadata to add the inflammation status in
#' @param sheet_inflammation_column The column in the sample sheet containing the inflammation statuses per lane
#' @returns the Seurat object with the inflammation status added
#' lpmcv2 <- add_inflammation_status(lpmcv2, sample_sheet)
add_inflammation_status <- function(seurat_object, sample_sheet, seurat_lane_column='lane', sheet_lane_column='lane', seurat_participant_column='soup_final_sample_assignment', sheet_participants_column='genoid', seurat_inflammation_column='inflammation_status', sheet_inflammation_column='inflammation_status') {
  # create a mapping of lane+sample to inflammation status
  mapping_per_lane_list <- list()
  for (i in 1:nrow(sample_sheet)) {
    # extract lane
    lane <- sample_sheet[i, sheet_lane_column]
    # extract the participants
    participant <- sample_sheet[i, sheet_participants_column]
    # and the inflammation condition
    condition <- sample_sheet[i, sheet_inflammation_column]
    # if not set, set to unknown
    if (is.null(condition)) {
      condition <- 'unknown'
    }
    # subset to the barcodes which have are this lane and participant
    barcodes_match <- rownames(seurat_object@meta.data[!is.na(seurat_object@meta.data[[seurat_lane_column]]) &
                                                         seurat_object@meta.data[[seurat_lane_column]] == lane &
                                                         !is.na(seurat_object@meta.data[[seurat_participant_column]]) &
                                                         seurat_object@meta.data[[seurat_participant_column]] == participant, ])
    # only add if there are matching barcodes
    if (length(barcodes_match) > 0) {
      # create dataframe
      df_lane_part <- data.frame(barcode = barcodes_match, condition = rep(condition, times = length(barcodes_match)))
      # set the colname to be the one we chose
      colnames(df_lane_part) <- c('barcode', seurat_inflammation_column)
      # then add to the list
      mapping_per_lane_list[[paste(lane, participant, sep = ':')]] <- df_lane_part
    }
  }
  # now merge all together
  mapping_all <- do.call('rbind', mapping_per_lane_list)
  # set the barcode as rownames
  rownames(mapping_all) <- mapping_all[['barcode']]
  # finally add to the object
  seurat_object <- AddMetaData(seurat_object, mapping_all[seurat_inflammation_column])
  return(seurat_object)
}


add_inflammation_status_each_object <- function(seurat_object_list, sample_sheet, seurat_lane_column='lane', sheet_lane_column='lane', seurat_participant_column='soup_final_sample_assignment', sheet_participants_column='genoid', seurat_inflammation_column='inflammation_status', sheet_inflammation_column='inflammation_status') {
  # let's go over each object
  for (object_name in names(seurat_object_list)) {
    # do the condition assignment
    seurat_object_list[[object_name]] <- add_inflammation_status(
      seurat_object_list[[object_name]],
      sample_sheet=sample_sheet, 
      seurat_lane_column=seurat_lane_column,
      sheet_lane_column=sheet_lane_column, 
      seurat_participant_column=seurat_participant_column, 
      sheet_participants_column=sheet_participants_column, 
      seurat_inflammation_column=seurat_inflammation_column, 
      sheet_inflammation_column=sheet_inflammation_column
    )
  }
  return(seurat_object_list)
}


split_by_column <- function(seurat_object_list, split_column) {
  # save split data
  split_seurat_object_list <- list()
  # let's go over each object
  for (object_name in names(seurat_object_list)) {
    # do the condition assignment
    seurat_object <- seurat_object_list[[object_name]]
    # get the possible values
    possible_values <- unique(seurat_object@meta.data[[split_column]])
    # exclude the NA values
    possible_values <- possible_values[!is.na(possible_values)]
    # check each value
    for (possible_value in possible_values) {
      # get seurat object for that possible value
      seurat_object_value <- seurat_object[, !is.na(seurat_object@meta.data[[split_column]]) &
                                             seurat_object@meta.data[[split_column]] == possible_value]
      # add this to the list
      split_seurat_object_list[[paste(possible_value, object_name, sep = '_')]] <- seurat_object_value
    }
  }
  return(split_seurat_object_list)
}


get_peaks_to_bed <- function(exp_per_group) {
  bed_per_group <- list()
  # check each of the columns
  for (identity in colnames(exp_per_group)) {
    print(identity)
    # subset to that identity
    peaks_ident <- exp_per_group[, c(identity), drop = F]
    # get the locations
    positions <- data.frame(do.call('rbind', strsplit(rownames(peaks_ident), '-')))
    # set the colnames properly to bed format
    colnames(positions) <- c('#chrom', 'start', 'end')
    # add the rownames themselves as the name
    positions[['name']] <- rownames(peaks_ident)
    # and the counts as score
    positions[['exp']] <- log10(peaks_ident[, c(identity)])
    # make the start and stop numeric
    positions[['start']] <- as.numeric(positions[['start']])
    positions[['end']] <- as.numeric(positions[['end']])
    # add bed to list
    bed_per_group[[identity]] <- positions
  }
  return(bed_per_group)
}



summarize_peak_info <- function(peaks_per_celltype, output_prepend='mo_peaks_', output_append='.bed') {
  # check each cell type
  for (cell_type in names(peaks_per_celltype)) {
    print(cell_type)
    # get this celltype object
    peaks_object <- peaks_per_celltype[[cell_type]]
    # calculate average expression
    avg_peaks <- AverageExpression(peaks_object)[[1]]
    # calculate average expression
    sum_peaks <- AggregateExpression(peaks_object)[[1]]
    # calculate pct exp
    counts_as_row_sparse <- as(peaks_object@assays$peaks@counts, "RsparseMatrix") # turn into row-wise sparse matrix
    non_zero_cell_nr <- as.vector(unlist(rowSums(counts_as_row_sparse != 0))) # do a rowsum on the T/F value you get from the zero-or-not comparison
    pct_peaks <- non_zero_cell_nr / ncol(peaks_object@assays$peaks@counts) # the percentage expressed is that number of cells divided by the total number
    # turn the exp per group into a bed
    peaks_bed <- get_peaks_to_bed(sum_peaks)[[1]]
    # add the other info
    peaks_bed[['avg']] <-  unlist(as.vector(avg_peaks[, 1]))
    peaks_bed[['ncell']] <-  ncol(peaks_object)
    peaks_bed[['pct_exp']] <- pct_peaks
    # and sort
    peaks_bed <- peaks_bed[order(peaks_bed[['#chrom']], peaks_bed[['start']], peaks_bed[['end']]), ]
    # create the output location
    output_loc <- paste(output_prepend, gsub(' ', '_', cell_type), output_append, sep = '')
    # write the result
    write.table(peaks_bed, output_loc, row.names = F, col.names = T, quote = F, sep = '\t')
  }
  return(0)
}


get_peaks_sharing_from_beds <- function(output_prepend='mo_peaks_', output_append='.bed', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), filter_column=NULL, filter_value=NULL, id_column='name') {
  # we'll do each identity, and store the result in per ident in a list
  peaks_per_identity <- list()
  # check each of the columns
  for (cell_type in cell_types) {
    # create the output location
    output_loc <- paste(output_prepend, cell_type, output_append, sep = '')
    # read the result
    peaks_bed <- read.table(output_loc, header = T, sep = '\t', comment.char='', check.names = F)
    # filter by filter if present
    if (!is.null(filter_column) & !is.null(filter_value)) {
      peaks_bed <- peaks_bed[peaks_bed[[filter_column]] > filter_value, ]
    }
    # extract the peaks
    peaks_per_identity[[cell_type]] <- peaks_bed[[id_column]]
  }
  return(peaks_per_identity)
}

plot_peak_sharing_from_beds <- function(output_prepend='mo_peaks_', output_append='.bed', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), filter_column=NULL, filter_value=NULL, id_column='name', use_color_dict=T, use_label_dict=T, precalculated_peaks=NULL) {
  # we will use peaks per ct
  peaks_per_ct <- NULL
  # if we have them already calculated, we will use the filtered peaks
  if (!is.null(precalculated_peaks)) {
    peaks_per_ct <- precalculated_peaks
  }
  # else get the peak sharing from the beds
  else {
    peaks_per_ct <- get_peaks_sharing_from_beds(output_prepend=output_prepend, output_append=output_append, cell_types=cell_types, filter_column=filter_column, filter_value=filter_value, id_column=id_column)
  }
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
  upset(fromList(peaks_per_ct), order.by = 'freq', nsets = length(peaks_per_ct), queries = queries, sets.bar.color=sets.bar.color, nintersects = length(peaks_per_ct)*length(peaks_per_ct))
  #upset(fromList(peaks_per_ct), order.by = 'freq', nsets = length(peaks_per_ct), sets.bar.color=sets.bar.color, nintersects = length(peaks_per_ct)*length(peaks_per_ct))
  #return(peaks_per_ct)
}


write_peak_mtx <- function(signac_object, output_folder) {
  # features
  features_gz <- gzfile(paste(output_folder, 'features.tsv.gz', sep = ''))
  write.table(data.frame(x = rownames(signac_object@assays$peaks@counts)), features_gz, row.names = F, col.names = F, quote = F)
  # barcodes
  barcodes_gz <- gzfile(paste(output_folder, 'barcodes.tsv.gz', sep = ''))
  write.table(data.frame(x = colnames(signac_object@assays$peaks@counts)), barcodes_gz, row.names = F, col.names = F, quote = F)
  # metadata
  metadata_gz <- gzfile(paste(output_folder, 'metadata.tsv.gz', sep = ''))
  write.table(cbind(data.frame(bc = rownames(signac_object@meta.data)), signac_object@meta.data), metadata_gz, row.names = F, col.names = T, quote = F, sep = '\t')
  # and finally the count matrix
  counts_gz <- (paste(output_folder, 'matrix.mtx', sep = ''))
  writeMM(signac_object@assays$peaks@counts, counts_gz)
  # needs to be bgzipped manually
  return(0)
}


write_peak_matrices <- function(signac_object_list, output_folder) {
  # check each object
  for (entry in names(signac_object_list)) {
    # get the object
    signac_object <- signac_object_list[[entry]]
    # get the full output directory
    full_output_folder <- paste(output_folder, '/', entry, '/', sep = '')
    # create the directory if not present already
    dir.create(full_output_folder, showWarnings = F, recursive = T)
    # write the actual stuff
    write_peak_mtx(signac_object, full_output_folder)
  }
  return(0)
}

####################
# Main Code        #
####################

# location of the condition assignment
condition_assignment_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_monocyte_based_condition_numbers.tsv'

# location of the cell type objects
cell_type_objects_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_percelltypemajor_1_80.rds'
cell_type_objects_wstatus_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_percelltypemajor_wstatus_1_80.rds'

# read the object
cell_type_objects <- readRDS(cell_type_objects_loc)

# read the conditions
condition_assignments <- read.table(condition_assignment_loc, header = T, sep = '\t')

# add barcodes back
for(cell_type in names(cell_type_objects)) {
  cell_type_objects[[cell_type]] <- read_barcode_and_lane(cell_type_objects[[cell_type]])
}

# get the assignment matrices
correlation_mapping_per_barcode_all <- read.table('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_sample_matched_vs_all.tsv', header = T, sep = '\t')
# set barcodes and remove data we already have
rownames(correlation_mapping_per_barcode_all) <- correlation_mapping_per_barcode_all[['barcode_lane']]
correlation_mapping_per_barcode_all[, c('lane', 'barcode_lane', 'barcode', 'barcode_original')] <- NULL
# now let's get the souporcell data specifically, which would be the same for all and per-lane
soup_only <- correlation_mapping_per_barcode_all[, setdiff(colnames(correlation_mapping_per_barcode_all), c('best_match_sample', 'second_match_sample', 'best_match_correlation', 'second_match_correlation'))]
# and the correlation data
correlations_unconfined <- correlation_mapping_per_barcode_all[, c('best_match_sample', 'second_match_sample', 'best_match_correlation', 'second_match_correlation')]
# add the confined sample
for(cell_type in names(cell_type_objects)) {
  cell_type_objects[[cell_type]] <- AddMetaData(cell_type_objects[[cell_type]], correlations_unconfined[, colnames(correlations_unconfined)])
}

# add the conditions of the original sheet
cell_type_objects <- add_inflammation_status_each_object(cell_type_objects, condition_assignments, seurat_participant_column='best_match_sample', sheet_participants_column = 'sample_final', seurat_inflammation_column = 'inflammation_sheet', sheet_inflammation_column = 'condition')
cell_type_objects <- add_inflammation_status_each_object(cell_type_objects, condition_assignments, seurat_participant_column='best_match_sample', sheet_participants_column = 'sample', seurat_inflammation_column = 'inflammation_prev', sheet_inflammation_column = 'cond_prev')
# rename CA in the original assignment, and set a consensus one
for (cell_type in names(cell_type_objects)) {
  cell_type_object <- cell_type_objects[[cell_type]]
  cell_type_object@meta.data[!is.na(cell_type_object@meta.data[['inflammation_sheet']]) &
                               cell_type_object@meta.data[['inflammation_sheet']] == '24hCa', 'inflammation_sheet'] <- '24hCA'
  # now also set the final inflammation assignment
  cell_type_object@meta.data[['inflammation_final']] <- cell_type_object@meta.data[['inflammation_sheet']]
  cell_type_object@meta.data[is.na(cell_type_object@meta.data[['inflammation_final']]), 'inflammation_final'] <- cell_type_object@meta.data[is.na(cell_type_object@meta.data[['inflammation_final']]), 'inflammation_prev']
  # add back to list
  cell_type_objects[[cell_type]] <- cell_type_object
}

# now split by inflammation status
cell_type_objects_condition <- split_by_column(cell_type_objects, 'inflammation_final')
# and make the beds
summarize_peak_info(cell_type_objects_condition, output_prepend = '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to80_')

# read the cell type annotations
celltype_annotations <- read.table('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cell_type_assignment/azimuth/10x_multiome_PBMCs/mo_azimuth_ct_10xmultiome.tsv', header = T, sep = '\t')
# subset to those
rna_infos <- celltype_annotations[['barcode']]
# also do this filtering for the matched ATAC and RNA
for(ct_cond in names(cell_type_objects_condition)) {
  # extract
  celltype_cond <- cell_type_objects_condition[[ct_cond]]
  # subset to ones we have the celltype annotations for, which are the matched RNA/ATAC ones
  celltype_cond <- celltype_cond[, colnames(celltype_cond) %in% rna_infos]
  # put in a new list
  celltype_cond_list <- list()
  celltype_cond_list[[ct_cond]] <- celltype_cond
  # write this dummy list
  summarize_peak_info(celltype_cond_list, output_prepend = '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to80_wrna_')
}

# do 
plot_peak_sharing_from_beds('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to64_UT_', filter_column = 'exp', filter_value = 1) # minimal ten counts per cell type
plot_peak_sharing_from_beds('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to64_UT_', filter_column = 'avg', filter_value = .1) # on average expressed in one out of ten cells
plot_peak_sharing_from_beds('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to64_UT_', filter_column = 'pct_exp', filter_value = .1) # expressed in at least 10% of cells
plot_peak_sharing_from_beds('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to64_24hCA_', filter_column = 'exp', filter_value = 1) # minimal ten counts per cell type
plot_peak_sharing_from_beds('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to64_24hCA_', filter_column = 'avg', filter_value = .1) # on average expressed in one out of ten cells
plot_peak_sharing_from_beds('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to64_24hCA_', filter_column = 'pct_exp', filter_value = .1) # expressed in at least 10% of cells

# now combine the UT and 24hCA
mo_peaks_lane1to64_UT_avgfilter <- get_peaks_sharing_from_beds(output_prepend = '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to64_UT_', filter_column = 'avg', filter_value = .1)
mo_peaks_lane1to64_24hCA_avgfilter <- get_peaks_sharing_from_beds(output_prepend = '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to64_24hCA_', filter_column = 'avg', filter_value = .1)
# but first rename the names in the named list
names(mo_peaks_lane1to64_UT_avgfilter) <- paste('UT', names(mo_peaks_lane1to64_UT_avgfilter), sep = '_')
names(mo_peaks_lane1to64_24hCA_avgfilter) <- paste('24hCA', names(mo_peaks_lane1to64_24hCA_avgfilter), sep = '_')
mo_peaks_lane1to64_both_avgfilter <- c(mo_peaks_lane1to64_UT_avgfilter, mo_peaks_lane1to64_24hCA_avgfilter)
plot_peak_sharing_from_beds(precalculated_peaks = mo_peaks_lane1to64_both_avgfilter, use_color_dict = T, use_label_dict = T)

# also do with the minimally expressed in 0.01 cells
mo_peaks_lane1to64_UT_minfilter <- get_peaks_sharing_from_beds(output_prepend = '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to64_UT_', filter_column = 'pct_exp', filter_value = .01)
mo_peaks_lane1to64_24hCA_minfilter <- get_peaks_sharing_from_beds(output_prepend = '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to64_24hCA_', filter_column = 'pct_exp', filter_value = .01)
names(mo_peaks_lane1to64_UT_minfilter) <- paste('UT', names(mo_peaks_lane1to64_UT_minfilter), sep = '_')
names(mo_peaks_lane1to64_24hCA_minfilter) <- paste('24hCA', names(mo_peaks_lane1to64_24hCA_minfilter), sep = '_')
mo_peaks_lane1to64_both_minfilter <- c(mo_peaks_lane1to64_UT_minfilter, mo_peaks_lane1to64_24hCA_minfilter)
plot_peak_sharing_from_beds(precalculated_peaks = mo_peaks_lane1to64_both_minfilter, use_color_dict = T, use_label_dict = T)

# save the object
saveRDS(cell_type_objects, cell_type_objects_wstatus_loc)

# create the peak matrices
write_peak_matrices(cell_type_objects, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/')

