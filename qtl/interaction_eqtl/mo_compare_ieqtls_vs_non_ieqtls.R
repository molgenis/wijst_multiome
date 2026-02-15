#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_compare_ieqtls_vs_non_ieqtls.R
# Function: compare the characteristics of the eQTLs that are also i-eQTLs versus those that are not
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(IRanges)
library(ggplot2)
library(cowplot)
library(stringr)


####################
# Functions        #
####################

get_qtls_per_celltype_limix <- function(qtl_output_loc, output_file='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz', gene_column='feature_id', significance_column='feature_q_value', significance_cutoff=0.05, nominal_cutoff_column='pval_nominal_threshold_global', nominal_significance_column='p_value', verbose=T) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(qtl_output_loc, full.names = F, recursive = F)
  # we will store the results in a list for now
  qtls_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the full path
    full_cell_type_path <- paste(qtl_output_loc, '/', cell_type, '/', output_file, sep = '')
    # log if requested
    if (verbose) {
      print(paste('reading', full_cell_type_path))
    }
    # read the file
    cell_type_output <- read.table(full_cell_type_path, sep = '\t', header = T)
    # filter the results on significance
    if (!is.null(significance_column)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering by significance column', nrow(cell_type_output)))
      }
      # filter
      cell_type_output <- cell_type_output[
        !is.na(cell_type_output[[significance_column]]) &
          cell_type_output[[significance_column]] < significance_cutoff, 
      ]
      if (verbose) {
        print(paste('variant+phenotype after filtering by significance column', nrow(cell_type_output)))
      }
    }
    if (!is.null(nominal_cutoff_column) & !is.null(nominal_significance_column)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering by nominal cutoff', nrow(cell_type_output)))
      }
      # filter
      cell_type_output <- cell_type_output[
        !is.na(cell_type_output[[nominal_cutoff_column]]) & 
          !is.na(cell_type_output[[nominal_significance_column]]) &
          cell_type_output[[nominal_significance_column]] <= cell_type_output[[nominal_cutoff_column]], 
      ]
      if (verbose) {
        print(paste('variant+phenotype after filtering by nominal cutoff', nrow(cell_type_output)))
      }
    }
    # add the celltype as a column
    cell_type_output[['cell_type']] <- cell_type
    # add to the list
    qtls_per_celltype[[cell_type]] <- cell_type_output
  }
  # turn into a dataframe
  return(qtls_per_celltype)
}


add_top_effect_annotation <- function(qtls_per_celltype, feature_column='feature_id', variant_column='snp_id', significance_column='p_value', decreasing=F) {
  # go through each of the cell types
  for (cell_type in names(qtls_per_celltype)) {
    # extract that table
    qtls_celltype <- qtls_per_celltype[[cell_type]]
    # order by the significancce
    qtls_celltype_ordered <- qtls_celltype[order(qtls_celltype[[significance_column]], decreasing = decreasing), ]
    # keep only the top effect
    qtls_celltype_top <- qtls_celltype_ordered[!duplicated(qtls_celltype_ordered[[feature_column]]), ]
    # and annotate in the original table if the variant-feature combination was the top one
    qtls_celltype[['is_top_variant']] <- paste(qtls_celltype[[variant_column]], qtls_celltype[[feature_column]]) %in% paste(qtls_celltype_top[[variant_column]], qtls_celltype_top[[feature_column]])
    # put that back in the list
    qtls_per_celltype[[cell_type]] <- qtls_celltype
  }
  return(qtls_per_celltype)
}


get_closest_flanks <- function(position_table, left_flank_column1, right_flank_column1, left_flank_column2, right_flank_column2) {
  # get the distance between left flanks
  dist_left_flank1_to_left_flank2 <- position_table[[left_flank_column1]] - position_table[[left_flank_column2]]
  # distance between the right flanks
  dist_right_flank1_to_right_flank2 <- position_table[[right_flank_column1]] - position_table[[right_flank_column2]]
  # distance between left flank 1 and right flank 2
  dist_left_flank1_to_right_flank2 <- position_table[[left_flank_column1]] - position_table[[right_flank_column2]]
  # distance between right flank 1 and left flank 2
  dist_right_flank1_to_left_flank2 <- position_table[[right_flank_column1]] - position_table[[left_flank_column2]]
  # put in a table for convenience sake
  distances_tbl <- data.table(
    'lf1_to_lf2' = dist_left_flank1_to_left_flank2, 
    'rf1_to_rf2' = dist_right_flank1_to_right_flank2, 
    'lf1_to_rf2' = dist_left_flank1_to_right_flank2, 
    'rf1_to_lf2' = dist_right_flank1_to_left_flank2
  )
  # add the minimum absolute distance
  distances_tbl[['min_dist']] <- apply(distances_tbl, 1, function(x) {
    return(min(abs(x)))
  })
  # but set this to zero if any of the flanks end in the bodies
  #              -----
  #                 ++++
  distances_tbl[(distances_tbl[['lf1_to_lf2']] < 0 & distances_tbl[['lf1_to_rf2']] > 0) |
                  #                   ----
                #                 ++++
                (distances_tbl[['rf1_to_lf2']] > 0 & distances_tbl[['lf1_to_rf2']] < 0) |
                  #                   ----
                #                 +++++++++
                (distances_tbl[['lf1_to_lf2']] > 0 & distances_tbl[['rf1_to_rf2']] < 0) |
                  #                 ---------
                #                   ++++
                (distances_tbl[['lf1_to_lf2']] < 0 & distances_tbl[['rf1_to_rf2']] > 0)
                , 'min_dist'] <- 0
  return(distances_tbl)
}

classify_variant_location <- function(qtl_tbl, variant_position_col='snp_position', feature_start_col='feature_start', feature_end_col='feature_end', strand_col='strand') {
  qtl_tbl[['variant_classification']] <- apply(qtl_tbl, 1, function(x) {
    # extract value
    snp_position <- as.numeric(x[variant_position_col])
    feature_start <- as.numeric(x[feature_start_col])
    feature_end <- as.numeric(x[feature_end_col])
    strand <- x[strand_col]
    # if we don't have position info
    if (is.na(snp_position) | is.na(feature_start) | is.na(feature_end)) {
      return(NA)
    }
    # check if the variant is in the gene body
    if ((snp_position >= feature_start) & (snp_position <= feature_end)) {
      return('in_gene_body')
    }
    else {
      if (is.na(strand)) {
        # without strand info, we cannot do anything
        return(NA)
      }
      else if (strand == '+') {
        # check if the variant is before the TSS
        if (snp_position < feature_start) {
          return('before_tss')
        }
        # or behind the gene body
        else if (snp_position > feature_end) {
          return('behind_gene_body')
        }
        else {
          return(NA)
        }
      }
      else if (strand == '-') {
        # if negative, the rules are the other way aroun
        if (snp_position < feature_start) {
          return('behind_gene_body')
        }
        else if (snp_position > feature_end) {
          return('before_tss')
        }
        else {
          return(NA)
        }
      }
      else {
        # if not plus or minus, we cannot do anything
        return(NA)
      }
    }
  })
  qtl_tbl[['variant_classification']] <- as.vector(unlist(qtl_tbl[['variant_classification']]))
  return(qtl_tbl)
}


add_introns_to_exons_info <- function(exons_annotations, feature_type_column='feature', gene_column='gene_id', start_column='start', end_column='end', strand_column='strand') {
  # filter to only keep exons
  exons_annotations <- exons_annotations[!is.na(exons_annotations[[feature_type_column]]) & exons_annotations[[feature_type_column]] == 'exon', ]
  # we'll need to do this separately for the positive and negative strand sets, so make a list to store these
  intron_annotations_per_strand <- list()
  # check each strand
  for (strand in c('-', '+')) {
    # use the start and end colunmns
    start_column_strand <- NULL
    end_column_strand <- NULL
    # these are different depending on the strand
    if (strand == '+') {
      start_column_strand <- start_column
      end_column_strand <- end_column
    }
    else if(strand == '-') {
      # start_column_strand <- end_column
      # end_column_strand <- start_column
      # apparently they are the same
      start_column_strand <- start_column
      end_column_strand <- end_column
    }
    # filter by strand
    exons_annotations_strand <- exons_annotations[!is.na(exons_annotations[[strand_column]]) & exons_annotations[[strand_column]] == strand, ]
    # order by gene and start position
    exons_annotations_strand <- exons_annotations_strand[
      order(
        exons_annotations_strand[[gene_column]], 
        exons_annotations_strand[[start_column_strand]], 
        exons_annotations_strand[[end_column_strand]]), 
    ]
    # take the end positions with the gene names
    intron_annotations <- cbind(
      # skipping the last row because that is the end of the last gene
      exons_annotations_strand[1 : (nrow(exons_annotations_strand) - 1), c(gene_column, end_column_strand)], 
      # with the starts of the next gene, skipping the first row because that is the start of the first gene
      exons_annotations_strand[2 : nrow(exons_annotations_strand), c(gene_column, start_column_strand)]
    )
    # set column names, as they are duplicated
    colnames(intron_annotations) <- c(gene_column, start_column, 'gene_id2', end_column)
    # remove where the gene names are different, as that is where one gene ends, and another starts
    intron_annotations <- intron_annotations[intron_annotations[[gene_column]] == intron_annotations[['gene_id2']], ]
    # remove any duplicates
    intron_annotations <- intron_annotations[!duplicated(intron_annotations), ]
    # remove where there is no space between the exons
    intron_annotations <- intron_annotations[intron_annotations[[end_column]] - intron_annotations[[start_column]] > 1, ]
    # remove the now superfluous gene_id2
    intron_annotations[['gene_id2']] <- NULL
    # add the strand info
    intron_annotations[[strand_column]] <- rep(strand, times = nrow(intron_annotations))
    # and that these are introns
    intron_annotations[[feature_type_column]] <- 'intron'
    # now add 1 to each start, so it does not overlap with the exons
    intron_annotations[[start_column]] <- intron_annotations[[start_column]] + 1
    # and substract 1 from each end, so it does not overlap with the exons
    intron_annotations[[end_column]] <- intron_annotations[[end_column]] - 1
    # add back the exon info as well
    intron_annotations <- rbind(
      intron_annotations, 
      exons_annotations_strand[, c(gene_column, start_column_strand, end_column_strand, strand_column, feature_type_column)]
    )
    # and order
    intron_annotations <- intron_annotations[
      order(
        intron_annotations[[gene_column]], 
        intron_annotations[[start_column_strand]], 
        intron_annotations[[end_column_strand]]), 
    ]
    # remove duplicates
    intron_annotations <- intron_annotations[!duplicated(intron_annotations), ]
    # add new column
    intron_annotations[['number_strand']] <- NA
    # add the exon and intron numbers
    intron_annotations[intron_annotations[[feature_type_column]] == 'exon', ][['number_strand']]   <- ave(1:nrow(intron_annotations[intron_annotations[[feature_type_column]] == 'exon', ]), intron_annotations[intron_annotations[[feature_type_column]] == 'exon', ][[gene_column]], FUN = function(x) { seq_along(x) })
    intron_annotations[intron_annotations[[feature_type_column]] == 'intron', ][['number_strand']]   <- ave(1:nrow(intron_annotations[intron_annotations[[feature_type_column]] == 'intron', ]), intron_annotations[intron_annotations[[feature_type_column]] == 'intron', ][[gene_column]], FUN = function(x) { seq_along(x) })
    # additionally, do the same, but take into consideration the strandedness
    if (strand == '+') {
      intron_annotations[['number']] <- intron_annotations[['number_strand']]
    }
    else {
      # where we reverse the positions so that the numbers are in the correct order in relation to being on the opposite strand
      intron_annotations <- intron_annotations[order(intron_annotations[[start_column_strand]], 
                                                     intron_annotations[[end_column_strand]], decreasing = T
      ), ]
      # but keep the gene order non-reversed
      intron_annotations <- intron_annotations[order(intron_annotations[[gene_column]]), ]
      # add that new number column
      intron_annotations[['number']] <- NA
      # then add that to the table
      intron_annotations[intron_annotations[[feature_type_column]] == 'exon', ][['number']]   <- ave(1:nrow(intron_annotations[intron_annotations[[feature_type_column]] == 'exon', ]), intron_annotations[intron_annotations[[feature_type_column]] == 'exon', ][[gene_column]], FUN = function(x) { seq_along(x) })
      intron_annotations[intron_annotations[[feature_type_column]] == 'intron', ][['number']]   <- ave(1:nrow(intron_annotations[intron_annotations[[feature_type_column]] == 'intron', ]), intron_annotations[intron_annotations[[feature_type_column]] == 'intron', ][[gene_column]], FUN = function(x) { seq_along(x) })
    }
    # now put into the list
    intron_annotations_per_strand[[strand]] <- intron_annotations
  }
  # merge the intron annotations
  intron_annotations_all <- do.call('rbind', intron_annotations_per_strand)
  # order by position again
  intron_annotations_all <- intron_annotations_all[
    order(
      intron_annotations_all[[gene_column]], 
      intron_annotations_all[[start_column]], 
      intron_annotations_all[[end_column]]), 
  ]
  return(intron_annotations_all)
}

add_gene_body_location <- function(qtl_output, intron_exon_annotation, feature_column='feature_id', variant_id_column='snp_id', variant_position_column='snp_position', variant_chromosome_column='snp_chromosome', variant_classification_column='variant_classification', exon_gene_column='feature', exon_type_column='feature', exon_annotation_column='annotation', exon_feature_column='gene_id', exon_start_column='start', exon_end_column='end', exon_number_column='number') {
  qtl_output <- data.frame(qtl_output)
  intron_exon_annotation <- data.frame(intron_exon_annotation)
  # filter the qtl output to only in-gene
  qtl_output_gene_body <- qtl_output[
    !is.na(qtl_output[[variant_classification_column]]) &
      qtl_output[[variant_classification_column]] == 'in_gene_body', 
  ]
  # subset to only the columns we need
  qtl_output_gene_body <- qtl_output_gene_body[
    , c(feature_column, variant_id_column, variant_position_column)
  ]
  # store a result per gene
  overlap_per_feature <- list()
  # check each gene
  for (feature in unique(qtl_output_gene_body[[feature_column]])) {
    # subset to that feature
    qtl_output_feature <- qtl_output_gene_body[
      !is.na(qtl_output_gene_body[[feature_column]]) &
        qtl_output_gene_body[[feature_column]] == feature, 
    ]
    # extract the positions for that gene from the exon info as well
    intron_exon_annotation_feature <- intron_exon_annotation[
      !is.na(intron_exon_annotation[[exon_feature_column]]) &
        intron_exon_annotation[[exon_feature_column]] == feature, 
    ]
    print(head(qtl_output_feature))
    print(head(intron_exon_annotation_feature))
    # make iranges objects
    ranges_qtl <- IRanges(start = qtl_output_feature[[variant_position_column]], end = qtl_output_feature[[variant_position_column]])
    ranges_exons <- IRanges(start = intron_exon_annotation_feature[[exon_start_column]], end = intron_exon_annotation_feature[[exon_end_column]])
    # get the overlapping entries
    hit_idx <- findOverlaps(ranges_qtl, ranges_exons, select = "first")
    # combine output by for each QTL cbinding the first overlapping row from the exon table
    qtl_feature_to_annotation <- cbind(
      qtl_output_feature,
      intron_exon_annotation_feature[hit_idx, , drop = FALSE]
    )
    # put in the list
    overlap_per_feature[[feature]] <- qtl_feature_to_annotation
  }
  # merge all the overlap tables
  overlaps_all <- do.call('rbind', overlap_per_feature)
  # add this to the qtl output
  qtl_output <- cbind(
    qtl_output, 
    overlaps_all[
      match(
        paste(qtl_output[[feature_column]], qtl_output[[variant_id_column]]), 
        paste(overlaps_all[[feature_column]], overlaps_all[[variant_id_column]])
      ), c(exon_type_column, exon_number_column, exon_annotation_column), drop = F
    ]
  )
  return(qtl_output)
}


add_utr_info <- function(exon_intron_annotation, cds_annotation, feature_id_column='gene_id', start_column='start', end_column='end') {
  # make a list where we have this info added
  exon_intro_with_utr_l <- list()
  # check each feature
  for (feature in unique(exon_intron_annotation[[feature_id_column]])) {
    # subset the exon-intro data
    exon_intron_annotation_feature <- exon_intron_annotation[!is.na(exon_intron_annotation[[feature_id_column]]) & exon_intron_annotation[[feature_id_column]] == feature, ]
    # get the start and stop of the cds annotation
    cds_annotation_feature <- cds_annotation[!is.na(cds_annotation[[feature_id_column]]) & cds_annotation[[feature_id_column]] == feature, ]
    # extract the start of the CDS
    cds_starts <- as.vector(cds_annotation_feature[, start_column])
    # keep only non NA and finite values
    cds_starts <- cds_starts[!is.na(cds_starts) & is.finite(cds_starts)]
    # same for ends
    cds_ends <- as.vector(cds_annotation_feature[, end_column])
    cds_ends <- cds_ends[!is.na(cds_ends) & is.finite(cds_ends)]
    # check if we have cds info
    if (!((length(cds_starts) > 0) & (length(cds_ends) > 0))) {
      exon_intron_annotation_feature[['utr']] <- NA
      # if we do, we can add UTR info
    } else {
      cds_start <- min(cds_starts)
      cds_end <- max(cds_ends)
      # check each intron/exon
      exon_intron_annotation_feature[['utr']] <- apply(exon_intron_annotation_feature, 1, function(x) {
        # extract the intron/exon start
        inex_start <- as.numeric(x[[start_column]])
        inex_stop <- as.numeric(x[[end_column]])
        # check if the stop is before the CDS
        if (inex_stop <  cds_start) {
          return('in_utr')
          # check if the start is after the CDS
        } else if (inex_stop < cds_start) {
          return('in_utr')
          # check if the intron/exon start overlaps the UTR
        } else if (inex_start < cds_start & inex_stop > cds_start) {
          return ('utr_overlap')
          # check if the intron/exon stop overlaps the UTR
        } else if (inex_start < cds_end & inex_stop > cds_end) {
          return ('utr_overlap')
          # otherwise there is no overlap
        } else {
          return ('no_utr')
        }
      })
    }
    # put feature in the list
    exon_intro_with_utr_l[[feature]] <- exon_intron_annotation_feature
  }
  # merge utr info
  exon_intro_with_utr <- do.call('rbind', exon_intro_with_utr_l)
  return(exon_intro_with_utr)
}

add_pct_cds_info <- function(intron_exon_utr_annotation, feature_id_column='gene_id', start_column='start', end_column='end', strand_column='strand', utr_column='utr', cds_bins=8) {
  # make a list where we have this info added
  exon_intro_with_utr_wpct_l <- list()
  # check each feature
  for (feature in unique(intron_exon_utr_annotation[[feature_id_column]])) {
    # extract feature table
    intex_utr_anno_ft <- intron_exon_utr_annotation[!is.na(intron_exon_utr_annotation[[feature_id_column]]) & intron_exon_utr_annotation[[feature_id_column]] == feature, ]
    # get the non-UTR annotations
    intex_utr_anno_ft_nonutr <- intex_utr_anno_ft[!is.na(intex_utr_anno_ft[[utr_column]]) & intex_utr_anno_ft[[utr_column]] == 'no_utr', ]
    # add new column
    intex_utr_anno_ft[['full_anno']] <- NA
    # set the end of the 5' UTR
    utr_5_end <- NULL
    # set the start of the 3' UTR
    utr_3_start <- NULL
    # check the types of strands
    strands <- unique(intex_utr_anno_ft[[strand_column]])
    # check for non NA
    strands <- strands[!is.na(strands)]
    # if the length is exactly 1, then we can proceed
    if (length(strands) == 1) {
      if (strands[[1]] == '+') {
        # get start and end
        utr_5_end <- min(intex_utr_anno_ft_nonutr[[start_column]])
        utr_3_start <- max(intex_utr_anno_ft_nonutr[[end_column]])
        # use that to annotate the UTR
        if (nrow(intex_utr_anno_ft[intex_utr_anno_ft[[end_column]] < utr_5_end, ]) > 0) {
          intex_utr_anno_ft[intex_utr_anno_ft[[end_column]] < utr_5_end, ][['full_anno']] <- '5 prime UTR'
        }
        if (nrow(intex_utr_anno_ft[intex_utr_anno_ft[[start_column]] > utr_3_start, ]) > 0) {
          intex_utr_anno_ft[intex_utr_anno_ft[[start_column]] > utr_3_start, ][['full_anno']] <- '3 prime UTR'
        }
        # calculate fraction where the end of the intron/exon is based on the end of the CDS
        intex_utr_anno_ft_nonutr[['full_anno']] <- ceiling((intex_utr_anno_ft_nonutr[[end_column]] - min(intex_utr_anno_ft_nonutr[[end_column]])) / (utr_3_start - min(intex_utr_anno_ft_nonutr[[end_column]])) * (cds_bins - 1)) + 1
        # add that to the original table
        intex_utr_anno_ft[!is.na(intex_utr_anno_ft[[utr_column]]) & intex_utr_anno_ft[[utr_column]] == 'no_utr',  ][['full_anno']] <- intex_utr_anno_ft_nonutr[['full_anno']]
        
      } else if(strands[[1]] == '-'){
        # get start and end
        utr_5_start <- max(intex_utr_anno_ft_nonutr[[end_column]])
        utr_3_end <- min(intex_utr_anno_ft_nonutr[[start_column]])
        if (nrow(intex_utr_anno_ft[intex_utr_anno_ft[[end_column]] < utr_3_end, ]) > 0) {
          intex_utr_anno_ft[intex_utr_anno_ft[[end_column]] < utr_3_end, ][['full_anno']] <- '3 prime UTR'
        }
        if (nrow(intex_utr_anno_ft[intex_utr_anno_ft[[start_column]] > utr_5_start, ])) {
          intex_utr_anno_ft[intex_utr_anno_ft[[start_column]] > utr_5_start, ][['full_anno']] <- '5 prime UTR'
        }
        # calculate fraction where the end of the intron/exon is based on the end of the CDS
        intex_utr_anno_ft_nonutr[['full_anno']] <- ceiling(
          (1 - (
            (intex_utr_anno_ft_nonutr[[start_column]] - min(intex_utr_anno_ft_nonutr[[start_column]])) / 
              (max(intex_utr_anno_ft_nonutr[[start_column]]) - min(intex_utr_anno_ft_nonutr[[start_column]]))
          )) * (cds_bins - 1)
        ) + 1
        # add that to the original table
        intex_utr_anno_ft[!is.na(intex_utr_anno_ft[[utr_column]]) & intex_utr_anno_ft[[utr_column]] == 'no_utr',  ][['full_anno']] <- intex_utr_anno_ft_nonutr[['full_anno']]
      }
    }
    else {
      intex_utr_anno_ft[['full_anno']] <- NA
    }
    # add back to the list
    exon_intro_with_utr_wpct_l[[feature]] <- intex_utr_anno_ft
  }
  # merge all
  exon_intro_with_utr_wpct <- do.call('rbind', exon_intro_with_utr_wpct_l)
  return(exon_intro_with_utr_wpct)
}

add_pct_trons_info <- function(intron_exon_utr_annotation, feature_id_column='gene_id', start_column='start', end_column='end', strand_column='strand', utr_column='utr', cds_bins=8) {
  # make a list where we have this info added
  exon_intro_with_utr_wpct_l <- list()
  intron_exon_utr_annotation <- data.frame(intron_exon_utr_annotation)
  # check each feature
  for (feature in unique(intron_exon_utr_annotation[[feature_id_column]])) {
    # extract feature table
    intex_utr_anno_ft <- intron_exon_utr_annotation[!is.na(intron_exon_utr_annotation[[feature_id_column]]) & intron_exon_utr_annotation[[feature_id_column]] == feature, ]
    # get the non-UTR annotations
    intex_utr_anno_ft_nonutr <- intex_utr_anno_ft[!is.na(intex_utr_anno_ft[[utr_column]]) & intex_utr_anno_ft[[utr_column]] == 'no_utr', ]
    # add new column
    intex_utr_anno_ft[['full_anno']] <- NA
    # check the types of strands
    strands <- unique(intex_utr_anno_ft[[strand_column]])
    # check for non NA
    strands <- strands[!is.na(strands)]
    if(nrow(intex_utr_anno_ft_nonutr) > 0) {
      # if the length is exactly 1, then we can proceed
      if (length(strands) == 1) {
        if (strands[[1]] == '+') {
          # get start and end
          utr_5_end <- min(intex_utr_anno_ft_nonutr[[start_column]])
          utr_3_start <- max(intex_utr_anno_ft_nonutr[[end_column]])
          # use that to annotate the UTR
          if (nrow(intex_utr_anno_ft[intex_utr_anno_ft[[end_column]] < utr_5_end, ]) > 0) {
            intex_utr_anno_ft[intex_utr_anno_ft[[end_column]] < utr_5_end, ][['full_anno']] <- '5 prime UTR'
          }
          if (nrow(intex_utr_anno_ft[intex_utr_anno_ft[[start_column]] > utr_3_start, ]) > 0) {
            intex_utr_anno_ft[intex_utr_anno_ft[[start_column]] > utr_3_start, ][['full_anno']] <- '3 prime UTR'
          }
          # order by start position
          intex_utr_anno_ft_nonutr <- intex_utr_anno_ft_nonutr[order(intex_utr_anno_ft_nonutr[[start_column]]), ]
          # then add an increasing number
          intex_utr_anno_ft_nonutr[['region_nr']] <- 1:nrow(intex_utr_anno_ft_nonutr)
          # then divide each value by the total
          intex_utr_anno_ft_nonutr[['full_anno']] <- ceiling(intex_utr_anno_ft_nonutr[['region_nr']] / nrow(intex_utr_anno_ft_nonutr) * (cds_bins))
          # add that info back to the original table
          intex_utr_anno_ft[!is.na(intex_utr_anno_ft[[utr_column]]) & intex_utr_anno_ft[[utr_column]] == 'no_utr', ][['full_anno']] <- intex_utr_anno_ft_nonutr[match(intex_utr_anno_ft[!is.na(intex_utr_anno_ft[[utr_column]]) & intex_utr_anno_ft[[utr_column]] == 'no_utr', ][[start_column]], intex_utr_anno_ft_nonutr[[start_column]]), ][['full_anno']]
        } else if(strands[[1]] == '-'){
          # get start and end
          utr_5_start <- max(intex_utr_anno_ft_nonutr[[end_column]])
          utr_3_end <- min(intex_utr_anno_ft_nonutr[[start_column]])
          if (nrow(intex_utr_anno_ft[intex_utr_anno_ft[[end_column]] < utr_3_end, ]) > 0) {
            intex_utr_anno_ft[intex_utr_anno_ft[[end_column]] < utr_3_end, ][['full_anno']] <- '3 prime UTR'
          }
          if (nrow(intex_utr_anno_ft[intex_utr_anno_ft[[start_column]] > utr_5_start, ])) {
            intex_utr_anno_ft[intex_utr_anno_ft[[start_column]] > utr_5_start, ][['full_anno']] <- '5 prime UTR'
          }
          # order by start position
          intex_utr_anno_ft_nonutr <- intex_utr_anno_ft_nonutr[order(intex_utr_anno_ft_nonutr[[start_column]], decreasing = T), ]
          # then add an increasing number
          intex_utr_anno_ft_nonutr[['region_nr']] <- 1:nrow(intex_utr_anno_ft_nonutr)
          # then divide each value by the total
          intex_utr_anno_ft_nonutr[['full_anno']] <- ceiling(intex_utr_anno_ft_nonutr[['region_nr']] / nrow(intex_utr_anno_ft_nonutr) * (cds_bins))
          # add that info back to the original table
          intex_utr_anno_ft[!is.na(intex_utr_anno_ft[[utr_column]]) & intex_utr_anno_ft[[utr_column]] == 'no_utr', ][['full_anno']] <- intex_utr_anno_ft_nonutr[match(intex_utr_anno_ft[!is.na(intex_utr_anno_ft[[utr_column]]) & intex_utr_anno_ft[[utr_column]] == 'no_utr', ][[start_column]], intex_utr_anno_ft_nonutr[[start_column]]), ][['full_anno']]
        }
      }
      else {
        intex_utr_anno_ft[['full_anno']] <- NA
      }
    }
    else {
      intex_utr_anno_ft[['full_anno']] <- NA
    }
    # add back to the list
    exon_intro_with_utr_wpct_l[[feature]] <- intex_utr_anno_ft
  }
  # merge all
  exon_intro_with_utr_wpct <- do.call('rbind', exon_intro_with_utr_wpct_l)
  return(exon_intro_with_utr_wpct)
}


#' Get Color List
#'
#' This function generates a list of colors corresponding to a given vector of names. It ensures each unique name is assigned a unique color.
#'
#' @param vector_of_names A character vector containing the names for which colors are to be generated.
#' @param use_sampling A logical value indicating whether to use sampling when generating colors (default is FALSE).
#' @param color_indices An optional vector of color indices to use for sampling.
#'
#' @return A named list of colors, where each name in the input vector is assigned a unique color.
#'
#' @examples
#' \dontrun{
#' color_list <- get_color_list(c("apple", "banana", "cherry"))
#' }
#'
get_color_list <- function(vector_of_names, use_sampling=F, color_indices=NULL) {
  # get the unique entries
  vector_unique <- unique(vector_of_names)
  # remove any NA
  vector_unique <- vector_unique[!is.na(vector_unique)]
  # get some colors
  colors_to_use <- NULL
  if (length(vector_unique) > 74) {
    colors_to_use <- roycols::sample_tons_of_colors(length(vector_unique), use_sampling = use_sampling, color_indices = color_indices)
  }
  else{
    colors_to_use <- roycols::sample_many_colours(length(vector_unique), use_sampling = use_sampling, color_indices = color_indices)
  }
  # turn into a list
  colors_to_use_list <- as.list(colors_to_use)
  # and add the names
  names(colors_to_use_list) <- vector_unique
  return(colors_to_use_list)
}


####################
# Main code        #
####################

# location of the gene annotations with exons
exons_annotation_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/ncbi_gencode/GCF_000001405.40/genomic.gtf.gz'
# read the annotation file
exons_annotation <- read.table(exons_annotation_loc, header = F, sep = '\t')
# set column names
colnames(exons_annotation) <- c('seqname','source','feature','start','end','score','strand','frame','attributes')
# add the gene ID
exons_annotation[['gene_id']] <- str_extract(exons_annotation$attributes, "(?<=gene_id )[^;]+")
# add biotype
exons_annotation[['gene_biotype']] <- str_extract(exons_annotation$attributes, "(?<=gene_biotype )[^;]+")
# save this info
exons_annotation_exons_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/ncbi_gencode/GCF_000001405.40/gencode_exons.tsv.gz'
write.table(exons_annotation, gzfile(exons_annotation_exons_loc), row.names = F, col.names = T, sep = '\t', quote = T)
# with a checksum
mdfiver::create_sha256_for_file(exons_annotation_exons_loc)
# reload
exons_annotation <- fread(exons_annotation_exons_loc, header = T, sep = '\t')

# subset to the transcripts
exons_annotation_transcripts <- exons_annotation[exons_annotation[['feature']] == 'transcript', ]
# add a TSS
exons_annotation_transcripts[['TSS']] <- exons_annotation_transcripts[['start']]
# but make this the end if we are on the negative strand
exons_annotation_transcripts[exons_annotation_transcripts[['strand']] == '-', ][['TSS']] <- exons_annotation_transcripts[exons_annotation_transcripts[['strand']] == '-', ][['end']]
# save this info
exons_annotatio_transcripts_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/ncbi_gencode/GCF_000001405.40/gencode_with_transcripts.tsv.gz'
write.table(exons_annotation_transcripts, gzfile(exons_annotatio_transcripts_loc), row.names = F, col.names = T, sep = '\t', quote = T)
# with a checksum
mdfiver::create_sha256_for_file(exons_annotatio_transcripts_loc)
# reload
exons_annotation_transcripts <- fread(exons_annotatio_transcripts_loc, header = T, sep = '\t')

# get exon/intron info
exon_intron_annotation <- add_introns_to_exons_info(exons_annotation)
# add the intron/exon number with the actual annotation
exon_intron_annotation[['annotation']] <- paste(exon_intron_annotation[['feature']], exon_intron_annotation[['number']])
# save this info
exons_annotation_anns_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/ncbi_gencode/GCF_000001405.40/gencode_with_annotations.tsv.gz'
write.table(exon_intron_utr_annotation, gzfile(exons_annotation_anns_loc), row.names = F, col.names = T, sep = '\t', quote = T)
# with a checksum
mdfiver::create_sha256_for_file(exons_annotation_anns_loc)
# reload
exon_intron_annotation <- fread(exons_annotation_anns_loc)

# get the CDS info as well
exons_annotation_cds <- exons_annotation[exons_annotation[['feature']] == 'CDS', ]
# add pct info
exon_intron_utr_annotation_pct <- add_pct_cds_info(exon_intron_utr_annotation)
# save this info
exons_annotatio_utr_pct_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/ncbi_gencode/GCF_000001405.40/gencode_with_annotations_bins.tsv.gz'
write.table(exon_intron_utr_annotation_pct, gzfile(exons_annotatio_utr_pct_loc), row.names = F, col.names = T, sep = '\t', quote = F)
mdfiver::create_sha256_for_file(exons_annotatio_utr_pct_loc)
# reload
exon_intron_utr_annotation_pct <- fread(exons_annotatio_utr_pct_loc, header = T, sep = '\t')
# also with smaller bins
exon_intron_utr_annotation_pct_5 <- add_pct_cds_info(exon_intron_utr_annotation, cds_bins = 5)

# add pct info
exon_intron_utr_annotation_pcttrons <- add_pct_trons_info(exon_intron_utr_annotation)
# save this info
exons_annotatio_utr_pcttrons_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/ncbi_gencode/GCF_000001405.40/gencode_with_annotations_bins_pcttrons.tsv.gz'
write.table(exon_intron_utr_annotation_pcttrons, gzfile(exons_annotatio_utr_pcttrons_loc), row.names = F, col.names = T, sep = '\t', quote = F)
mdfiver::create_sha256_for_file(exons_annotatio_utr_pcttrons_loc)
# reload this info
exon_intron_utr_annotation_pcttrons <- fread(exons_annotatio_utr_pcttrons_loc, header = T, sep = '\t')
# also with smaller bins
exon_intron_utr_annotation_pcttrons_5 <- add_pct_trons_info(exon_intron_utr_annotation, cds_bins = 5)

# get the location of the file annotating the QTL variants and their SCREEN annotation
variant_to_screen_region_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_screen_overlap.tsv.gz'
# read
variant_to_screen_region <- fread(variant_to_screen_region_loc, header = T, sep = '\t')
# get the annotation for each region
screen_region_anno_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/encode_cres/v4/GRCh38-cCREs.bed'
# read that annotation
screen_region_anno <- fread(screen_region_anno_loc, sep = '\t', header = F)
# set column names
colnames(screen_region_anno) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'function')
# add identifier as we have them in other data
screen_region_anno[['feature']] <- paste(screen_region_anno[['chromosome']], screen_region_anno[['start']], screen_region_anno[['end']], sep = '-')
# make a table which is the variant to screen region
variant_to_screen_anno <- cbind(
  variant_to_screen_region, 
  screen_region_anno[match(variant_to_screen_region[['overlapping_feature']], screen_region_anno[['feature']]), c('function'), drop = F]
)

# location of the eQTL output
qtl_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/'
# get all the QTL output
qtl_output <- get_qtls_per_celltype_limix(qtl_output_loc)
# add the top effect information
qtl_output <- add_top_effect_annotation(qtl_output)
# merge all the results
qtl_output_all <- do.call('rbind', qtl_output)
# add 'chr' to the chromosome
qtl_output_all[['snp_chromosome']] <- paste0('chr', qtl_output_all[['snp_chromosome']])

# get the cpeaks overlaps for each variant
qtl_variants_all_cpeaks_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz'
qtl_variants_all_cpeaks <- fread(qtl_variants_all_cpeaks_loc, header = T, sep = '\t')
# add overlapping feature to QTL
qtl_output_all[['region']] <- qtl_variants_all_cpeaks[match(qtl_output_all[['snp_id']], qtl_variants_all_cpeaks[['snp_id']]), ][['overlapping_feature']]
# filter on significance
qtl_output_all_sig <- qtl_output_all[qtl_output_all[['feature_q_value']] < 0.05 &
                                       qtl_output_all[['p_value']] < qtl_output_all[['pval_nominal_threshold_global']], ]

# location of the i-eqtl output
iqtl_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/output/ut_and_24hca_significant/L1/'
# get all the QTL output
iqtl_output <- get_qtls_per_celltype_limix(iqtl_output_loc, output_file = 'inflammation_final/iqtl_results_all_eigenmt_qval.tsv.gz', gene_column='feature', significance_column='feature_q_value', significance_cutoff=0.05, nominal_cutoff_column=NULL, nominal_significance_column=NULL)
# merge all the results
iqtl_output_all <- do.call('rbind', iqtl_output)
# filter on signifiacnce
iqtl_output_all_sig <- iqtl_output_all[iqtl_output_all[['feature_q_value']] < 0.05 &
                                         iqtl_output_all[['feature_bf_eigen']] < 0.05, ]
# rename the columns
colnames(iqtl_output_all_sig) <- c('i_beta', 'i_beta_se', 'i_empirical_feature_p_value', 'i_p_value', 'snp_id', 'feature_id', 'i_n_tests_feature', 'i_feature_bf_eigen', 'i_total_bf_eigen', 'i_feature_q_value', 'i_cell_type')
# merge the iqtl to the qtl output
qtl_output_all_sig_wi <- merge(qtl_output_all_sig, iqtl_output_all_sig, on = c('snp_id', 'feature_id', 'cell_type'), all.x = T)

# get the distances
qtl_distances <- get_closest_flanks(qtl_output_all_sig_wi, 'feature_start', 'feature_end', 'snp_position', 'snp_position')
# add those distances
qtl_output_all_sig_wi[['distance']] <- qtl_distances[['min_dist']]
# add column denoting the interaction effect
qtl_output_all_sig_wi[['interaction_direction']] <- 'none'
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['i_beta']]) & sign(qtl_output_all_sig_wi[['i_beta']]) == 1, ][['interaction_direction']] <- 'positive'
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['i_beta']]) & sign(qtl_output_all_sig_wi[['i_beta']]) == -1, ][['interaction_direction']] <- 'negative'
# add the end to start distance
qtl_output_all_sig_wi[['feature_end_start_dist']] <- qtl_output_all_sig_wi[['feature_end']] - qtl_output_all_sig_wi[['feature_start']]

# add to the QTL info the strand
qtl_output_all_sig_wi[['strand']] <- exons_annotation_transcripts[match(qtl_output_all_sig_wi[['feature_id']], exons_annotation_transcripts[['gene_id']]), ][['strand']]
# add to the QTL info the TSS
qtl_output_all_sig_wi[['TSS']] <- exons_annotation_transcripts[match(qtl_output_all_sig_wi[['feature_id']], exons_annotation_transcripts[['gene_id']]), ][['TSS']]
# subset to the genes for the exon annotations
exons_annotation_genes <- exons_annotation[exons_annotation[['feature']] == 'gene', ]
# add the biotype
qtl_output_all_sig_wi[['gene_biotype']] <- exons_annotation_genes[match(qtl_output_all_sig_wi[['feature_id']], exons_annotation_genes[['gene_id']]), ][['gene_biotype']]

# add the location of the variant
qtl_output_all_sig_wi <- classify_variant_location(qtl_output_all_sig_wi)
# now add a corrected distance
qtl_output_all_sig_wi[['distance_directional']] <- NA
# based on the distance already calculated
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['variant_classification']]), 'distance_directional'] <- qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['variant_classification']]), 'distance']
# but flipping direction for when in front of gene
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['variant_classification']]) & qtl_output_all_sig_wi[['variant_classification']] == 'before_tss', 'distance_directional'] <- -1 * qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['variant_classification']]) & qtl_output_all_sig_wi[['variant_classification']] == 'before_tss', 'distance_directional']
# add the gen-body size
qtl_output_all_sig_wi[['gb_size']] <- abs(qtl_output_all_sig_wi[['feature_start']] - qtl_output_all_sig_wi[['feature_end']])

# calculate distance to tss
qtl_output_all_sig_wi[['tss_dist']] <- qtl_output_all_sig_wi[['TSS']] - qtl_output_all_sig_wi[['snp_position']]
# where the it was on the negative strand, the distance is in the other direction
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['strand']]) & qtl_output_all_sig_wi[['strand']] == -1, ][['tss_dist']] <- -1 * qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['strand']]) & qtl_output_all_sig_wi[['strand']] == -1, ][['tss_dist']]
# add info on if the beta is the same direction as the i_beta
qtl_output_all_sig_wi[['g_i_beta_same']] <- !is.na(qtl_output_all_sig_wi[['beta']]) & !is.na(qtl_output_all_sig_wi[['i_beta']]) & sign(qtl_output_all_sig_wi[['beta']]) == qtl_output_all_sig_wi[['i_beta']]
# make a classification column
qtl_output_all_sig_wi[['interaction_beta_direction']] <- 'none'
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['i_beta']]) & 
                        sign(qtl_output_all_sig_wi[['i_beta']]) == 1 &
                        sign(qtl_output_all_sig_wi[['beta']]) == 1, 'interaction_beta_direction'] <- '+/+'
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['i_beta']]) & 
                        sign(qtl_output_all_sig_wi[['i_beta']]) == -1 &
                        sign(qtl_output_all_sig_wi[['beta']]) == 1, 'interaction_beta_direction'] <- '+/-'
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['i_beta']]) & 
                        sign(qtl_output_all_sig_wi[['i_beta']]) == 1 &
                        sign(qtl_output_all_sig_wi[['beta']]) == -1, 'interaction_beta_direction'] <- '-/+'
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['i_beta']]) & 
                        sign(qtl_output_all_sig_wi[['i_beta']]) == -1 &
                        sign(qtl_output_all_sig_wi[['beta']]) == -1, 'interaction_beta_direction'] <- '-/-'

# check if there is a difference
kruskal.test(tss_dist ~ interaction_direction, data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['tss_dist']]), ])
# for protein coding specifically
kruskal.test(tss_dist ~ interaction_direction, data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['tss_dist']]) & !is.na(qtl_output_all_sig_wi[['gene_biotype']]) & qtl_output_all_sig_wi[['gene_biotype']] == 'protein_coding', ])

# keep top over all cell types
qtl_output_all_sig_wi_wi <- qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['i_beta']]), ]
qtl_output_all_sig_wi_wi_tc[['i_zscore']] <- qtl_output_all_sig_wi_wi_tc[['i_beta']] / qtl_output_all_sig_wi_wi_tc[['i_beta_se']]
qtl_output_all_sig_wi_wi_tc[['zscore']] <- qtl_output_all_sig_wi_wi_tc[['beta']] / qtl_output_all_sig_wi_wi_tc[['beta_se']]
#qtl_output_all_sig_wi_wi_tc <- qtl_output_all_sig_wi[order(qtl_output_all_sig_wi_wi_tc[['i_zscore']]), ]
qtl_output_all_sig_wi_wi_tc <- qtl_output_all_sig_wi[order(qtl_output_all_sig_wi[['p_value']]), ]
qtl_output_all_sig_wi_wi_tc <- qtl_output_all_sig_wi[order(qtl_output_all_sig_wi[['i_feature_bf_eigen']]), ]
qtl_output_all_sig_wi_wi_tc <- qtl_output_all_sig_wi_wi_tc[!duplicated(qtl_output_all_sig_wi[['feature_id']]), ]

# check if there is a difference
kruskal.test(distance_directional ~ interaction_direction, data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['distance_directional']]), ])
# p-value < 2.2e-16
# do post-hoc test for positive vs negative
wilcox.test(distance_directional ~ interaction_direction, data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['distance_directional']]) & qtl_output_all_sig_wi[['interaction_direction']] != 'none', ])
# p-value = 0.0059

# check if the MAF values are significantly different as well
wilcox.test(
  x = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['distance_directional']]) & sign(qtl_output_all_sig_wi[['distance_directional']]) == -1, ][['maf']], 
  y = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['distance_directional']]) & sign(qtl_output_all_sig_wi[['distance_directional']]) == 1, ][['maf']]
)
# p-value 0.06335

# check if there is a difference for protein coding or not
kruskal.test(distance_directional ~ gene_biotype, data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['distance_directional']]), ])
# p-value = 0.05437
# check if the maf is related to whether or not it is protein coding
kruskal.test(maf ~ gene_biotype, data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['distance_directional']]), ])
# p-value = 0.01058
# next, try to model multiple variable at once by first checking if the distance at all plays a role
qtl_output_all_sig_wi_factorized <- qtl_output_all_sig_wi
summary(
  lm(
    formula = as.formula('distance_directional ~ 
                          interaction_direction + 
                          gb_size +
                          maf +
                          beta + 
                          gene_biotype'), 
    data = qtl_output_all_sig_wi_factorized[qtl_output_all_sig_wi_factorized[['is_top_variant']] & !is.na(qtl_output_all_sig_wi_factorized[['distance_directional']]), ]
  )
)
# when modelling all of these, only yes vs no interaction is significant
# interaction_directionnone           3.342e+04  8.928e+03   3.744 0.000187 ***
# interaction_directionpositive       1.915e+04  1.197e+04   1.600 0.109849
# exclusing the 'none' category
summary(
  lm(
    formula = as.formula('distance_directional ~ 
                          interaction_direction + 
                          gb_size +
                          maf +
                          beta + 
                          gene_biotype'), 
    data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['distance_directional']]) & qtl_output_all_sig_wi[['interaction_direction']] != 'none', ]
  )
)
# not significant at the numerical level
#                                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                        -5.401e+04  1.862e+04  -2.901 0.004037 ** 
#   interaction_directionpositive       1.463e+03  8.829e+03   0.166 0.868535    
# gb_size                            -1.659e-01  6.214e-02  -2.669 0.008080 ** 
#   maf                                 4.576e+04  3.183e+04   1.438 0.151675    
# beta                                2.287e+04  5.928e+03   3.857 0.000145 ***
#   gene_biotypeprotein_coding          3.929e+04  1.327e+04   2.960 0.003356 ** 
#   gene_biotypetranscribed_pseudogene  1.905e+05  4.793e+04   3.975 9.11e-05 ***
#   gene_biotypeV_segment               6.577e+03  6.656e+04   0.099 0.921365    
# now run with the distance as a categorical value
qtl_output_all_sig_wi_factorized[['distance_directional_category']] <- as.factor(sign(qtl_output_all_sig_wi_factorized[['distance_directional']]))
summary(
  glm(
    formula = as.formula('distance_directional_category ~ 
                          interaction_direction + 
                          gb_size + 
                          maf +
                          beta + 
                          gene_biotype'), 
    data = qtl_output_all_sig_wi_factorized[qtl_output_all_sig_wi_factorized[['is_top_variant']] & !is.na(qtl_output_all_sig_wi_factorized[['distance_directional_category']]), ], 
    family = 'binomial'
  )
)
# shows a significant effect of positive/negative
#                                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                        -2.748e+00  1.453e+00  -1.891  0.05866 .  
# interaction_directionnone           2.630e+00  2.587e-01  10.166  < 2e-16 ***
#   interaction_directionpositive       8.899e-01  3.149e-01   2.826  0.00471 ** 
#   gb_size                             1.069e-06  3.971e-07   2.691  0.00712 ** 
#   maf                                 1.034e+00  4.169e-01   2.480  0.01315 *  
#   beta                                1.895e-01  7.263e-02   2.609  0.00909 ** 
#   gene_biotypelncRNA                  3.460e-01  1.428e+00   0.242  0.80856    
# gene_biotypencRNA_pseudogene       -1.502e+01  8.827e+02  -0.017  0.98643    
# gene_biotypeprotein_coding          6.869e-01  1.422e+00   0.483  0.62911    
# gene_biotypetranscribed_pseudogene  1.863e+00  1.628e+00   1.145  0.25224    
# gene_biotypeV_segment              -1.444e+01  4.560e+02  -0.032  0.97474    
# gene_biotypeV_segment_pseudogene   -1.471e+01  8.827e+02  -0.017  0.98670   

# this also has UTR info
exon_intron_utr_annotation <- exon_intron_annotation

# now add the intron/exon info
qtl_output_all_sig_wi_gb <- add_gene_body_location(qtl_output_all_sig_wi, exon_intron_annotation)
# add to the ieqtl data
qtl_output_all_sig_wi_gb[['screen_overlap_func']] <- variant_to_screen_anno[match(qtl_output_all_sig_wi_gb[['snp_id']], variant_to_screen_anno[['snp_id']]), ][['function']]

# save QTL output info
ieqtl_annotated_loc <- '~/tables/mo_ieqtl_annotated.rds'
saveRDS(qtl_output_all_sig_wi_gb, ieqtl_annotated_loc)
mdfiver::create_sha256_for_file(ieqtl_annotated_loc)
# read back
qtl_output_all_sig_wi_gb <- readRDS(ieqtl_annotated_loc)

# add UTR info as well
qtl_output_all_sig_wi_gb[['utr']] <- exon_intron_utr_annotation[match(paste(qtl_output_all_sig_wi_gb[['feature']], qtl_output_all_sig_wi_gb[['number']]), paste(exon_intron_utr_annotation[['feature']], exon_intron_utr_annotation[['number']])), ][['utr']]
qtl_output_all_sig_wi_gb[['full_anno']] <- exon_intron_utr_annotation_pct[match(paste(qtl_output_all_sig_wi_gb[['feature']], qtl_output_all_sig_wi_gb[['number']]), paste(exon_intron_utr_annotation_pct[['feature']], exon_intron_utr_annotation_pct[['number']])), ][['full_anno']]
qtl_output_all_sig_wi_gb[['tron_bin']] <- exon_intron_utr_annotation_pcttrons[match(paste(qtl_output_all_sig_wi_gb[['feature']], qtl_output_all_sig_wi_gb[['number']]), paste(exon_intron_utr_annotation_pcttrons[['feature']], exon_intron_utr_annotation_pcttrons[['number']])), ][['full_anno']]

# just calculate the distance to the gene body


# plot the distances per group
p_snp_gene_interaction_distances <- ggplot(data = qtl_output_all_sig_wi, mapping = aes(x = distance, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and gene') + 
  ylab('Density') + 
  ggtitle('Distance between variant and gene\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_gene_interaction_distances
# the same, but only for the top variant per gene
p_snp_gene_interaction_distances_top <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = distance, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and gene') + 
  ylab('Density') + 
  ggtitle('Distance between top variant and gene\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_gene_interaction_distances_top
# now without variant in genes
p_snp_gene_interaction_distances_distnozero <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['distance']] > 0, ], mapping = aes(x = distance, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and gene\n(excluding variants in genes)') + 
  ylab('Density') + 
  ggtitle('Distance between variant and gene\nper interaction category (dist>0)') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_gene_interaction_distances_distnozero
# and top variants outside of genes
p_snp_gene_interaction_distances_distnozero_top <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['distance']] > 0 & qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = distance, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and gene\n(excluding variants in genes)') + 
  ylab('Density') + 
  ggtitle('Distance between top variant and gene\nper interaction category (dist>0)') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_gene_interaction_distances_distnozero_top
# show all together
p_snp_gene_interaction_all <- plot_grid(
  p_snp_gene_interaction_distances, 
  p_snp_gene_interaction_distances_top, 
  p_snp_gene_interaction_distances_distnozero, 
  p_snp_gene_interaction_distances_distnozero_top
)
p_snp_gene_interaction_all
# and top variants outside of genes, but also mono
p_snp_gene_interaction_distances_distnozero_top_mono <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['distance']] > 0 & qtl_output_all_sig_wi[['is_top_variant']] & qtl_output_all_sig_wi[['cell_type']] == 'monocyte', ], mapping = aes(x = distance, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and gene for\n(excluding variants in genes)') + 
  ylab('Density') + 
  ggtitle('Distance between top variant and gene\nper interaction category (dist>0)\nfor monocytes only') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_gene_interaction_distances_distnozero_top_mono


# plot the distances per group
p_snp_tss_interaction_distances <- ggplot(data = qtl_output_all_sig_wi, mapping = aes(x = tss_dist, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and TSS') + 
  ylab('Density') + 
  ggtitle('Distance between variant and TSS\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_tss_interaction_distances
# plot the distances per group for top effects
p_snp_tss_interaction_distances_top <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = tss_dist, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between top variant and TSS') + 
  ylab('Density') + 
  ggtitle('Distance between top variant and TSS\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_tss_interaction_distances_top
# plot the distances per group for top effects using absolute distances
p_snp_tss_interaction_distances_top_abs <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = abs(tss_dist), fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Absolute istance between top variant and TSS') + 
  ylab('Density') + 
  ggtitle('Absolute distance between top variant and TSS\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_tss_interaction_distances_top_abs

# show all plots
plot_grid(
  p_snp_tss_interaction_distances, 
  p_snp_tss_interaction_distances_top, 
  p_snp_tss_interaction_distances_top_abs, 
  nrow = 2, 
  ncol = 2
)
# save 
ggsave('~/plots/mo_top_eqtl_distances_interactions.pdf', plot = p_snp_tss_interaction_distances_top, width = 5, height = 5)

p_snp_tss_beta_distances <- ggplot(data = qtl_output_all_sig_wi, mapping = aes(x = tss_dist, fill = interaction_beta_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and TSS') + 
  ylab('Density') + 
  ggtitle('Distance between variant and TSS\nper beta category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', '+/+' = '#386CB0', '-/-' = '#F0027F', '+/-' = 'lightblue', '+/-' = 'lightpink')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_tss_beta_distances_top <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = tss_dist, fill = interaction_beta_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between top variant and TSS') + 
  ylab('Density') + 
  ggtitle('Distance between top variant and TSS\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', '+/+' = '#386CB0', '-/-' = '#F0027F', '+/-' = 'lightblue', '+/-' = 'lightpink')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))

# now with directed interaction directions
p_snp_gene_interaction_distances_directional <- ggplot(data = qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['distance_directional']]), ], mapping = aes(x = distance_directional, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and gene') + 
  ylab('Density') + 
  ggtitle('Distance between variant and gene\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# top as well
p_snp_gene_interaction_distances_directional_top <- ggplot(data = qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['distance_directional']]) & qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = distance_directional, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and gene') + 
  ylab('Density') + 
  ggtitle('Distance between variant and gene\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# and show the MAF
p_snp_maf_directional_top <- ggplot(data = qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['distance_directional']]) & qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = maf, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('MAF of variants') + 
  ylab('Density') + 
  ggtitle('MAF of variants\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# per biotype as well
p_snp_maf_biotype_top <- ggplot(data = qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['gene_biotype']]) & qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = maf, fill = gene_biotype)) +
  geom_density(alpha = 0.5) +
  xlab('MAF of variants') + 
  ylab('Density') + 
  ggtitle('MAF of variants\nper biotype category') + 
  labs(fill = 'bioptype') +
  scale_fill_manual(values = roycols::get_color_list(unique(qtl_output_all_sig_wi[['gene_biotype']]))) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# and type
p_snp_maf_qltdirection_top <- ggplot(data = qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['variant_classification']]) & qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = maf, fill = variant_classification)) +
  geom_density(alpha = 0.5) +
  xlab('MAF of variants') + 
  ylab('Density') + 
  ggtitle('MAF of variants\nper variant-gene direction') + 
  labs(fill = 'variant\ndirection') +
  scale_fill_manual(values = roycols::get_color_list(unique(qtl_output_all_sig_wi[['variant_classification']]))) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))


plot_grid(
  # p_snp_gene_interaction_distances_directional_top, 
  p_snp_maf_directional_top, 
  p_snp_maf_biotype_top, 
  p_snp_maf_qltdirection_top
)

# check how the variants are located
variant_gb_numbers <- data.frame(table(qtl_output_all_sig_wi_gb[, c('feature', 'number', 'interaction_direction')]))
# order
variant_gb_numbers <- variant_gb_numbers[order(variant_gb_numbers[['Freq']], decreasing = T), ]
# rename columns
colnames(variant_gb_numbers) <- c('intron_exon', 'intron_exon_nr', 'interaction', 'n_obs')
# remove zeroes
variant_gb_numbers <- variant_gb_numbers[variant_gb_numbers[['n_obs']] > 0, ]
# add intron exon as explicit column
variant_gb_numbers[['intron_exon_with_number']] <- paste(variant_gb_numbers[['intron_exon']], variant_gb_numbers[['intron_exon_nr']])
# add the proportion as well
variant_gb_numbers[['proportion_obs']] <- variant_gb_numbers[['n_obs']] / sum(variant_gb_numbers[['n_obs']])
# order by number, then intro exon
variant_gb_numbers <- variant_gb_numbers[order(variant_gb_numbers[['intron_exon_nr']], variant_gb_numbers[['intron_exon']]), ]
# then make the readable intron/exon and number a factor so that is the order
variant_gb_numbers[['intron_exon_with_number']] <- factor(variant_gb_numbers[['intron_exon_with_number']], levels = variant_gb_numbers[['intron_exon_with_number']])
# make into a plot
p_nrs <- ggplot(
  data = variant_gb_numbers, 
  mapping = aes(x = intron_exon_with_number, y = proportion_obs, fill = intron_exon_with_number) 
) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = roycols::get_color_list(unique(variant_gb_numbers[['intron_exon_with_number']]))) +
  xlab('Intron/Exon') + 
  ylab('Proportion of variants') +
  ggtitle('Location of variants') + 
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# save the plot
ggsave('~/plots/mo_variant_position_eqtls.pdf', width = 14, height = 7, plot = p_nrs)

# make this into a table
intron_exon_nrs <- data.table(table(qtl_output_all_sig_wi_gb[!is.na(qtl_output_all_sig_wi_gb$full_anno), ]$full_anno))
colnames(intron_exon_nrs) <- c('loc', 'n')
# the loc into an factor
intron_exon_nrs[['loc']] <- factor(intron_exon_nrs[['loc']], levels = c('5 prime UTR', 0:8, '3 prime UTR'))
# make that into a plot
p_locs <- ggplot(data = intron_exon_nrs, mapping = aes(x = loc, y = n)) +
  geom_bar(stat = 'identity') +
  xlab('Variant location') + 
  ylab('N') +
  ggtitle('Location of variants') + 
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  # label sizes
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20)
  )
# finally also check exon/intron
prime5_nrs <- data.frame(table(qtl_output_all_sig_wi_gb[!is.na(qtl_output_all_sig_wi_gb[['full_anno']]) & qtl_output_all_sig_wi_gb[['full_anno']] == '5 prime UTR', 'feature']))
colnames(prime5_nrs) <- c('feature', 'nr')
prime5_nrs[['frac']] <- prime5_nrs[['nr']] / sum(prime5_nrs[['nr']])
# add a zero because there is apparently no introns
prime5_nrs <- rbind(prime5_nrs, data.frame('feature' = c('intron'), 'nr' = c(0), 'frac' = c(0)))
prime3_nrs <- data.frame(table(qtl_output_all_sig_wi_gb[!is.na(qtl_output_all_sig_wi_gb[['full_anno']]) & qtl_output_all_sig_wi_gb[['full_anno']] == '3 prime UTR', 'feature']))
colnames(prime3_nrs) <- c('feature', 'nr')
prime3_nrs[['frac']] <- prime3_nrs[['nr']] / sum(prime3_nrs[['nr']])
prime3_nrs <- rbind(prime3_nrs, data.frame('feature' = c('intron'), 'nr' = c(0), 'frac' = c(0)))
no_utr_nrs <- data.frame(table(qtl_output_all_sig_wi_gb[!is.na(qtl_output_all_sig_wi_gb[['full_anno']]) & qtl_output_all_sig_wi_gb[['full_anno']] != '5 prime UTR' & qtl_output_all_sig_wi_gb[['full_anno']] != '3 prime UTR', 'feature']))
colnames(no_utr_nrs) <- c('feature', 'nr')
no_utr_nrs[['frac']] <- no_utr_nrs[['nr']] / sum(no_utr_nrs[['nr']])
# into plots
prime5_p <- ggplot(data = prime5_nrs, mapping = aes(x = feature, y = frac, fill = feature)) +
  geom_bar(stat = 'identity') +
  xlab('Variant location') + 
  ylab('Fraction of variants') +
  ggtitle('Location of variants') + 
  theme(legend.position = 'none') + 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  # label sizes
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20)
  ) +
  scale_fill_manual(values = list('exon' = 'darkblue', 'intron' = 'darkred'))
prime3_p <- ggplot(data = prime3_nrs, mapping = aes(x = feature, y = frac, fill = feature)) +
  geom_bar(stat = 'identity') +
  xlab('Variant location') + 
  ylab('Fraction of variants') +
  ggtitle('Location of variants') + 
  theme(legend.position = 'none') + 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  # label sizes
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20)
  ) +
  scale_fill_manual(values = list('exon' = 'darkblue', 'intron' = 'darkred'))
noutr_p <- ggplot(data = no_utr_nrs, mapping = aes(x = feature, y = frac, fill = feature)) +
  geom_bar(stat = 'identity') +
  xlab('Variant location') + 
  ylab('Fraction of variants') +
  ggtitle('Location of variants') + 
  theme(legend.position = 'none') + 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  # label sizes
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20)
  ) +
  scale_fill_manual(values = list('exon' = 'darkblue', 'intron' = 'darkred'))
# finally put in a full plot
p_locs_combined <- plot_grid(
  p_locs + xlab(''), 
  plot_grid(
    prime5_p, 
    noutr_p + ylab(''), 
    prime3_p + ylab(''), 
    nrow = 1, 
    ncol = 3
  ),
  nrow = 2, 
  ncol = 1, 
  rel_heights = c(2,1)
)
ggsave('~/plots/mo_ievariants_locs.pdf', p_locs_combined, width = 9, height = 6)


# do the same, but now with the intron/exons numbered
intron_exon_nrs <- data.table(table(qtl_output_all_sig_wi_gb[!is.na(qtl_output_all_sig_wi_gb$tron_bin), ]$tron_bin))
colnames(intron_exon_nrs) <- c('loc', 'n')
# the loc into an factor
intron_exon_nrs[['loc']] <- factor(intron_exon_nrs[['loc']], levels = c('5 prime UTR', 0:8, '3 prime UTR'))
# make that into a plot
p_locs <- ggplot(data = intron_exon_nrs, mapping = aes(x = loc, y = n)) +
  geom_bar(stat = 'identity') +
  xlab('Variant location') + 
  ylab('N') +
  ggtitle('Location of variants') + 
  theme(legend.position = 'none') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  # label sizes
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20)
  )
# the the exon/intron fractions
prime5_nrs <- data.frame(table(qtl_output_all_sig_wi_gb[!is.na(qtl_output_all_sig_wi_gb[['tron_bin']]) & qtl_output_all_sig_wi_gb[['tron_bin']] == '5 prime UTR', 'feature']))
colnames(prime5_nrs) <- c('feature', 'nr')
prime5_nrs[['frac']] <- prime5_nrs[['nr']] / sum(prime5_nrs[['nr']])
# add a zero because there is apparently no introns
prime5_nrs <- rbind(prime5_nrs, data.frame('feature' = c('intron'), 'nr' = c(0), 'frac' = c(0)))
prime3_nrs <- data.frame(table(qtl_output_all_sig_wi_gb[!is.na(qtl_output_all_sig_wi_gb[['tron_bin']]) & qtl_output_all_sig_wi_gb[['tron_bin']] == '3 prime UTR', 'feature']))
colnames(prime3_nrs) <- c('feature', 'nr')
prime3_nrs[['frac']] <- prime3_nrs[['nr']] / sum(prime3_nrs[['nr']])
prime3_nrs <- rbind(prime3_nrs, data.frame('feature' = c('intron'), 'nr' = c(0), 'frac' = c(0)))
no_utr_nrs <- data.frame(table(qtl_output_all_sig_wi_gb[!is.na(qtl_output_all_sig_wi_gb[['tron_bin']]) & qtl_output_all_sig_wi_gb[['tron_bin']] != '5 prime UTR' & qtl_output_all_sig_wi_gb[['tron_bin']] != '3 prime UTR', 'feature']))
colnames(no_utr_nrs) <- c('feature', 'nr')
no_utr_nrs[['frac']] <- no_utr_nrs[['nr']] / sum(no_utr_nrs[['nr']])
# into plots
prime5_p <- ggplot(data = prime5_nrs, mapping = aes(x = feature, y = frac, fill = feature)) +
  geom_bar(stat = 'identity') +
  xlab('Variant location') + 
  ylab('Fraction of variants') +
  ggtitle('Location of variants') + 
  theme(legend.position = 'none') + 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  # label sizes
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20)
  ) +
  scale_fill_manual(values = list('exon' = 'darkblue', 'intron' = 'darkred'))
prime3_p <- ggplot(data = prime3_nrs, mapping = aes(x = feature, y = frac, fill = feature)) +
  geom_bar(stat = 'identity') +
  xlab('Variant location') + 
  ylab('Fraction of variants') +
  ggtitle('Location of variants') + 
  theme(legend.position = 'none') + 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  # label sizes
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20)
  ) +
  scale_fill_manual(values = list('exon' = 'darkblue', 'intron' = 'darkred'))
noutr_p <- ggplot(data = no_utr_nrs, mapping = aes(x = feature, y = frac, fill = feature)) +
  geom_bar(stat = 'identity') +
  xlab('Variant location') + 
  ylab('Fraction of variants') +
  ggtitle('Location of variants') + 
  theme(legend.position = 'none') + 
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  # label sizes
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    plot.title = element_text(size = 20)
  ) +
  scale_fill_manual(values = list('exon' = 'darkblue', 'intron' = 'darkred'))
# finally put in a full plot
p_locs_combined <- plot_grid(
  p_locs + xlab(''), 
  plot_grid(
    prime5_p, 
    noutr_p + ylab(''), 
    prime3_p + ylab(''), 
    nrow = 1, 
    ncol = 3
  ),
  nrow = 2, 
  ncol = 1, 
  rel_heights = c(2,1)
)
ggsave('~/plots/mo_ievariants_locs.pdf', p_locs_combined, width = 9, height = 6)

# now instead plot the distances
p_distance_wiedirection <- ggplot(data = qtl_output_all_sig_wi_wi_tc, mapping = aes(y = interaction_direction, x = distance, fill = interaction_direction)) +
  geom_boxplot(outlier.shape = NA) + 
  # and add jitter
  geom_jitter(size = 0.5, alpha = 0.5) +
  # and labels
  xlab('interaction') + 
  ylab('distance to gene') +
  ggtitle('Gene distance per interaction direction') + 
  theme(axis.text.x=element_blank(), 
                 axis.ticks = element_blank()) + 
  theme(legend.position = 'none') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = get_color_list(qtl_output_all_sig_wi_gb[['interaction_direction']]))
# with the forward and back
p_distance_wiedirection_both <- ggplot(data = qtl_output_all_sig_wi_wi_tc, mapping = aes(y = interaction_direction, x = distance_directional, fill = interaction_direction)) +
  geom_boxplot(outlier.shape = NA) + 
  # and add jitter
  geom_jitter(size = 0.5, alpha = 0.5) +
  # and labels
  xlab('interaction') + 
  ylab('distance to gene') +
  ggtitle('Gene distance per interaction direction') + 
  theme(axis.text.x=element_blank(), 
        axis.ticks = element_blank()) + 
  theme(legend.position = 'none') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = get_color_list(qtl_output_all_sig_wi_gb[['interaction_direction']]))
