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

####################
# Main code        #
####################

# location of the eQTL output
qtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/'
# get all the QTL output
qtl_output <- get_qtls_per_celltype_limix(qtl_output_loc)
# add the top effect information
qtl_output <- add_top_effect_annotation(qtl_output)
# merge all the results
qtl_output_all <- do.call('rbind', qtl_output)
# add 'chr' to the chromosome
qtl_output_all[['snp_chromosome']] <- paste0('chr', qtl_output_all[['snp_chromosome']])

# get the cpeaks overlaps for each variant
qtl_variants_all_cpeaks_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz'
qtl_variants_all_cpeaks <- fread(qtl_variants_all_cpeaks_loc, header = T, sep = '\t')
# add overlapping feature to QTL
qtl_output_all[['region']] <- qtl_variants_all_cpeaks[match(qtl_output_all[['snp_id']], qtl_variants_all_cpeaks[['snp_id']]), ][['overlapping_feature']]
# filter on significance
qtl_output_all_sig <- qtl_output_all[qtl_output_all[['feature_q_value']] < 0.05 &
                                       qtl_output_all[['p_value']] < qtl_output_all[['pval_nominal_threshold_global']], ]

# location of the i-eqtl output
iqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/output/ut_and_24hca_significant/L1/'
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

# # location of the GTF
# gtf_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/refdata-cellranger-arc-GRCh38-2024-A/genes/genes.gtf.gz'
# gtf <- read.table(gtf_loc, header = F, sep = '\t')

# # and gencode
# gencode_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/ucsc_gencode/gencode43.tsv.gz'
# gencode <- read.table(gencode_loc, header = F, sep = '\t')
# colnames(gencode) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts')
# # add this info
# qtl_output_all_sig_wi[['strand']] <- gencode[match(
#   paste0('chr', qtl_output_all_sig_wi[['feature_chromosome']], '-', qtl_output_all_sig_wi[['feature_start']], '-', qtl_output_all_sig_wi[['feature_end']]),
#   paste0(gencode[['chrom']], '-', gencode[['chromStart']], '-', gencode[['chromEnd']])
# ), ][['strand']]

# # and the annotation of SCENIC+
# scenic_anno_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/genome_annotation.tsv'
# scenic_anno <- fread(scenic_anno_loc, header = T, sep = '\t')
# colnames(scenic_anno) <-c('chrom', 'start', 'end', 'strand', 'gs','Transcription_Start_Site','Transcript_type')
# # add to the QTL info the TSS
# qtl_output_all_sig_wi[['TSS']] <- scenic_anno[match(qtl_output_all_sig_wi[['feature_id']], scenic_anno[['gs']]), ][['Transcription_Start_Site']]
# # add to the QTL info the strand
# qtl_output_all_sig_wi[['strand']] <- scenic_anno[match(qtl_output_all_sig_wi[['feature_id']], scenic_anno[['gs']]), ][['strand']]

# location of the gene annotations with exons
exons_annotation_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/ncbi_gencode/GCF_000001405.40/genomic.gtf.gz'
# read the annotation file
exons_annotation <- read.table(exons_annotation_loc, header = F, sep = '\t')
# set column names
colnames(exons_annotation) <- c('seqname','source','feature','start','end','score','strand','frame','attributes')
# add the gene ID
exons_annotation[['gene_id']] <- str_extract(exons_annotation$attributes, "(?<=gene_id )[^;]+")
# add biotype
exons_annotation[['gene_biotype']] <- str_extract(exons_annotation$attributes, "(?<=gene_biotype )[^;]+")
# subset to the transcripts
exons_annotation_transcripts <- exons_annotation[exons_annotation[['feature']] == 'transcript', ]
# add a TSS
exons_annotation_transcripts[['TSS']] <- exons_annotation_transcripts[['start']]
# but make this the end if we are on the negative strand
exons_annotation_transcripts[exons_annotation_transcripts[['strand']] == '-', ][['TSS']] <- exons_annotation_transcripts[exons_annotation_transcripts[['strand']] == '-', ][['end']]
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

# get tss info
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/genome_annotation.tsv'
gene_anno <- fread(gene_anno_loc, header = T, sep = '\t')
# # rename columns to be the same as in limix
# colnames(gene_anno) <-c('chrom', 'start', 'end', 'strand', 'gs','Transcription_Start_Site','Transcript_type')
# # add TSS info to the QTL output as well
# qtl_output_all_sig_wi[['Transcription_Start_Site']] <- gene_anno[match(qtl_output_all_sig_wi[['feature_id']], gene_anno[['gs']]), ][['Transcription_Start_Site']]
# # get extra annotations for the pseudobulk output
# strand_information_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eQTA/LimixExpAnnotationFile.incStrand.txt'
# strand_information <- fread(strand_information_loc, header = T, sep = '\t')
# # add to the pseudobulk
# qtl_output_all_sig_wi[['strand']] <- strand_information[match(qtl_output_all_sig_wi[['feature_id']], strand_information[['feature_id']]), ][['strand']]
# # set the TSS
# qtl_output_all_sig_wi[['TSS']] <- qtl_output_all_sig_wi[['feature_start']]
# # to the end if the strand was negative
# qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['strand']]) & qtl_output_all_sig_wi[['strand']] == -1, ][['TSS']] <- qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['strand']]) & qtl_output_all_sig_wi[['strand']] == -1, ][['feature_end']]
# # set to NA if strand info was NA
# qtl_output_all_sig_wi[is.na(qtl_output_all_sig_wi[['strand']]), ][['TSS']] <- NA
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
# check if there is a difference
kruskal.test(tss_dist ~ interaction_direction, data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['tss_dist']]), ])
# for protein coding specifically
kruskal.test(tss_dist ~ interaction_direction, data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['tss_dist']]) & !is.na(qtl_output_all_sig_wi[['gene_biotype']]) & qtl_output_all_sig_wi[['gene_biotype']] == 'protein_coding', ])

# show all plots
plot_grid(
  p_snp_tss_interaction_distances, 
  p_snp_tss_interaction_distances_top, 
  p_snp_tss_interaction_distances_top_abs, 
  nrow = 2, 
  ncol = 2
)

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
kruskal.test(maf ~ gene_biotype, data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['distance_directional']]), ])
# p-value = 0.01058

plot_grid(
  # p_snp_gene_interaction_distances_directional_top, 
  p_snp_maf_directional_top, 
  p_snp_maf_biotype_top, 
  p_snp_maf_qltdirection_top
)

# get exon/intron info
exon_intron_annotation <- add_introns_to_exons_info(exons_annotation)
# add the intron/exon number with the actual annotation
exon_intron_annotation[['annotation']] <- paste(exon_intron_annotation[['feature']], exon_intron_annotation[['number']])

# now add the intron/exon info
qtl_output_all_sig_wi_gb <- add_gene_body_location(qtl_output_all_sig_wi, exon_intron_annotation)
# check how the variants are located
variant_gb_numbers <- data.frame(table(qtl_output_all_sig_wi_gb[, c('feature', 'number', 'interaction_direction')]))
# order
variant_gb_numbers <- variant_gb_numbers[order(variant_gb_numbers[['Freq']], decreasing = T), ]
# rename columns
colnames(variant_gb_numbers) <- c('intron_exon', 'intron_exon_nr', 'interaction', 'n_obs')
# remove zeroes
variant_gb_numbers <- variant_gb_numbers[variant_gb_numbers[['n_obs']] > 0, ]
