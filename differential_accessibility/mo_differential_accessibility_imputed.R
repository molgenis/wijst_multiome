#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_differential_accessibility_imputed.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

# read the object
library(Seurat)
library(Signac)
# plotting
library(ggplot2)
library(cowplot)
# limma DE dependencies
library(variancePartition)
library(edgeR)
library(BiocParallel)
library(doParallel)
# convert count matrices
library(Matrix)
#library(Matrix.utils) # NOT IN CONTAINER! als grr or remotes::install_github("cvarrichio/Matrix.utils")
library(optparse) # NOT IN CONTAINER


####################
# Functions        #
####################

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


get_peak_info <- function(peaks_object) {
  # calculate average expression
  avg_peaks <- AverageExpression(peaks_object)[['peaks']]
  # calculate average expression
  sum_peaks <- AggregateExpression(peaks_object)[['peaks']]
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
  return(peaks_bed)
}


plot_qc_metrics <- function(qc_data, plot_groups, plot_prepend, plot_values=c('nr', 'complexity')) {
  # make a plot list
  plot_list <- list()
  # do each group
  for (plot_group in plot_groups) {
    # check each plot value
    for (plot_value in plot_values) {
      # make plot
      p <- ggplot(data = NULL, mapping = aes(x = qc_data[[plot_value]]), fill = qc_data[[plot_group]]) +
        geom_density(alpha = 0.2) + 
        xlab(plot_value)
      # do plot
      p
      # save the plot
      ggsave(filename = paste(plot_prepend, '_', plot_group, '_', plot_value, '.pdf', sep = ''))
    }
  }
}


do_dream <- function(geneExpr, aggregate_metadata, cell_numbers, form, condition_combinations, verbose = T) {
  tryCatch(
    {
      # estimate weights using linear mixed model of dream
      vobjDream = voomWithDreamWeights( counts = geneExpr, formula = form, data = aggregate_metadata, weights = cell_numbers[['nr']] ) # the cell numbers are in the same order as the metadata, and as such can be passed like this
      
      # do each combination
      for(combination_name in names(condition_combinations)){
        # grab the combination
        combination <- condition_combinations[[combination_name]]
        
        if(verbose){
          print(paste('doing combination:', combination_name, sep = ' ', collapse = ' '))
        }
        tryCatch(
          {
            
            # define and then cbind contrasts
            L = getContrast( vobjDream, form, aggregate_metadata, paste(combination_name, combination, sep = ''))
            
            # fit contrast
            fit = dream( vobjDream, form, aggregate_metadata, L)
            
            # grab the exact fit
            limma_result <- topTable(fit, coef='L1', number=length(fit$F.p.value))
            
            
            # add bonferroni adjustment
            #limma_result[['p.bonferroni']] <- p.adjust(limma_result[['P.Value']])
            
            # add some statistics
            result_stats_list <- list()
            # check each condition
            result_stats_list[['combination']] <- data.frame(combination=rep(paste(combination_name, paste(combination, collapse='-'), sep = '.'), times = nrow(limma_result)))
            for (condition in combination) {
              # get the cell numbers for the condition
              cell_numbers_condition <- cell_numbers[aggregate_metadata[[combination_name]] == condition, ]
              result_stats_list[[paste('nsample', condition, sep = '_')]] <- data.frame(nsample=rep(nrow(cell_numbers_condition), times = nrow(limma_result)))
              result_stats_list[[paste('ncell', condition, sep = '_')]] <- data.frame(ncells=rep(
                paste(as.character(min(cell_numbers_condition[['nr']])),
                      as.character(quantile(cell_numbers_condition[['nr']])[['25%']]),
                      as.character(quantile(cell_numbers_condition[['nr']])[['50%']]),
                      as.character(quantile(cell_numbers_condition[['nr']])[['75%']]),
                      as.character(max(cell_numbers_condition[['nr']])),
                      sep = ';'
                ), times = nrow(limma_result)))
            }
            # merge all
            result_stats <- do.call('cbind', result_stats_list)
            colnames(result_stats) <- names(result_stats_list)
            
            # add combination as first column
            limma_result <- cbind(result_stats, limma_result)
            
            # finally also add the feature as an explicit column
            limma_result <- cbind(data.frame(feature = rownames(limma_result)), limma_result)
            
            # return the result
            return(limma_result)
          }, error=function(cond) {
            print(paste('analysis failed in', combination))
            message(cond)
          }
        )
      }
    }, error=function(cond) {
      print(paste('model build failed'))
      message(cond)
    }
  )
}


#' perform pseudobulk limma
#' 
#' @param seurat_object The Seurat object to add the sample assignment to
#' @param output_loc the location to write the output tables to
#' @param condition_combinations named list of vectors which to pairwise compare e.g. list('condition_final' = c('24hCA', 'UT'))
#' @param aggregates the columns to aggregate on, so pseudobulk per donor and inflammation is c('donor', 'inflammation')
#' @param fixed_effects vector of fixed effects to include in the model
#' @param random_effects vector or random effects to include in the model
#' @param minimal_cells the minimal number of cells that needs to be present for a sample to not be excluded
#' @param min_peaks the minimal number of UMIs a cell must have to be used for the pseudobulk
#' @param minimal_complexity the minimal number of cells transcripts a pseudobulk needs to be based on for a sample to not be excluded
#' @param plot_metrics plot the QC metrics that were used for filtering
#' @param verbose whether to print progress messages
#' @param permute do a permutation round instead
#' @param seed use supplied seed if permuting. If not supplied a random seed is selected (recommended)
#' @returns 0 if successful
#' dream_pairwise(pbmc, './bulk_test/', list('inflammation'=c('AI','NI')))
dream_pairwise_mt <- function(seurat_object, output_loc, condition_combinations, aggregates=c('sample_final', 'condition_final'), fixed_effects=c('condition_final'), random_effects=c('sample_final'), minimal_cells=0, min_peaks=200, minimal_complexity=5000, min_pct=0.01, verbose=T, plot_metrics=T, nthreads=5, permute=F, seed=NULL){
  # set assay
  DefaultAssay(seurat_object) <- 'peaks'
  # grab the countmatrix
  countMatrix <- NULL
  # subset the object based on the minimal number of umis if requested
  if (min_peaks > 0) {
    seurat_object <- seurat_object[, seurat_object@meta.data[['nCount_peaks']] >= min_peaks]
  }
  # grab the peak info
  peak_info <- get_peak_info(seurat_object)
  # filter on the percentage of cells having the peak
  if (min_pct > 0) {
    # get the names that survive the filter
    peak_names_filtered <- peak_info[peak_info[['pct_exp']] > min_pct, 'name']
    # filter by rownames
    seurat_object <- seurat_object[peak_names_filtered, ]
  }
  # depending on the version
  if (grepl('^3|4.0', seurat_object@version)) {
    countMatrix <- seurat_object@assays$peaks@counts
  }
  else if (grepl('^4.9|5', seurat_object@version)) {
    countMatrix <- seurat_object@assays$peaks@counts
  }
  # get the unique combinations of aggregates and effects
  unique_aggragation_names <- unique(c(aggregates, fixed_effects, random_effects))
  # add this as a new column to the metadata
  seurat_object@meta.data[['aggregate']] <- seurat_object@meta.data[[unique_aggragation_names[1]]]
  if (length(unique_aggragation_names) > 1) {
    for (i in 2:length(unique_aggragation_names)) {
      seurat_object@meta.data[['aggregate']] <- paste(seurat_object@meta.data[['aggregate']], seurat_object@meta.data[[unique_aggragation_names[i]]], sep = '-')
    }
  }
  # get the metadata
  metadata <- seurat_object@meta.data
  # set that as the ident
  Idents(seurat_object) <- 'aggregate'
  # fetch unique metadata
  aggregate_metadata <- unique(metadata[, c(unique_aggragation_names, 'aggregate')])
  # get the aggregated count matrix
  aggregate_countMatrix <- AggregateExpression(seurat_object)[['peaks']]
  # order the metadata by the count matrix order
  rownames(aggregate_metadata) <- aggregate_metadata[['aggregate']]
  aggregate_metadata <- aggregate_metadata[colnames(aggregate_countMatrix), ]
  
  # next get the cell numbers for each observation, these will follow the order of the original aggregated metadata
  cell_numbers <- data.frame(table(metadata[['aggregate']]))
  # set colnames
  colnames(cell_numbers) <- c('aggregate', 'nr')
  # order the same as the metadata
  cell_numbers <- cell_numbers[match(aggregate_metadata[['aggregate']], cell_numbers[['aggregate']]), ]
  
  # filter by the number of cells if requested
  if (minimal_cells > 0) {
    # get the indices of where the cell numbers are above this
    indices_above_threshold <- which(cell_numbers[['nr']] >= minimal_cells)
    # report how many
    if (verbose) {
      message(paste('of', as.character(nrow(cell_numbers)), 'entries, ', length(indices_above_threshold), 'contained more cells than the', as.character(minimal_cells), 'threshold'))
    }
    # do the actual filtering
    cell_numbers <- cell_numbers[indices_above_threshold, ]
    aggregate_countMatrix <- aggregate_countMatrix[, indices_above_threshold]
    aggregate_metadata <- aggregate_metadata[indices_above_threshold, ]
  }
  # calculate the complexity first
  complexity <- data.frame((colSums(aggregate_countMatrix)))
  # filter by the complexity if requested
  if (minimal_complexity > 0) {
    # set better column names
    colnames(complexity) <- c('complexity')
    # get the indices of where the complexity is above the threshold
    indices_above_complexity <- which(complexity[['complexity']] >= minimal_complexity)
    # report how many
    if (verbose) {
      message(paste('of', as.character(nrow(complexity)), 'entries, ', length(indices_above_complexity), 'contained more fragments than the', as.character(minimal_complexity), 'threshold'))
    }
    # and filter
    cell_numbers <- cell_numbers[indices_above_complexity, ]
    aggregate_countMatrix <- aggregate_countMatrix[, indices_above_complexity]
    aggregate_metadata <- aggregate_metadata[indices_above_complexity, ]
    complexity <- complexity[indices_above_complexity, ]
  }
  
  # filter genes by number of counts
  isexpr = rowSums(cpm(aggregate_countMatrix)>0.1) >= 5
  
  # Standard usage of limma/voom
  geneExpr = DGEList( aggregate_countMatrix[isexpr,] )
  geneExpr = calcNormFactors( geneExpr )
  
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(4, "SOCK", progressbar=TRUE)
  register(param)
  
  # show head of the tables if we are being verbose
  if(verbose){
    message('aggregated metadata head:')
    print(head(aggregate_metadata))
    message('aggregated counts head:')
    print(head(geneExpr[['counts']]))
    message('cell numbers head')
    #print(head(cell_numbers[order(cell_numbers[['nr']]), ]))
    print(head(cell_numbers))
    message('complexity head:')
    print(head(data.frame('aggregate' = cell_numbers[['aggregate']], 'complexity' = complexity)))
  }
  # plot the metrics if they were requested
  if (plot_metrics) {
    # add the cell counts just to the aggregated metadata
    qc_metrics_table <- cbind(aggregate_metadata, cell_numbers[, c('nr'), drop = F])
    # add the complexity as well
    qc_metrics_table <- cbind(qc_metrics_table, data.frame('complexity' = complexity))
    # make the plots
    plot_qc_metrics(qc_data = qc_metrics_table, plot_groups = fixed_effects, plot_prepend = paste(output_loc, names(condition_combinations)[[1]], sep = ''), plot_values=c('nr', 'complexity'))
  }
  
  # paste together the model
  model_formula <- '~ 0 '
  for(fixed_effect in fixed_effects){
    model_formula <- paste(model_formula, fixed_effect, sep = ' + ')
  }
  for(random_effect in random_effects){
    model_formula <- paste(model_formula, ' + (1|', random_effect, ')', sep = '')
  }
  
  if(verbose){
    message(paste('formula:', model_formula))
  }
  
  # and turn into a formula
  form <- as.formula(model_formula)
  
  # if we are doing a permutation run, we'll have to change some metadata. We'll use a different variable, so we keep the original as a backup
  aggregate_metadata_to_use <- aggregate_metadata
  # the permutation seed we'll save, though it will only be used if there is an actual permutation
  permutation_seed <- NA
  # let's set up some things for the permutation
  if (permute) {
    # if we have a seed, we'll use that, otherwise we need to get one
    if (!(is.null(seed))) {
      permutation_seed <- seed
      if(verbose){
        message(paste('permuting with set seed of', as.character(permutation_seed), sep = ' ', collapse = ' '))
      }
    }
    else {
      # get a seed
      set.seed(NULL)
      permutation_seed <- round(.Machine$integer.max * runif(1), digits = 0)
      if(verbose){
        message(paste('permuting with random seed of', as.character(permutation_seed), sep = ' ', collapse = ' '))
      }
    }
    # set this seed
    set.seed(permutation_seed)
    if (length(fixed_effects > 1)) {
      warning(paste('when permuting with multithreading enabled, only the first fixed effect can be permuted, in this case using', fixed_effects[1]))
    }
    # then permute the combination we are looking at
    aggregate_metadata_to_use[[fixed_effects[1]]] <- sample(aggregate_metadata_to_use[[fixed_effects[1]]], size = nrow(aggregate_metadata_to_use), replace = F)
  }
  
  # get the number of rows
  nrow_matrix <- nrow(geneExpr)
  # get the number of chunks
  n_chunks <- nthreads
  # round up for the number of threads
  n_chunks <- ceiling(n_chunks)
  # get the number of rows per chunk
  nrow_chunk <- nrow(geneExpr) / n_chunks
  # now do parallel processing of chunks
  res_per_chunk <- foreach(i = 1:n_chunks) %dopar% {
    # calculate the chunk start and stop
    chunk_start <- (i - 1) * nrow_chunk + 1
    chunk_stop <- i * nrow_chunk
    # check if we are exceeding the number of rows
    if (chunk_stop > nrow_matrix) {
      # because then we will stop there
      chunk_stop <- nrow_matrix
    }
    # extract these rows
    geneExprChunk <- geneExpr[chunk_start : chunk_stop, ]
    # do the limma run for this chunk
    limma_result_chunk <- do_dream(geneExprChunk, aggregate_metadata_to_use, cell_numbers, form, condition_combinations, verbose)
    return(limma_result_chunk)
  }
  # merge chuncks
  limma_result <- do.call('rbind', res_per_chunk)
  
  # redo the B&H, as it was done per chunk, while it needs to be done on all data
  limma_result[['adj.P.Val']] <- p.adjust(limma_result[['P.Value']], method = 'BH')
  
  # add bonferroni correction
  limma_result[['p.bonferroni']] <- p.adjust(limma_result[['P.Value']], method = 'bonferroni')
  
  # now add the permutation status and seed
  limma_result[['permuted']] <- permute
  limma_result[['seed']] <- permutation_seed
  
  # set an output location
  limma_output_loc <- (paste(output_loc, names(condition_combinations)[[1]], '.tsv.gz', sep = ''))
  
  # which will be a bit different when we are permuting
  if (permute) {
    limma_output_loc <- (paste(output_loc, names(condition_combinations)[[1]], '.', permutation_seed, '.tsv.gz', sep = ''))
    # we'll also order by the name of the feature
    limma_result <- limma_result[order(limma_result[['feature']]), ]
  }
  else {
    # otherwise we order by the significance
    limma_result <- limma_result[order(limma_result[['P.Value']]), ]
  }
  
  # also write the model we used
  limma_formula_loc <- paste(output_loc, names(condition_combinations)[[1]], '.formula', sep = '')
  
  if(verbose){
    message(paste('writing result', limma_output_loc, sep = ' '))
  }
  
  # write the result
  write.table(limma_result, gzfile(limma_output_loc), sep = '\t', row.names = F)
  # and the formula
  write.table(model_formula, limma_formula_loc, row.names = F, col.names = F)
  # write md5
  mdfiver::create_md5_for_file(limma_output_loc)
  
  # now write just the nominally significant values as well, but only if we don't permute
  if (!permute) {
    limma_result <- limma_result[limma_result[['P.Value']] < 0.05, ]
    limma_output_nomsig_loc <- gzfile(paste(output_loc, names(condition_combinations)[[1]], '.nominal_significant.tsv.gz', sep = ''))
    write.table(limma_result, limma_output_nomsig_loc, sep = '\t', row.names = F)
    mdfiver::create_md5_for_file(paste(output_loc, names(condition_combinations)[[1]], '.nominal_significant.tsv.gz', sep = ''))
  }
  return(0)
}


#' perform pseudobulk limma
#' 
#' @param seurat_object The Seurat object to add the sample assignment to
#' @param output_loc the location to write the output tables to
#' @param condition_combinations named list of vectors which to pairwise compare e.g. list('condition_final' = c('24hCA', 'UT'))
#' @param aggregates the columns to aggregate on, so pseudobulk per donor and inflammation is c('donor', 'inflammation')
#' @param fixed_effects vector of fixed effects to include in the model
#' @param random_effects vector or random effects to include in the model
#' @param minimal_cells the minimal number of cells that needs to be present for a sample to not be excluded
#' @param min_peaks the minimal number of UMIs a cell must have to be used for the pseudobulk
#' @param minimal_complexity the minimal number of cells transcripts a pseudobulk needs to be based on for a sample to not be excluded
#' @param verbose whether to print progress messages
#' @param plot_metrics plot the QC metrics that were used for filtering
#' @param permute do a permutation round instead
#' @param seed use supplied seed if permuting. If not supplied a random seed is selected (recommended)
#' @returns 0 if successful
#' dream_pairwise(pbmc, './bulk_test/', list('inflammation'=c('AI','NI')))
dream_pairwise <- function(seurat_object, output_loc, condition_combinations, aggregates=c('sample_final', 'condition_final'), fixed_effects=c('condition_final'), random_effects=c('sample_final'), minimal_cells=0, min_peaks=200, minimal_complexity=5000, min_pct=0.01, verbose=T, plot_metrics=T, permute=F, seed=NULL, nthreads=20){
  # set assay
  DefaultAssay(seurat_object) <- 'peaks'
  # grab the countmatrix
  countMatrix <- NULL
  # subset the object based on the minimal number of umis if requested
  if (min_peaks > 0) {
    seurat_object <- seurat_object[, seurat_object@meta.data[['nCount_peaks']] >= min_peaks]
  }
  # grab the peak info
  peak_info <- get_peak_info(seurat_object)
  # filter on the percentage of cells having the peak
  if (min_pct > 0) {
    # get the names that survive the filter
    peak_names_filtered <- peak_info[peak_info[['pct_exp']] > min_pct, 'name']
    # filter by rownames
    seurat_object <- seurat_object[peak_names_filtered, ]
  }
  # depending on the version
  if (grepl('^3|4.0', seurat_object@version)) {
    countMatrix <- seurat_object@assays$peaks@counts
  }
  else if (grepl('^4.9|5', seurat_object@version)) {
    countMatrix <- seurat_object@assays$peaks@counts
  }
  # get the unique combinations of aggregates and effects
  unique_aggragation_names <- unique(c(aggregates, fixed_effects, random_effects))
  # add this as a new column to the metadata
  seurat_object@meta.data[['aggregate']] <- seurat_object@meta.data[[unique_aggragation_names[1]]]
  if (length(unique_aggragation_names) > 1) {
    for (i in 2:length(unique_aggragation_names)) {
      seurat_object@meta.data[['aggregate']] <- paste(seurat_object@meta.data[['aggregate']], seurat_object@meta.data[[unique_aggragation_names[i]]], sep = '-')
    }
  }
  # get the metadata
  metadata <- seurat_object@meta.data
  # set that as the ident
  Idents(seurat_object) <- 'aggregate'
  # fetch unique metadata
  aggregate_metadata <- unique(metadata[, c(unique_aggragation_names, 'aggregate')])
  # get the aggregated count matrix
  aggregate_countMatrix <- AggregateExpression(seurat_object)[['peaks']]
  # order the metadata by the count matrix order
  rownames(aggregate_metadata) <- aggregate_metadata[['aggregate']]
  aggregate_metadata <- aggregate_metadata[colnames(aggregate_countMatrix), ]
  
  # create the groups to aggregate on, here it's on sample and inflammation status usually
  #groups <- metadata[, unique(c(aggregates, fixed_effects, random_effects))]
  # create an aggregated counts matrix
  #aggregate_countMatrix <- t(aggregate.Matrix(t(countMatrix), groupings = groups, fun = 'sum'))
  # create aggregated metadata
  #aggregate_metadata <- unique(metadata[, unique(c(aggregates, fixed_effects, random_effects))])
  # set the rownames of the aggregate metadata
  #rownames_to_set_agg_metadata <- aggregate_metadata[[unique(c(aggregates, fixed_effects, random_effects))[1]]]
  #for(i in 2:length(unique(c(aggregates, fixed_effects, random_effects)))){
  #  rownames_to_set_agg_metadata <- paste(rownames_to_set_agg_metadata, aggregate_metadata[[unique(c(aggregates, fixed_effects, random_effects))[i]]], sep='_')
  #}
  #rownames(aggregate_metadata) <- rownames_to_set_agg_metadata
  # set in the same order as the count matrix
  #aggregate_metadata <- aggregate_metadata[colnames(aggregate_countMatrix), ]
  
  # next get the cell numbers for each observation, these will follow the order of the original aggregated metadata
  #cell_numbers <- get_nr_cells_aggregate_combination(aggregate_metadata[, unique(c(aggregates, fixed_effects, random_effects))], metadata)
  cell_numbers <- data.frame(table(metadata[['aggregate']]))
  # set colnames
  colnames(cell_numbers) <- c('aggregate', 'nr')
  # order the same as the metadata
  cell_numbers <- cell_numbers[match(aggregate_metadata[['aggregate']], cell_numbers[['aggregate']]), ]
  
  # filter by the number of cells if requested
  if (minimal_cells > 0) {
    # get the indices of where the cell numbers are above this
    indices_above_threshold <- which(cell_numbers[['nr']] >= minimal_cells)
    # report how many
    if (verbose) {
      message(paste('of', as.character(nrow(cell_numbers)), 'entries, ', length(indices_above_threshold), 'contained more cells than the', as.character(minimal_cells), 'threshold'))
    }
    # do the actual filtering
    cell_numbers <- cell_numbers[indices_above_threshold, ]
    aggregate_countMatrix <- aggregate_countMatrix[, indices_above_threshold]
    aggregate_metadata <- aggregate_metadata[indices_above_threshold, ]
  }
  # calculate the complexity first
  complexity <- data.frame((colSums(aggregate_countMatrix)))
  # filter by the complexity if requested
  if (minimal_complexity > 0) {
    # print complexity if we are using it as a filter
    print(head(complexity))
    # set better column names
    colnames(complexity) <- c('complexity')
    # get the indices of where the complexity is above the threshold
    indices_above_complexity <- which(complexity[['complexity']] >= minimal_complexity)
    # report how many
    if (verbose) {
      message(paste('of', as.character(nrow(complexity)), 'entries, ', length(indices_above_complexity), 'contained more fragments than the', as.character(minimal_complexity), 'threshold'))
    }
    # and filter
    cell_numbers <- cell_numbers[indices_above_complexity, ]
    aggregate_countMatrix <- aggregate_countMatrix[, indices_above_complexity]
    aggregate_metadata <- aggregate_metadata[indices_above_complexity, ]
    complexity <- complexity[indices_above_complexity, ]
  }
  
  # filter genes by number of counts
  isexpr = rowSums(cpm(aggregate_countMatrix)>0.1) >= 5
  
  # Standard usage of limma/voom
  geneExpr = DGEList( aggregate_countMatrix[isexpr,] )
  geneExpr = calcNormFactors( geneExpr )
  
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(nthreads, "SOCK", progressbar=TRUE)
  register(param)
  
  # show head of the tables if we are being verbose
  if(verbose){
    message('aggregated metadata head:')
    print(head(aggregate_metadata))
    message('aggregated counts head:')
    print(head(geneExpr[['counts']]))
    message('cell numbers head')
    #print(head(cell_numbers[order(cell_numbers[['nr']]), ]))
    print(head(cell_numbers))
  }
  # plot the metrics if they were requested
  if (plot_metrics) {
    # add the cell counts just to the aggregated metadata
    qc_metrics_table <- cbind(aggregate_metadata, cell_numbers[, c('nr'), drop = F])
    # add the complexity as well
    qc_metrics_table <- cbind(qc_metrics_table, data.frame('complexity' = complexity[[1]]))
    # make the plots
    plot_qc_metrics(qc_data = qc_metrics_table, plot_groups = fixed_effects, plot_prepend = paste(output_loc, names(condition_combinations)[[1]], sep = ''), plot_values=c('nr', 'complexity'))
  }
  
  # paste together the model
  model_formula <- '~ 0 '
  for(fixed_effect in fixed_effects){
    model_formula <- paste(model_formula, fixed_effect, sep = ' + ')
  }
  for(random_effect in random_effects){
    model_formula <- paste(model_formula, ' + (1|', random_effect, ')', sep = '')
  }
  
  if(verbose){
    message(paste('formula:', model_formula))
  }
  
  # and turn into a formula
  form <- as.formula(model_formula)
  # tryCatch(
  #   {
  # estimate weights using linear mixed model of dream
  vobjDream = voomWithDreamWeights( counts = geneExpr, formula = form, data = aggregate_metadata, weights = cell_numbers[['nr']] ) # the cell numbers are in the same order as the metadata, and as such can be passed like this
  
  # do each combination
  for(combination_name in names(condition_combinations)){
    # grab the combination
    combination <- condition_combinations[[combination_name]]
    
    if(verbose){
      message(paste('doing combination:', combination_name, 'with coefficients', paste(combination, collapse = ','), sep = ' ', collapse = ' '))
    }
    
    # if we are doing a permutation run, we'll have to change some metadata. We'll use a different variable, so we keep the original as a backup
    aggregate_metadata_to_use <- aggregate_metadata
    # the permutation seed we'll save, though it will only be used if there is an actual permutation
    permutation_seed <- NA
    # let's set up some things for the permutation
    if (permute) {
      # if we have a seed, we'll use that, otherwise we need to get one
      if (!is.null(seed)) {
        permutation_seed <- seed
        if(verbose){
          message(paste('permuting with set seed of', as.character(permutation_seed), sep = ' ', collapse = ' '))
        }
      }
      else {
        # get a seed
        set.seed(NULL)
        permutation_seed <- round(.Machine$integer.max * runif(1), digits = 0)
        if(verbose){
          message(paste('permuting with random seed of', as.character(permutation_seed), sep = ' ', collapse = ' '))
        }
      }
      # set this seed
      set.seed(permutation_seed)
      # then permute the combination we are looking at
      aggregate_metadata_to_use[[combination_name]] <- sample(aggregate_metadata_to_use[[combination_name]], size = nrow(aggregate_metadata_to_use), replace = F)
    }
    
    # tryCatch(
    #   {
    # define and then cbind contrasts
    L = getContrast( vobjDream, form, aggregate_metadata_to_use, paste(combination_name, combination, sep=''))
    
    # fit contrast
    fit = dream( vobjDream, form, aggregate_metadata_to_use, L)
    
    # grab the exact fit
    limma_result <- topTable(fit, coef='L1', number=length(fit$F.p.value))
    
    # add bonferroni adjustment
    limma_result[['p.bonferroni']] <- p.adjust(limma_result[['P.Value']], method = 'bonferroni')
    
    # add some statistics
    result_stats_list <- list()
    # check each condition
    result_stats_list[['combination']] <- data.frame(combination=rep(paste(combination_name, paste(combination, collapse='-'), sep = '.'), times = nrow(limma_result)))
    for (condition in combination) {
      # get the cell numbers for the condition
      cell_numbers_condition <- cell_numbers[cell_numbers[[combination_name]] == condition, ]
      result_stats_list[[paste('nsample', condition, sep = '_')]] <- data.frame(nsample=rep(nrow(cell_numbers_condition), times = nrow(limma_result)))
      result_stats_list[[paste('ncell', condition, sep = '_')]] <- data.frame(ncells=rep(
        paste(as.character(min(cell_numbers_condition[['nr']])),
              as.character(quantile(cell_numbers_condition[['nr']])[['25%']]),
              as.character(quantile(cell_numbers_condition[['nr']])[['50%']]),
              as.character(quantile(cell_numbers_condition[['nr']])[['75%']]),
              as.character(max(cell_numbers_condition[['nr']])),
              sep = ';'
        ), times = nrow(limma_result)))
    }
    # merge all
    result_stats <- do.call('cbind', result_stats_list)
    colnames(result_stats) <- names(result_stats_list)
    
    # now add the permutation status and seed
    result_stats[['permuted']] <- permute
    result_stats[['seed']] <- permutation_seed
    
    # add combination as first column
    limma_result <- cbind(result_stats, limma_result)
    
    # finally also add the feature as an explicit column
    limma_result <- cbind(data.frame(feature = rownames(limma_result)), limma_result)
    
    # set an output location
    limma_output_loc <- (paste(output_loc, combination_name, '.tsv.gz', sep = ''))
    
    # which will be a bit different when we are permuting
    if (permute) {
      limma_output_loc <- (paste(output_loc, combination_name, '.', permutation_seed, '.tsv.gz', sep = ''))
      # we'll also order by the name of the feature
      limma_result <- limma_result[order(limma_result[['feature']]), ]
    }
    else {
      # otherwise we order by the significance
      limma_result <- limma_result[order(limma_result[['P.Value']]), ]
    }
    
    # also write the model we used
    limma_formula_loc <- paste(output_loc, combination_name, '.formula', sep = '')
    
    if(verbose){
      message(paste('writing result', limma_output_loc))
    }
    
    # write the result
    write.table(limma_result, gzfile(limma_output_loc), sep = '\t', row.names = F)
    # and the formula
    write.table(model_formula, limma_formula_loc, row.names = F, col.names = F)
    # write md5
    mdfiver::create_md5_for_file(limma_output_loc)
    
    # now write just the nominally significant values as well, if we don't permute
    if (!permute) {
      limma_result <- limma_result[limma_result[['P.Value']] < 0.05, ]
      limma_output_nomsig_loc <- (paste(output_loc, names(condition_combinations)[[1]], '.nominal_significant.tsv.gz', sep = ''))
      write.table(limma_result, gzfile(limma_output_nomsig_loc), sep = '\t', row.names = F)
      mdfiver::create_md5_for_file(paste(output_loc, names(condition_combinations)[[1]], '.nominal_significant.tsv.gz', sep = ''))
    }
    #   }, error=function(cond) {
    #     print(paste('analysis failed in', combination))
    #     message(cond)
    #   }
    # )
  }
  #   }, error=function(cond) {
  #     print(paste('model build failed'))
  #     message(cond)
  #   }
  # )
  return(0)
}


#' perform pseudobulk limma per celltype
#' 
#' @param seurat_object The Seurat object to add the sample assignment to
#' @param output_loc the location to write the output tables to
#' @param condition_combinations named list of vectors which to pairwise compare e.g. list('condition_final' = c('24hCA', 'UT'))
#' @param celltype_column the column in the metadata that contains the cell type to perform limma on
#' @param cell_type_to_use the cell types to perform limma on. If NULL, all cell types will be used
#' @param aggregates the columns to aggregate on, so pseudobulk per donor and inflammation is c('donor', 'inflammation')
#' @param fixed_effects vector of fixed effects to include in the model
#' @param random_effects vector or random effects to include in the model
#' @param minimal_cells the minimal number of cells that needs to be present for a sample to not be excluded
#' @param min_peaks the minimal number of UMIs a cell must have to be used for the pseudobulk
#' @param minimal_complexity the minimal number of cells transcripts a pseudobulk needs to be based on for a sample to not be excluded
#' @param verbose whether to print progress messages
#' @param nthreads number of threads
#' @param permute whether to do a permutation-round instead
#' @param seed use supplied seed if permuting. If not supplied a random seed is selected (recommended)
#' @returns 0 if successful
#' do_limma_dream_pairwise_per_celltype(pbmc, './ct_test/')
do_limma_dream_pairwise_per_celltype <- function(seurat_object, output_loc, condition_combinations=list('condition_final' =  c('24hCA', 'UT')), celltype_column='cell_type_final', cell_types_to_use=NULL, aggregates=c('sample_final', 'condition_final'), fixed_effects=c('condition_final'), random_effects=c('sample_final'), minimal_cells=0, min_peaks=200, minimal_complexity=5000, verbose=T, nthreads=5, permute=F, seed=NULL){
  # use the cell types supplied, or all if none are supplied
  cell_types <- cell_types_to_use
  if(is.null(cell_types_to_use)){
    cell_types <- unique(seurat_object@meta.data[[celltype_column]])
  }
  # remove NA ones
  cell_types <- cell_types[!is.na(cell_types)]
  # check each cell type
  for(cell_type in cell_types){
    if (verbose) {
      print(paste('doing', cell_type))
    }
    # get a POSIX-safe celltype name
    cell_type_safe <- gsub(' |/', '_', cell_type)
    cell_type_safe <- gsub('-', '_negative', cell_type_safe)
    cell_type_safe <- gsub('\\+', '_positive', cell_type_safe)
    cell_type_safe <- gsub('\\)', '', cell_type_safe)
    cell_type_safe <- gsub('\\(', '', cell_type_safe)
    # make that into the output prepend
    output_loc_celltype <- paste(output_loc, cell_type_safe, '_', sep = '')
    # subset to the cell type
    seurat_object_celltype <- seurat_object[, seurat_object@meta.data[[celltype_column]] == cell_type]
    # perform the analysis
    # if (nthreads == 1) {
    #  dream_pairwise(seurat_object = seurat_object_celltype, output_loc = output_loc_celltype, condition_combinations = condition_combinations, aggregates = aggregates, fixed_effects = fixed_effects, random_effects = random_effects, minimal_cells = minimal_cells, min_peaks = min_peaks, minimal_complexity = minimal_complexity, verbose = verbose, permute = permute, seed = seed, nthreads = nthreads)
    # }
    # else if (nthreads > 1) {
    dream_pairwise_mt(seurat_object = seurat_object_celltype, output_loc = output_loc_celltype, condition_combinations = condition_combinations, aggregates = aggregates, fixed_effects = fixed_effects, random_effects = random_effects, minimal_cells = minimal_cells, min_peaks = min_peaks, minimal_complexity = minimal_complexity, verbose = verbose, nthreads = nthreads, permute = permute, seed = seed)
    # }
    # else {
    #   stop(paste('nthreads should be a positive number, now is', as.character(nthreads)))
    # }
  }
  return(0)
}

#' remove illegal characters from the cell type annotation
#' 
#' @param cell_types vector of cell types
#' @returns vector of celltypes with safe name
#' make_celltypes_safe(c('CD4+T', 'NK-dim'))
make_celltypes_safe <- function(cell_types){
  # get a safe file name
  cell_type_safes <- gsub(' |/', '_', cell_types)
  cell_type_safes <- gsub('-', '_negative', cell_type_safes)
  cell_type_safes <- gsub('\\+', '_positive', cell_type_safes)
  cell_type_safes <- gsub('\\)', '', cell_type_safes)
  cell_type_safes <- gsub('\\(', '', cell_type_safes)
  return(cell_type_safes)
}

#' remove illegal characters from the cell type annotation
#' 
#' @param seurat_object the seurat object to add a lower classification of cell type to
#' @param reclassification_mapping a table that has the original cell type, and what it should be changed to
#' @param mapping_original_column the column in the mapping table with the original cell types
#' @param mapping_reclass_column the column in the mapping table with what to change the cell type to
#' @param metadata_original_column the column in the metadata to get the original cell types from
#' @param metadata_reclassification_column the metadata column to add with the new mapping
#' @returns vector of celltypes with safe name
#' pmbc <- add_lower_classification(pmbc, mapping_table, 'high_res_ct', 'low_res_ct', 'cell.type', 'cell.type.low')
add_lower_classification <- function(seurat_object, reclassification_mapping, mapping_original_column, mapping_reclass_column, metadata_original_column, metadata_reclassification_column){
  # add the new column
  seurat_object@meta.data[[metadata_reclassification_column]] <- NA
  # get each cell type in the data
  metadata_original_cts <- unique(seurat_object@meta.data[[metadata_original_column]])
  # and the originals in the mapping
  reclassification_original_cts <- unique(reclassification_mapping[[mapping_original_column]])
  # we can only map what is present in both
  originals_both <- intersect(metadata_original_cts, reclassification_original_cts)
  # check what is missing
  only_metadata <- setdiff(metadata_original_cts, reclassification_original_cts)
  only_mapping <- setdiff(reclassification_original_cts, metadata_original_cts)
  # warn what is missing
  if(length(only_metadata) > 0){
    print('some celltypes only in metadata')
    print(only_metadata)
  }
  if(length(only_mapping) > 0){
    print('some celltypes only in remapping ')
    print(only_mapping)
  }
  # check each cell type
  for(celltype_original in originals_both){
    # get the appropriate remapping
    celltype_remapped <- reclassification_mapping[reclassification_mapping[[mapping_original_column]] == celltype_original, mapping_reclass_column]
    # now remap in the metadata
    seurat_object@meta.data[seurat_object@meta.data[[metadata_original_column]] == celltype_original, metadata_reclassification_column] <- celltype_remapped
  }
  return(seurat_object)
}

#' create a column which combines the values of multiple other columns
#' 
#' @param dataframe dataframe that has the original column, to which to add the aggregate column
#' @param columns the columns to combine into a new column
#' @returns the original dataframe, with a new column 'all_aggregates', which is the combination of the supplied columns
#' new_df <- combine_columns(df, c('V1', 'V3'))
combine_columns <- function(dataframe, columns) {
  # check each column
  for (column in columns) {
    # check if we were already aggregating
    if ('all_aggregates' %in% colnames(dataframe)) {
      # add the new column
      dataframe[['all_aggregates']] <- paste(dataframe[['all_aggregates']], dataframe[[column]], sep = '-')
    }
    # otherwise we need to start our aggregation
    else{
      dataframe[['all_aggregates']] <- dataframe[[column]]
    }
  }
  return(dataframe)
}

#' get the number of cells describing a combination of values from the seurat metadata
#' 
#' @param aggregate_df the dataframe with the combinations that were aggregated over
#' @param seurat_metadata the metadata of a seurat object from which to get cell numbers for the combinations
#' @returns the original dataframe, with a 'nr' column, denoting the number of cells for that combination
#' cell_numbers <- get_nr_cells_aggregate_combination(aggretates, pbmc_meta.data)
get_nr_cells_aggregate_combination <- function(aggregate_df, seurat_metadata) {
  # grab the columns we aggregated on
  aggregate_columns <- colnames(aggregate_df)
  # now use table to get the number of entries of these aggregates in the metadata
  aggregate_numbers <- data.frame(table(seurat_metadata[, aggregate_columns, drop = F]))
  # add a new column for the combination of the aggregates
  aggregate_df <- combine_columns(aggregate_df, aggregate_columns)
  # in our numbers as well
  aggregate_numbers <- combine_columns(aggregate_numbers, aggregate_columns)
  # now join the frequencies from the aggragate numbers onto the aggregate columns
  aggregate_df[['nr']] <- aggregate_numbers[match(aggregate_df[['all_aggregates']], aggregate_numbers[['all_aggregates']]), 'Freq']
  # now we are sure that the cell numbers are in the same order as the aggregated metadata
  # remove the column we created
  aggregate_df[['all_aggregates']] <- NULL
  # and return the result
  return(aggregate_df)
}


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 900 * 1000 * 1024^2)

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# the cell type
cell_type <- 'monocyte'

# location of the object
signac_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_monocyte_wstatus_1_80_20240709.rds'
signac_celltype <- readRDS(signac_object_loc)

# and parameters
nthreads <- 4
topic_ann_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/topic_annotations/mo_topic_20_otsu.tsv'
limma_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/output/topics20_otsu_imputed/'
celltype_column <- 'cell_type'
min_cells <- 10
min_cell_umis <- 200
min_pseudo_umis <- 10000
permute <- F

# set number of parallel threads
parallel::mcaffinity(1:nthreads)
registerDoParallel(cores=nthreads)

# replace underscore with dash for the lane
signac_celltype@meta.data[['lane']] <- gsub('_', '-', signac_celltype@meta.data[['lane']])
# add the day to the metadata
signac_celltype@meta.data[['day']] <- gsub('_lane\\d+', '', signac_celltype@meta.data[['lane']])

# filter where we don't have the sex
signac_celltype <- signac_celltype[, !is.na(signac_celltype@meta.data[['sex']])]
# or the age
signac_celltype <- signac_celltype[, !is.na(signac_celltype@meta.data[['age']])]
# or the inflammation status
signac_celltype <- signac_celltype[, !is.na(signac_celltype@meta.data[['condition_final']]) & signac_celltype@meta.data[['condition_final']] != 'unknown']
# add c to conditions so we won't have issues with the aggregation
signac_celltype@meta.data[['condition_final']] <- paste('c', signac_celltype@meta.data[['condition_final']], sep = '')

# read the topic annotation file
topic_ann <- read.table(topic_ann_loc, header = T, sep = '\t', row.names = 1)
# keep only the topics
topic_ann <- topic_ann[, grep('Topic', colnames(topic_ann))]
# keep only cells we have topic annotation info for
signac_celltype <- signac_celltype[, colnames(signac_celltype) %in% rownames(topic_ann)]
# add the topic annotation data
signac_celltype <- AddMetaData(signac_celltype, topic_ann)

# these are the topics we'll consider
topics <- colnames(topic_ann)
# we'll do only a specific topic if requested
# if (!is.null(opt[['topic']])) {
  # topics <- c(opt[['topic']])
# }
topics <- paste('Topic', c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20), sep = '')


# the location of the input matrices
imputed_matrix_loc <- paste('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/imputed_pycistopic_matrices/', cell_type, '/', sep = '')
# the prepend
imputed_matrix_prepend <- 'matrix_'
# the append
imputed_matrix_append <- '.mtx.gz'

# list the files in the directory
matrix_files <- list.files(imputed_matrix_loc, pattern = paste(imputed_matrix_prepend, '(\\d+)_(\\d+)', imputed_matrix_append, sep = ''))

# get the location of the barcodes and features files
features_loc <- paste(imputed_matrix_loc, 'regiondata.tsv.gz', sep = '')
barcodes_loc <- paste(imputed_matrix_loc, 'barcodes.tsv.gz', sep = '')

# read the features file
features <- read.table(features_loc, header = F, sep = '\t')
# read the barcodes
barcodes <- read.table(barcodes_loc, header = F)$V1

# get the feature names
feature_names <- features$V4

# go through each matrix
for (matrix_file in matrix_files) {
  # extract the region
  regions_string <- stringr::str_extract(matrix_file, '(\\d+)_(\\d+)')
  # split by underscore
  regions_vector <- regions_string[[1]]
  # read the accompanying features file
  features_matrix_loc <- paste(imputed_matrix_loc, 'features_', regions_string, '.tsv.gz', sep = '')
  features_matrix <- read.table(features_matrix_loc)$V1
  # read the matrix
  matrix_regions <- Matrix::readMM(paste(imputed_matrix_loc, matrix_file, sep = ''))
  # set the features and barcodes
  colnames(matrix_regions) <- barcodes
  rownames(matrix_regions) <- features_matrix
  # extract the metadata
  signac_metadata <- signac_celltype@meta.data
  signac_fragments <- Fragments(signac_celltype)
  # create the chromatin assay
  chrom_assay <- CreateChromatinAssay(
    counts = matrix_regions,
    sep = c(":", "-"),
    fragments = signac_fragments,
    min.cells = 10,
    min.features = 200
  )
  # create object
  seurat_object_regions <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = signac_metadata[barcodes, ],
    project = 'wijst_multiome'
  )
  # set the annotations to the object now
  #Annotation(seurat_object) <- annotations
  # create the output directory
  output_dir_full <- paste(limma_output_loc, '/', cell_type, '/', regions_string, '/', sep = '')
  dir.create(output_dir_full, recursive = T, showWarnings = F)
  # do the bulk analysis
  for(topic in topics){
    condition_combinations <- list()
    condition_combinations[[topic]] <- c('True', 'False')
    do_limma_dream_pairwise_per_celltype(seurat_object_regions, 
                                         output_loc = output_dir_full, 
                                         condition_combinations = condition_combinations,
                                         celltype_column = celltype_column, 
                                         aggregates = c(topic, 'lane', 'sample_final'), 
                                         fixed_effects = c(topic, 'age', 'sex'), 
                                         random_effects = c('sample_final', 'lane'),
                                         minimal_cells = min_cells,
                                         min_peaks = min_cell_umis, 
                                         minimal_complexity = min_pseudo_umis, 
                                         nthreads = nthreads, 
                                         permute = permute)
  }
}

