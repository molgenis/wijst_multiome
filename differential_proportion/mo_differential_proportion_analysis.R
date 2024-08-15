############################################################################################################################
# Authors: Roy Oelen
# Name: mo_differential_proportion_analysis.R
# Function: perform differential proportion analysis with speckle
############################################################################################################################

####################
# libraries        #
####################

library(speckle) #BiocManager::install(c("CellBench", "BiocStyle", "scater", "org.Mm.eg.db"))
library(SingleCellExperiment)
library(CellBench)
library(limma)
library(ggplot2)
library(scater)
library(patchwork)
library(edgeR)
library(statmod)

####################
# Functions        #
####################

#' perform differential proportion analysis
#' 
#' @param seurat_object The Seurat object to perform differential proportion on, supply this or straight-up metadata
#' @param metadata The metadata to perform differential proportion on, supply this or the Seurat object
#' @param cell_type_column the column in the metadata describing the celltype
#' @param variable_of_interest the variable we are interested in for DPA
#' @param participant_column the column denoting the participant the cells came from
#' @param covariates vector of covariates to correct for
#' @returns the table of the fit for the coefficient of interest
#' fit_ct <- do_differential_proportion_analysis(mo, cell_type_column = 'cell_type_safe', covariates=c('sex', 'age', 'day'))
do_differential_proportion_analysis <- function(seurat_object=NULL, metadata=NULL, cell_type_column='cell_type_final', variable_of_interest='inflammation_status', participant_column='biopsy_id', covariates=c('day')) {
  # initialize variable
  meta.data <- NULL
  # extract metadata from object if required
  if (!is.null(metadata) & is.null(seurat_object)) {
    meta.data <- metadata
  }
  else if (is.null(metadata) & is.null(seurat_object)) {
    stop('metadata and Seurat object are both NULL, either needs to be non-NULL')
  }
  else if (!is.null(metadata) & !is.null(seurat_object)) {
    warning('both metadata and Seurat object are non-NULL, will use metadata')
    meta.data <- metadata
  }
  else if (is.null(metadata) & !is.null(seurat_object)) {
    meta.data <- seurat_object@meta.data
  }
  # collect the variables we need to combine
  combined_variables <- c(variable_of_interest, participant_column, covariates)
  # add as a column
  for (variable in combined_variables) {
    # as a new one if this is the first variable
    if (! ('combined_variable' %in% colnames(meta.data)) ) {
      meta.data[['combined_variable']] <- meta.data[[variable]]
    }
    # or add if we were already building
    else{
      meta.data[['combined_variable']] <- paste(meta.data[['combined_variable']], meta.data[[variable]], sep = '_')
    }
  }
  # convert to proportions
  props <- getTransformedProps(clusters = meta.data[[cell_type_column]], sample = meta.data[['combined_variable']])
  # get the unique combinations of covariates, variable of interest, and participant column
  pair_group <- unique(meta.data[, c('combined_variable', combined_variables)])
  # sort the pair group by this combination
  pair_group <- pair_group[order(pair_group[['combined_variable']]), ]
  # extract the 'pair', which is the donor
  pair <- pair_group[[participant_column]]
  # extract the 'group', which is the condition to compare
  group <- pair_group[[variable_of_interest]]
  # create the formula
  formula_to_use <- paste('~', variable_of_interest, '+', paste(covariates, collapse='+'), sep = '')
  des.tech <- model.matrix(as.formula(formula_to_use), data=pair_group)
  # create the correlations, with replicates (the donors)
  dupcor <- duplicateCorrelation(props$TransformedProps, design=des.tech,
                                 block=pair)
  # fit a regression model using limma
  fit1 <- lmFit(props$TransformedProps, design=des.tech, block=pair, 
                correlation=dupcor$consensus)
  # calculate statistics on the model
  fit1 <- eBayes(fit1)
  # get statistically significant proportions
  fit1_significant <- decideTests(fit1)
  # show summary
  summary(fit1_significant)
  # and specifically which ones
  topTable(fit1,coef=2, number = length(unique(meta.data[[cell_type_column]])))
  return(fit1)
}

####################
# Main Code        #
####################

# we can also just read the metadata itself
mo_metadata <- read.table('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_celllevel_metadata_20240620.tsv.gz', sep = '\t', header = T)

# filter for empty rows
mo_metadata <- mo_metadata[!is.na(mo_metadata[['age']]) &
                             !is.na(mo_metadata[['sex']]) &
                             !is.na(mo_metadata[['lane']]) &
                             !is.na(mo_metadata[['condition_final']]) &
                             !is.na(mo_metadata[['sample_final']]), ]

# location of results
result_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/differential_proportion/speckle/results/'

# perform differential proportion analysis using the high resolution cell types
fit_ct <- do_differential_proportion_analysis(metadata = mo_metadata[!is.na(mo_metadata[['predicted.mo_10x_cell_type']]), ], 
                                              cell_type_column = 'predicted.mo_10x_cell_type', 
                                              covariates=c('sex', 'age', 'lane'), 
                                              variable_of_interest = 'condition_final', 
                                              participant_column = 'sample_final')

# get statistically significant proportions
fit_ct_significant <- decideTests(fit_ct)
# show summary
summary(fit_ct_significant)
# and specifically which ones
topTable(fit_ct,coef=2, number = length(unique(mo_metadata[!is.na(mo_metadata[['predicted.mo_10x_cell_type']]),'predicted.mo_10x_cell_type'])))

# save result
fit_ct_df <- data.frame(topTable(fit_ct, coef=2, number = length(unique(mo_metadata[!is.na(mo_metadata[['predicted.mo_10x_cell_type']]),'predicted.mo_10x_cell_type']))))
fit_ct_df <- cbind(data.frame(celltype=rownames(fit_ct_df)), fit_ct_df)
write.table(fit_ct_df,
            paste(result_loc, 'high_condition_20240620.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)

# perform differential proportion analysis using the L1 resolution cell types
fit_mhct <- do_differential_proportion_analysis(metadata = mo_metadata[!is.na(mo_metadata[['celltype_imputed_lowerres']]), ], 
                                                cell_type_column = 'celltype_imputed_lowerres', covariates=c('sex', 'age', 'lane'), 
                                                variable_of_interest = 'condition_final', 
                                                participant_column = 'sample_final')

# get statistically significant proportions
fit_mhct_significant <- decideTests(fit_mhct)
# show summary
summary(fit_mhct_significant)
# and specifically which ones
topTable(fit_mhct, coef=2, number = length(unique(mo_metadata[['celltype_imputed_lowerres']])))

# save result
fit_mhct_df <- data.frame(topTable(fit_mhct, coef=2, number = length(unique(mo_metadata[['celltype_imputed_lowerres']]))))
fit_mhct_df <- cbind(data.frame(celltype=rownames(fit_mhct_df)), fit_mhct)
write.table(fit_mhct_df,
            paste(result_loc, 'lowerres_condition_20240620.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)

# perform differential proportion analysis using the high resolution cell types
fit_ct_covid <- do_differential_proportion_analysis(metadata = mo_metadata[!is.na(mo_metadata[['predicted.mo_10x_cell_type']]), ], 
                                              cell_type_column = 'predicted.mo_10x_cell_type', 
                                              covariates=c('sex', 'age', 'lane'), 
                                              variable_of_interest = 'LONG_COVID_final', 
                                              participant_column = 'sample_final')

# get statistically significant proportions
fit_ct_covid_significant <- decideTests(fit_ct_covid)
# show summary
summary(fit_ct_covid_significant)
# and specifically which ones
topTable(fit_ct_covid,coef=2, number = length(unique(mo_metadata[!is.na(mo_metadata[['predicted.mo_10x_cell_type']]),'predicted.mo_10x_cell_type'])))

# save result
fit_ct_covid_df <- data.frame(topTable(fit_ct_covid, coef=2, number = length(unique(mo_metadata[!is.na(mo_metadata[['predicted.mo_10x_cell_type']]),'predicted.mo_10x_cell_type']))))
fit_ct_covid_df <- cbind(data.frame(celltype=rownames(fit_ct_covid_df)), fit_ct_covid_df)
write.table(fit_ct_covid_df,
            paste(result_loc, 'high_covid_20240620.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)


# perform differential proportion analysis using the L1 resolution cell types
fit_mhct_covid <- do_differential_proportion_analysis(metadata = mo_metadata[!is.na(mo_metadata[['celltype_imputed_lowerres']]), ], 
                                                cell_type_column = 'celltype_imputed_lowerres', covariates=c('sex', 'age', 'lane'), 
                                                variable_of_interest = 'LONG_COVID_final', 
                                                participant_column = 'sample_final')

# get statistically significant proportions
fit_mhct_covid_significant <- decideTests(fit_mhct_covid)
# show summary
summary(fit_mhct_covid_significant)
# and specifically which ones
topTable(fit_mhct_covid, coef=2, number = length(unique(mo_metadata[['celltype_imputed_lowerres']])))

# save result
fit_mhct_covid_df <- data.frame(topTable(fit_mhct_covid, coef=2, number = length(unique(mo_metadata[['celltype_imputed_lowerres']]))))
fit_mhct_covid_df <- cbind(data.frame(celltype=rownames(fit_mhct_covid_df)), fit_mhct_covid)
write.table(fit_mhct_covid_df,
            paste(result_loc, 'lowerres_covid_20240620.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)
