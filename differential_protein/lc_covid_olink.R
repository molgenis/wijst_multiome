#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: lc_covid_olink.R
# Function: run differential protein expression on the olink data
############################################################################################################################

####################
# libraries        #
####################

# Load OlinkAnalyze
library(OlinkAnalyze)
# for plotting
library(ggplot2)
library(cowplot)
# we'll do correlations in parallel
library(foreach)
library(doParallel)
# limma DE dependencies
library(variancePartition)
library(edgeR)
library(BiocParallel)
# manual regression
library(lme4)
library(lmerTest)

####################
# Functions        #
####################

plot_plate_layout <- function(olink_data, well_id_column='WellID', plate_id_column='PlateID', sample_column='SampleID') {
  # add the row
  olink_data[['row']] <- gsub('\\d+', '', olink_data[[well_id_column]])
  olink_data[['column']] <- gsub('[A-Z]+', '', olink_data[[well_id_column]])
  # order the rows in reverse order
  olink_data[['row']] <- factor(olink_data[['row']], levels = unique(olink_data[['row']])[rev(order(unique(olink_data[['row']])))])
  # also order the columns 
  olink_data[['column']] <- factor(olink_data[['column']], levels = as.character(1:max(as.numeric(olink_data[['column']]))))
  # subset to the data we can plot
  olink_data_samples <- olink_data[, c(plate_id_column, sample_column, 'row', 'column')]
  # make a tile plot
  p <- ggplot(data=olink_data_samples, aes(x=column, y=row, fill = .data[[plate_id_column]])) + 
    geom_tile() + geom_text(data = olink_data, mapping = aes(label=.data[[sample_column]])) +
    facet_grid(rows = vars(olink_data[[plate_id_column]]))
  # some options
  return(p)
}

correlate_samples <- function(olink_data, sample_column='SampleID', protein_name_column='OlinkID', value_column='NPX', cor_method = 'spearman') {
  # get the unique samples
  samples_present <- unique(olink_data[[sample_column]])
  # remove empty ones
  samples_present <- samples_present[!is.na(samples_present)]
  # get the number of samples
  n_samples <- length(samples_present)
  # check each row in the square matrix
  res_per_sample_on_row <- foreach(i = 1:n_samples, .combine = 'rbind') %dopar% {
    # this is how many correlations we'll calculate
    n_cors <- n_samples - i + 1
    # reserve memory
    res_frame <- data.frame(matrix(NA, nrow = n_cors, ncol = 6))
    # set the column names
    colnames(res_frame) <- c('sample1', 'sample2', 'nproteins_s1', 'nproteins_s2', 'nproteins_both', 'correlation')
    # we need to keep an index to know what values to fill in for the res_frame
    res_frame_i <- 1
    # correlate with other samples
    for (i2 in i:n_samples) {
        # extract sample 1 and sample 2
        sample_1 <- samples_present[i]
        sample_2 <- samples_present[i2]
        # subset to values for these two samples
        sample_1_data <- olink_data[
          !is.na(olink_data[[sample_column]]) &
            olink_data[[sample_column]] == sample_1, 
            c(protein_name_column, value_column), 
        ]
        sample_2_data <- olink_data[
          !is.na(olink_data[[sample_column]]) &
            olink_data[[sample_column]] == sample_2, 
            c(protein_name_column, value_column)
        ]
        # let's check how many proteins we have for each of these
        proteins_sample_1_wna <- unique(sample_1_data[[protein_name_column]])
        proteins_sample_2_wna <- unique(sample_2_data[[protein_name_column]])
        # drop NA as well
        proteins_sample_1 <- proteins_sample_1_wna[!is.na(sample_1_data[match(proteins_sample_1_wna, sample_1_data[[protein_name_column]]), value_column])]
        proteins_sample_2 <- proteins_sample_2_wna[!is.na(sample_2_data[match(proteins_sample_2_wna, sample_2_data[[protein_name_column]]), value_column])]
        # warn if dropping due to NA
        if (length(proteins_sample_1) < length(proteins_sample_1_wna)) {
          warning(paste('dropping proteins when from', sample_1, 'due to NA values:', paste(setdiff(proteins_sample_1_wna, proteins_sample_1)), '\n'))
        }
        if (length(proteins_sample_2) < length(proteins_sample_2_wna)) {
          warning(paste('dropping proteins when from', sample_2, 'due to NA values:', paste(setdiff(proteins_sample_2_wna, proteins_sample_2)), '\n'))
        }
        # let's check the inner join of those two sets of proteins
        proteins_sample_both <- intersect(proteins_sample_1, proteins_sample_2)
        # if the sizes are different, we will drop proteins, let them know which ones
        if (length(proteins_sample_both) < length(proteins_sample_1)) {
          warning(paste('dropping proteins when correlating', sample_1, 'and', sample_2, 'from', sample_1, ':', paste(setdiff(proteins_sample_1, proteins_sample_both)), '\n'))
        }
        if (length(proteins_sample_both) < length(proteins_sample_2)) {
          warning(paste('dropping proteins when correlating', sample_1, 'and', sample_2, 'from', sample_2, ':', paste(setdiff(proteins_sample_2, proteins_sample_both)), '\n'))
        }
        # subset to those proteins
        prot_values_sample1 <- sample_1_data[match(proteins_sample_both, sample_1_data[[protein_name_column]]), value_column]
        prot_values_sample2 <- sample_2_data[match(proteins_sample_both, sample_2_data[[protein_name_column]]), value_column]
        # calculate the correlation
        cor_both <- cor(prot_values_sample1, prot_values_sample2, method = cor_method, use = 'complete.obs')
        # make result into table c('sample1', 'sample2', 'nproteins_s1', 'nproteins_s2', 'nproteins_both', 'correlation')
        res <- c(sample_1, sample_2, length(proteins_sample_1), length(proteins_sample_2), length(proteins_sample_both), cor_both)
        res_frame[res_frame_i, ] <- res
        # and update the result frame index
        res_frame_i <- res_frame_i + 1
    }
    return(res_frame)
  }
  # make correlation numeric again
  res_per_sample_on_row[['correlation']] <- as.numeric(res_per_sample_on_row[['correlation']])
  return(res_per_sample_on_row)
}


create_confusion_matrix <- function(assignment_table, truth_column, prediction_column, truth_column_label=NULL, prediction_column_label=NULL, angle_labels=T, confusion_table=NULL, freq_column='freq', show_text=T, axis_text_x_size=3, axis_text_y_size=3){
  # create confusion table from assignments
  if (is.null(confusion_table)) {
    confusion_table <- create_confusion_table(assignment_table, truth_column, prediction_column)
  }
  # unless we already have the confusion table
  else {
    # then we just need to harmonize the value
    confusion_table$freq <- confusion_table[[freq_column]]
    confusion_table$truth <- confusion_table[[truth_column]]
    confusion_table$prediction <- confusion_table[[prediction_column]]
  }
  # round the frequency off to a sensible cutoff
  confusion_table$freq <- round(confusion_table$freq, digits=2)
  # turn into plot
  p <- ggplot(data=confusion_table, aes(x=truth, y=prediction, fill=freq)) + geom_tile() + scale_fill_gradient2(low='red', high='blue', mid = 'white', limits = c(-1, 1))
  # add text if requested
  if (show_text) {
    p <- p + geom_text(aes(label=freq))
  }
  # some options
  if(!is.null(truth_column_label)){
    p <- p + xlab(truth_column_label)
  }
  if(!is.null(prediction_column_label)){
    p <- p + ylab(prediction_column_label)
  }
  if (angle_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  # set text sizes
  p <- p + theme(axis.text.x = element_text(color = "grey20", size = axis_text_x_size, face = "plain"),
        axis.text.y = element_text(size = axis_text_y_size, face = "plain"))
  return(p)
}


olink_to_limma_format <- function(olink_data, metadata_columns_of_interest, sample_column='SampleID', olink_protein_name_column='OlinkID', olink_protein_value_column='NPX') {
  # we'll make a standard dataframe from this
  olink_data <- data.frame(olink_data)
  # extract the metadata we care about
  metadata_of_interest <- unique(olink_data[, unique(c(metadata_columns_of_interest, sample_column))])
  # get the samplew
  samples <- unique(olink_data[[sample_column]])
  # get the proteins
  proteins <- unique(olink_data[[olink_protein_name_column]])
  # create a matrix
  expression_matrix <- matrix(NA, nrow = length(proteins), ncol = length(samples), dimnames = list(proteins, samples))
  # start filling the matrix
  for (i in 1:nrow(olink_data)) {
    # add each entry
    expression_matrix[olink_data[i, olink_protein_name_column], olink_data[i, sample_column]] <- olink_data[i, olink_protein_value_column]
  }
  # make sure they are in the same order
  metadata_of_interest <- metadata_of_interest[match(colnames(expression_matrix), metadata_of_interest[[sample_column]]), ]
  # save the these in a list
  return(list('expression' = expression_matrix, 'metadata' = metadata_of_interest))
}


plot_olink_expression <- function(olink, protein_name, protein_name_column='OlinkID', protein_value_column='NPX', group_column='case_control', olinkid_to_uid_loc=NULL, uniprotid_to_gs_loc=NULL, violin=F, paper_style=T, legendless=T, pointless=F, use_label_dict=F, plot_order=NULL){
  # subset to the protein we care about
  olink <- olink[olink[[protein_name_column]] == protein_name, ]
  # add extra column
  olink$protein_expression <- olink[[protein_value_column]]
  # order the x
  if (!is.null(plot_order)) {
    olink[[group_column]] <- factor(olink[[group_column]], levels = plot_order)
  }
  # create the plot
  p <- NULL
  if(violin){
    p <- ggplot(data=olink, mapping=aes(x=.data[[group_column]], y=protein_expression, fill=.data[[group_column]])) +
      geom_violin() +
      geom_jitter(size = 0.5, alpha = 0.5)
  }
  else{
    p <- ggplot(data=olink, mapping=aes(x=.data[[group_column]], y=protein_expression, fill=.data[[group_column]])) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 0.5, alpha = 0.5)
  }
  # create title
  title <- paste('protein expression of', protein_name)
  if(!is.null(olinkid_to_uid_loc)){
    # read mapping of olink ID to uniprot ID
    olinkid_to_uid <- read.table(olinkid_to_uid_loc, sep = '\t', header = T, stringsAsFactors = F)
    # get the uid
    uid <- olinkid_to_uid[olinkid_to_uid$OlinkID == protein_name, 'Uniprot.ID']
    # add to title
    title <- paste(title, '-', uid, sep = '')
    # if there is a gene symbol mapping, do that one as well
    if(!is.null(uniprotid_to_gs_loc)){
      uniprotid_to_gs <- read.table(uniprotid_to_gs_loc, sep = '\t', header = T, stringsAsFactors = F)
      gs <- uniprotid_to_gs[uniprotid_to_gs$From == uid, 'To']
      # add to title
      title <- paste(title, '(', gs, ')')
    }
  }
  # add paper style if requested
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  # remove points if requested
  if (pointless) {
    p <- p + theme(axis.text.x=element_blank(),
                   axis.ticks = element_blank(),
                   axis.title.x = element_blank())
  }
  # remove legend if requested
  if (legendless) {
    p <- p + theme(legend.position = 'none')
  }
  # add title
  p <- p + ggtitle(title)
  # and y lab
  p <- p + ylab('protein expression')
  return(p)
}


do_inverse_normal <- function(protein_data, protein_name_column='OlinkID', protein_value_column='NPX', new_protein_value_column_name='INE') {
  # add this new normalized value as a column
  protein_data[[new_protein_value_column_name]] <- NA
  # check each protein
  for (protein_name in unique(protein_data[[protein_name_column]])) {
    # get indices of this protein
    protein_indices <- which(!is.na(protein_data[[protein_name_column]]) & protein_data[[protein_name_column]] == protein_name)
    # extract the values for that
    protein_expression_raw <- protein_data[protein_indices, protein_value_column]
    # normalize these
    protein_expression_norm <- qnorm((rank(protein_expression_raw, na.last = "keep") -0.5) / sum(!is.na(protein_expression_raw)))
    # and add these values
    protein_data[protein_indices, new_protein_value_column_name] <- protein_expression_norm
  }
  return(protein_data)
}


do_regression <- function(protein_data, protein_name_column='OlinkID', protein_value_column='NPX', fixed_effects=c('age', 'sex', 'pandemic', 'case_control'), random_effects=c('SampleID'), value_of_interest='case_control') {
  # paste together the model
  model_formula <- paste(protein_value_column, '~ 0')
  for(fixed_effect in fixed_effects){
    model_formula <- paste(model_formula, fixed_effect, sep = ' + ')
  }
  for(random_effect in random_effects){
    model_formula <- paste(model_formula, ' + (1|', random_effect, ')', sep = '')
  }
  # and turn into a formula
  form <- as.formula(model_formula)
  
  # check each protein
  res_per_protein <- list()
  for (protein in unique(protein_data[[protein_name_column]])) {
    # subset to that protein
    protein_data_protein <- protein_data[!is.na(protein_data[[protein_name_column]]) & protein_data[[protein_name_column]] == protein, ]
    # subset to complete data
    protein_data_protein <- protein_data_protein[complete.cases(protein_data_protein[, c(fixed_effects, random_effects)]), ]
    tryCatch({
      if (nrow(protein_data_protein) > 0) {
        # do analysis
        if (length(random_effects) > 0) {
          res_per_protein[[protein]] <- lmer(form, data = protein_data_protein)
        }
        else {
          res_per_protein[[protein]] <- glm(form, data = protein_data_protein)
        }
      }
    }, error=function(cond) {
      print(paste('model build failed'))
      message(cond)
    })
  }
  
  # convert to dataframe
  df_per_protein <- list()
  for (protein in names(res_per_protein)) {
    # extract the result
    res_protein <- res_per_protein[[protein]]
    # turn into dataframe
    protein_df <- data.frame(summary(res_protein)$coefficients)
    # add the rownames as explicit column
    protein_df <- cbind(data.frame(term = rownames(protein_df)), protein_df)
    # add the protein
    protein_df <- cbind(data.frame(protein = rep(protein, times = nrow(protein_df))), protein_df)
    # put in list
    df_per_protein[[protein]] <- protein_df
  }
  # save in big table
  df_proteins <- do.call('rbind', df_per_protein)
  return(df_proteins)
}


####################
# Main Code        #
####################

# location of the protein data

protein_data_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/olink_protein/OLINK2023-037_INF_EXTENDED_NPX_2024-07-10.csv'
# location of the sample mapping
sample_mapping_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/olink_protein/OV21_00402_linkage_file_olink_20240719.tsv.gz'
# the full id table
mo_full_id_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_full_id_table.tsv.gz'
# the age/sex for the mo participant
mo_age_sex_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_age_sex.tsv.gz'
# the mapping of the covid assignments
lc_assignments_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/olink_protein/mo_case_dob_sex.tsv'
# the next info
mo_sample_sheet_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_age_sex_ids.tsv.gz'

# load the data
protein_data <- read_NPX(protein_data_loc)
# load the sample data
sample_mapping <- read.table(sample_mapping_loc, header = T, sep = '\t')

# load sample sheet
mo_sample_sheet <- read.table(mo_sample_sheet_loc, header = T, sep = '\t')

# load the sample data
sample_mapping <- read.table(sample_mapping_loc, header = T, sep = '\t')
# make sex lowercase
sample_mapping[['Gender']] <- tolower(sample_mapping[['Gender']])

# add assignment data
lc_assignments <- read.table(lc_assignments_loc, header = T, sep = '\t')
# making sex lowercase
lc_assignments[['sex']] <- tolower(lc_assignments[['sex']])
# and add MO to the name
lc_assignments[['olink_id']] <- paste('MO', lc_assignments[['olink_id']], sep = '')
# add the MO olink identifier to hte sample mapping file
sample_mapping[['olink_id']] <- lc_assignments[match(paste(sample_mapping[['dob']], sample_mapping[['Gender']]), paste(lc_assignments[['date_of_birth']], lc_assignments[['sex']])), 'olink_id']

# load age/sex data
mo_age_sex <- read.table(mo_age_sex_loc, header = T, sep = '\t')

# rename the plate
sample_mapping[['plate']] <- gsub('OV210402EDPL0', 'PlateLayout_Plate', sample_mapping$Release.Rack)
# also plate1 and plate2 are actually plate 4 and 5
sample_mapping[sample_mapping[['plate']] == 'PlateLayout_Plate1', 'plate'] <- 'PlateLayout_Plate4'
sample_mapping[sample_mapping[['plate']] == 'PlateLayout_Plate2', 'plate'] <- 'PlateLayout_Plate5'
# remove the 0 in the well name
sample_mapping[['Release.Position']] <- gsub('(0)(\\d+)$', '\\2', sample_mapping[['Release.Position']])
# add combination of well and plate
sample_mapping[['well_plate']] <- paste(sample_mapping[['plate']], sample_mapping[['Release.Position']], sep =  '_')

# the entries with just a number are MO ones, rename them
protein_data[['SampleID']] <- gsub('^(\\d+)$', 'MO\\1', protein_data[['SampleID']])
# also the duplicated one
protein_data[['SampleID']] <- gsub('(82_1)|(82_2)', 'MO82', protein_data[['SampleID']])
# get more descriptive ID for those
#protein_data_realid <- mo_realid_assignments[match(protein_data[['SampleID']], mo_realid_assignments[['sample_final']]), 'realid_final']
#protein_data_realid <- lc_assignments[match(protein_data[['SampleID']], lc_assignments[['olink']]), 'realid']
# replace the sample name where possible
#protein_data[!is.na(protein_data_realid), 'realid'] <- protein_data_realid[!is.na(protein_data_realid)]
# if the ID already starts with ll-next, we can just copy that
protein_data[!is.na(protein_data[['SampleID']]) & startsWith(protein_data[['SampleID']], 'LL-NEXT'), 'realid'] <- protein_data[!is.na(protein_data[['SampleID']]) & startsWith(protein_data[['SampleID']], 'LL-NEXT'), 'SampleID']

# get the samples from lifelines
protein_data_sex2 <- sample_mapping[match(protein_data[['SampleID']], sample_mapping[['SampleID']]), 'Gender']
protein_data_age2 <- sample_mapping[match(protein_data[['SampleID']], sample_mapping[['SampleID']]), 'Age']
# again add sex and age
protein_data[!is.na(protein_data_sex2), 'sex'] <- protein_data_sex2[!is.na(protein_data_sex2)]
protein_data[!is.na(protein_data_age2), 'age'] <- protein_data_age2[!is.na(protein_data_age2)]
# and the mo ones
protein_data_sex4 <- mo_age_sex[match(protein_data[['SampleID']], mo_age_sex[['sample']]), 'sex']
protein_data_age4 <- mo_age_sex[match(protein_data[['SampleID']], mo_age_sex[['sample']]), 'age']
protein_data[!is.na(protein_data_sex4), 'sex'] <- protein_data_sex4[!is.na(protein_data_sex4)]
protein_data[!is.na(protein_data_age4), 'age'] <- protein_data_age4[!is.na(protein_data_age4)]

# everything that is not MO, is a control
protein_data[!(grepl('(MO\\d+)', protein_data[['SampleID']])), 'case_control'] <- 'control'
# set pandemic status as well
protein_data[protein_data[['PlateID']] %in% c('PlateLayout_Plate1', 'PlateLayout_Plate2', 'PlateLayout_Plate3'), 'pandemic'] <- 'postpandemic'
protein_data[protein_data[['PlateID']] %in% c('PlateLayout_Plate4', 'PlateLayout_Plate5'), 'pandemic'] <- 'prepandemic'
# the NEXT samples on the other plates are still post-pandemic
protein_data[(grepl('(LL-NEXT)', protein_data[['SampleID']])), 'pandemic'] <- 'postpandemic'

# get combination of well and plate
protein_data_well_plate <- paste(protein_data[['PlateID']], protein_data[['WellID']], sep = '_')

# get the status where we can
protein_data_case <- lc_assignments[match(protein_data[['SampleID']], lc_assignments[['olink_id']]), 'case_controle']
# and add that information
protein_data[!is.na(protein_data_case), 'case_control'] <- protein_data_case[!is.na(protein_data_case)]
# add the combination of case/control
protein_data[['casecontrolpandemic']] <- NA
protein_data[protein_data[['pandemic']] == 'postpandemic' & protein_data[['case_control']] == 'case', 'casecontrolpandemic'] <- 'case'
protein_data[protein_data[['pandemic']] == 'postpandemic' & protein_data[['case_control']] == 'control', 'casecontrolpandemic'] <- 'control-recent'
protein_data[protein_data[['pandemic']] == 'prepandemic' & protein_data[['case_control']] == 'control', 'casecontrolpandemic'] <- 'control-old'

# add the final sample id
protein_data[['sample_final']] <- protein_data[['SampleID']]
# add inferred olink id
protein_data[['olink_id']] <- sample_mapping[!is.na(sample_mapping[['olink_id']]) & sample_mapping[['olink_id']] != 'MO' & sample_mapping[['olink_id']] != '', ][match(protein_data[['SampleID']], sample_mapping[!is.na(sample_mapping[['olink_id']]) & sample_mapping[['olink_id']] != 'MO' & sample_mapping[['olink_id']] != '', 'SampleID']), 'olink_id']
# use that instead of the sample id where possible
protein_data[!is.na(protein_data[['olink_id']]), 'sample_final'] <- protein_data[!is.na(protein_data[['olink_id']]), 'olink_id']

# plot the plates to see how they line up
plot_plate_layout(protein_data) + scale_fill_manual(values = roycols::get_color_list(protein_data[['PlateID']]))
plot_plate_layout(protein_data, sample_column = 'case_control') + scale_fill_manual(values = roycols::get_color_list(protein_data[['PlateID']]))
plot_plate_layout(protein_data, sample_column = 'pandemic') + scale_fill_manual(values = roycols::get_color_list(protein_data[['PlateID']]))
plot_plate_layout(protein_data, sample_column = 'sex') + scale_fill_manual(values = roycols::get_color_list(protein_data[['PlateID']]))
plot_plate_layout(protein_data, sample_column = 'age') + scale_fill_manual(values = roycols::get_color_list(protein_data[['PlateID']]))

# we may have matched samples in some cases, so we'll combine some columns to make those samples unique
protein_data[['sample_unique']] <- paste(protein_data[['sample_final']], protein_data[['PlateID']], protein_data[['WellID']])
# correlate the samples to one another
protein_correlations_samples <- correlate_samples(protein_data, sample_column = 'sample_unique')
# also add the reverse sample comparison (we have A-B, but we also want B-A there explicitly)
protein_correlations_samples <- rbind(protein_correlations_samples,
                                      data.frame(sample1 = protein_correlations_samples[['sample2']],
                                                 sample2 = protein_correlations_samples[['sample1']],
                                                 nproteins_s1 = protein_correlations_samples[['nproteins_s2']], 
                                                 nproteins_s2 = protein_correlations_samples[['nproteins_s1']],
                                                 nproteins_both = protein_correlations_samples[['nproteins_both']], 
                                                 correlation = protein_correlations_samples[['correlation']]))
# unfortunately that duplicated the A-A ones as well, so let's remove the duplicates
protein_correlations_samples <- protein_correlations_samples[!duplicated(paste(protein_correlations_samples[['sample1']], protein_correlations_samples[['sample2']])), ]
# let's plot a tile, to see if there is any pattern at all
create_confusion_matrix(assignment_table = NULL, confusion_table = protein_correlations_samples, freq_column = 'correlation', truth_column = 'sample1', prediction_column = 'sample2', truth_column_label='sample', prediction_column_label='sample', angle_labels=T, show_text = F)
# also plot specifically the different controls agains one another
plot_grid(
  create_confusion_matrix(assignment_table = NULL, confusion_table = protein_correlations_samples[protein_correlations_samples[['sample1']] %in% unique(data.frame(protein_data)[data.frame(protein_data)[['casecontrolpandemic']] == 'control-old', 'sample_unique']) &
                                                                                                    protein_correlations_samples[['sample2']] %in% unique(data.frame(protein_data)[data.frame(protein_data)[['casecontrolpandemic']] == 'control-recent', 'sample_unique']), ], freq_column = 'correlation', truth_column = 'sample1', prediction_column = 'sample2', truth_column_label='sample', prediction_column_label='sample', angle_labels=T, show_text = F) +
    ggtitle('control-old vs control-recent'),
  create_confusion_matrix(assignment_table = NULL, confusion_table = protein_correlations_samples[protein_correlations_samples[['sample1']] %in% unique(data.frame(protein_data)[data.frame(protein_data)[['casecontrolpandemic']] == 'case', 'sample_unique']) &
                                                                                                    protein_correlations_samples[['sample2']] %in% unique(data.frame(protein_data)[data.frame(protein_data)[['casecontrolpandemic']] == 'control-recent', 'sample_unique']), ], freq_column = 'correlation', truth_column = 'sample1', prediction_column = 'sample2', truth_column_label='sample', prediction_column_label='sample', angle_labels=T, show_text = F) + 
    ggtitle('case vs control-recent'),
  create_confusion_matrix(assignment_table = NULL, confusion_table = protein_correlations_samples[protein_correlations_samples[['sample1']] %in% unique(data.frame(protein_data)[data.frame(protein_data)[['casecontrolpandemic']] == 'case', 'sample_unique']) &
                                                                                                    protein_correlations_samples[['sample2']] %in% unique(data.frame(protein_data)[data.frame(protein_data)[['casecontrolpandemic']] == 'control-old', 'sample_unique']), ], freq_column = 'correlation', truth_column = 'sample1', prediction_column = 'sample2', truth_column_label='sample', prediction_column_label='sample', angle_labels=T, show_text = F) + 
    ggtitle('case vs control-old'),
  nrow = 2,
  ncol = 2
  
)


# add inverse normal transformation
protein_data <- do_inverse_normal(protein_data)
# save in big table
df_proteins <- do_regression(protein_data, random_effects=c('sample_final', 'PlateID'))
# subset to the variate we care about
df_proteins_casecontrol <- df_proteins[df_proteins[['term']] == 'case_controlcontrol', ]
# do B&H correction
df_proteins_casecontrol[['BH']] <- p.adjust(df_proteins_casecontrol[['Pr...t..']], method = 'BH')
df_proteins_casecontrol[['p.bonferroni']] <- p.adjust(df_proteins_casecontrol[['Pr...t..']], method = 'bonferroni')
df_proteins_casecontrol <- df_proteins_casecontrol[order(df_proteins_casecontrol[['Pr...t..']]), ]
# add Z-score
df_proteins_casecontrol[['z']] <- df_proteins_casecontrol[['Estimate']] / df_proteins_casecontrol[['Std..Error']]

# try inverse normal as well
df_proteins_ine <- do_regression(protein_data, random_effects=c('sample_final', 'PlateID'), protein_value_column = 'INE')
df_proteins_ine_casecontrol <- df_proteins_ine[df_proteins_ine[['term']] == 'case_controlcontrol', ]
df_proteins_ine_casecontrol[['BH']] <- p.adjust(df_proteins_ine_casecontrol[['Pr...t..']], method = 'BH')
df_proteins_ine_casecontrol[['p.bonferroni']] <- p.adjust(df_proteins_ine_casecontrol[['Pr...t..']], method = 'bonferroni')
df_proteins_ine_casecontrol <- df_proteins_ine_casecontrol[order(df_proteins_ine_casecontrol[['Pr...t..']]), ]

# fit model for only the first three plates
df_proteins_plate123 <- do_regression(protein_data[protein_data[['PlateID']] %in% c('PlateLayout_Plate1', 'PlateLayout_Plate2', 'PlateLayout_Plate3'), ], random_effects=c('sample_final', 'PlateID'), fixed_effects=c('age', 'sex', 'case_control'))
df_proteins_plate123_casecontrol <- df_proteins_plate123[df_proteins_plate123[['term']] == 'case_controlcontrol', ]
df_proteins_plate123_casecontrol[['BH']] <- p.adjust(df_proteins_plate123_casecontrol[['Pr...t..']], method = 'BH')
df_proteins_plate123_casecontrol[['p.bonferroni']] <- p.adjust(df_proteins_plate123_casecontrol[['Pr...t..']], method = 'bonferroni')
df_proteins_plate123_casecontrol <- df_proteins_plate123_casecontrol[order(df_proteins_plate123_casecontrol[['Pr...t..']]), ]
df_proteins_plate123_casecontrol[['z']] <- df_proteins_plate123_casecontrol[['Estimate']] / df_proteins_plate123_casecontrol[['Std..Error']]

# let's also try to do plate one and two separately
df_proteins_plate1 <- do_regression(protein_data[protein_data[['PlateID']] == 'PlateLayout_Plate1', ], random_effects=c('sample_final'), fixed_effects=c('age', 'sex', 'case_control'))
df_proteins_plate1_casecontrol <- df_proteins_plate1[df_proteins_plate1[['term']] == 'case_controlcontrol', ]
df_proteins_plate1_casecontrol[['BH']] <- p.adjust(df_proteins_plate1_casecontrol[['Pr...t..']], method = 'BH')
df_proteins_plate1_casecontrol[['p.bonferroni']] <- p.adjust(df_proteins_plate1_casecontrol[['Pr...t..']], method = 'bonferroni')
df_proteins_plate1_casecontrol <- df_proteins_plate1_casecontrol[order(df_proteins_plate1_casecontrol[['Pr...t..']]), ]
df_proteins_plate2 <- do_regression(protein_data[protein_data[['PlateID']] == 'PlateLayout_Plate2', ], random_effects=c(), fixed_effects=c('age', 'sex', 'case_control'))
df_proteins_plate2_casecontrol <- df_proteins_plate2[df_proteins_plate2[['term']] == 'case_controlcontrol', ]
df_proteins_plate2_casecontrol[['BH']] <- p.adjust(df_proteins_plate2_casecontrol[['Pr...t..']], method = 'BH')
df_proteins_plate2_casecontrol[['p.bonferroni']] <- p.adjust(df_proteins_plate2_casecontrol[['Pr...t..']], method = 'bonferroni')
df_proteins_plate2_casecontrol <- df_proteins_plate2_casecontrol[order(df_proteins_plate2_casecontrol[['Pr...t..']]), ]
df_proteins_plate3 <- do_regression(protein_data[protein_data[['PlateID']] == 'PlateLayout_Plate3', ], random_effects=c(), fixed_effects=c('age', 'sex', 'case_control'))
df_proteins_plate3_casecontrol <- df_proteins_plate3[df_proteins_plate3[['term']] == 'case_controlcontrol', ]
df_proteins_plate3_casecontrol[['BH']] <- p.adjust(df_proteins_plate3_casecontrol[['Pr...t..']], method = 'BH')
df_proteins_plate3_casecontrol[['p.bonferroni']] <- p.adjust(df_proteins_plate3_casecontrol[['Pr...t..']], method = 'bonferroni')
df_proteins_plate3_casecontrol <- df_proteins_plate3_casecontrol[order(df_proteins_plate3_casecontrol[['Pr...t..']]), ]
# merge the results of these runs
df_proteins_plates_casecontrol <- merge(df_proteins_plate1_casecontrol[, c('protein', 'Estimate', 'Std..Error', 't.value', 'Pr...t..')], df_proteins_plate2_casecontrol[, c('protein', 'Estimate', 'Std..Error', 't.value', 'Pr...t..')], by = 'protein')
# set colnames
colnames(df_proteins_plates_casecontrol) <- c('protein', 'Estimate.1', 'Std..Error.1', 't.value.1', 'P.1', 'Estimate.2',  'Std.Error.2', 't.value.2', 'P.2')
# merge third table
df_proteins_plates_casecontrol <- merge(df_proteins_plates_casecontrol, df_proteins_plate3_casecontrol[, c('protein', 'Estimate', 'Std..Error', 't.value', 'Pr...t..')], by = 'protein')
# set colnames
colnames(df_proteins_plates_casecontrol) <- c('protein', 'Estimate.1', 'Std.Error.1', 't.value.1', 'P.1', 'Estimate.2',  'Std.Error.2', 't.value.2', 'P.2', 'Estimate.3', 'Std.Error.3', 't.value.3', 'P.3')
# join the p of the full analysis on there
df_proteins_plates_casecontrol[['P.all']] <- df_proteins_casecontrol[match(df_proteins_plates_casecontrol[['protein']], df_proteins_casecontrol[['protein']]), 'Pr...t..']
# and the one of the first three plates
df_proteins_plates_casecontrol[['P.123']] <- df_proteins_plate3_casecontrol[match(df_proteins_plates_casecontrol[['protein']], df_proteins_plate3_casecontrol[['protein']]), 'Pr...t..']
# add z scores
df_proteins_plates_casecontrol[['z.1']] <- df_proteins_plates_casecontrol[['Estimate.1']] / df_proteins_plates_casecontrol[['Std.Error.1']]
df_proteins_plates_casecontrol[['z.2']] <- df_proteins_plates_casecontrol[['Estimate.2']] / df_proteins_plates_casecontrol[['Std.Error.2']]
df_proteins_plates_casecontrol[['z.3']] <- df_proteins_plates_casecontrol[['Estimate.3']] / df_proteins_plates_casecontrol[['Std.Error.3']]
# check the z score of using plate 4 and 5 or not
ggplot(data = merge(df_proteins_casecontrol, df_proteins_plate123_casecontrol, by = 'protein'), mapping = aes(x = z.x, y = z.y)) + 
  geom_point() +
  xlab('Z score all plates') +
  ylab('Z score plates 1,2,3')
# and the -log10 p values
ggplot(data = merge(df_proteins_casecontrol, df_proteins_plate123_casecontrol, by = 'protein'), mapping = aes(x = -log10(Pr...t...x), y = -log10(Pr...t...y))) + 
  geom_point() +
  xlab('-log10 p all plates') +
  ylab('-log10 p plates 1,2,3')

# let's plot per control and case
plot_olink_expression(olink = protein_data, protein_name = 'OID20524', plot_order = c('control-old', 'control-recent', 'case'), group_column = 'casecontrolpandemic') + scale_fill_manual(values = roycols::get_color_list(protein_data$casecontrolpandemic)) + xlab('sample status')
# check how much the plates differ
plot_olink_expression(olink = protein_data, protein_name = 'OID20524', group_column = 'PlateID') + scale_fill_manual(values = roycols::get_color_list(protein_data$PlateID)) + xlab('plate')
# check how case/control differs in each plate separately
plot_grid(
  plot_olink_expression(olink = protein_data[protein_data[['PlateID']] == 'PlateLayout_Plate1' & !((grepl('^(PC)|(NC)|(SC)', protein_data$sample_final))), ], protein_name = 'OID20524', plot_order = c('control-old', 'control-recent', 'case'), group_column = 'casecontrolpandemic') + scale_fill_manual(values = roycols::get_color_list(protein_data$casecontrolpandemic)) + xlab('sample status') + ylim(c(-2,10)) + ggtitle('OID20524 plate 1'),
  plot_olink_expression(olink = protein_data[protein_data[['PlateID']] == 'PlateLayout_Plate2' & !((grepl('^(PC)|(NC)|(SC)', protein_data$sample_final))), ], protein_name = 'OID20524', plot_order = c('control-old', 'control-recent', 'case'), group_column = 'casecontrolpandemic') + scale_fill_manual(values = roycols::get_color_list(protein_data$casecontrolpandemic)) + xlab('sample status') + ylim(c(-2,10)) + ggtitle('OID20524 plate 2'),
  plot_olink_expression(olink = protein_data[protein_data[['PlateID']] == 'PlateLayout_Plate3' & !((grepl('^(PC)|(NC)|(SC)', protein_data$sample_final))), ], protein_name = 'OID20524', plot_order = c('control-old', 'control-recent', 'case'), group_column = 'casecontrolpandemic') + scale_fill_manual(values = roycols::get_color_list(protein_data$casecontrolpandemic)) + xlab('sample status') + ylim(c(-2,10)) + ggtitle('OID20524 plate 3'),
  plot_olink_expression(olink = protein_data[protein_data[['PlateID']] == 'PlateLayout_Plate4' & !((grepl('^(PC)|(NC)|(SC)', protein_data$sample_final))), ], protein_name = 'OID20524', plot_order = c('control-old', 'control-recent', 'case'), group_column = 'casecontrolpandemic') + scale_fill_manual(values = roycols::get_color_list(protein_data$casecontrolpandemic)) + xlab('sample status') + ylim(c(-2,10)) + ggtitle('OID20524 plate 4'),
  plot_olink_expression(olink = protein_data[protein_data[['PlateID']] == 'PlateLayout_Plate5' & !((grepl('^(PC)|(NC)|(SC)', protein_data$sample_final))), ], protein_name = 'OID20524', plot_order = c('control-old', 'control-recent', 'case'), group_column = 'casecontrolpandemic') + scale_fill_manual(values = roycols::get_color_list(protein_data$casecontrolpandemic)) + xlab('sample status') + ylim(c(-2,10)) + ggtitle('OID20524 plate 5'),
  nrow = 3,
  ncol = 2
)

# plot the significant ones
plot_grid(
  plot_olink_expression(olink = protein_data[protein_data$pandemic == 'postpandemic', ], protein_name = 'OID20477', plot_order = c('control', 'case')) + scale_fill_manual(values = roycols::get_color_list(protein_data$case_control)) + xlab('sample status') + ggtitle('Interleukin-17C'),
  plot_olink_expression(olink = protein_data[protein_data$pandemic == 'postpandemic', ], protein_name = 'OID20504', plot_order = c('control', 'case')) + scale_fill_manual(values = roycols::get_color_list(protein_data$case_control)) + xlab('sample status') + ggtitle('Integral membrane protein 2A'),
  plot_olink_expression(olink = protein_data[protein_data$pandemic == 'postpandemic', ], protein_name = 'OID20524', plot_order = c('control', 'case')) + scale_fill_manual(values = roycols::get_color_list(protein_data$case_control)) + xlab('sample status') + ggtitle('Dual adapter for phosphotyrosine and\n3-phosphotyrosine and 3-phosphoinositide'),
  plot_olink_expression(olink = protein_data[protein_data$pandemic == 'postpandemic', ], protein_name = 'OID20563', plot_order = c('control', 'case')) + scale_fill_manual(values = roycols::get_color_list(protein_data$case_control)) + xlab('sample status') + ggtitle('Interleukin-6'),
  plot_olink_expression(olink = protein_data[protein_data$pandemic == 'postpandemic', ], protein_name = 'OID20577', plot_order = c('control', 'case')) + scale_fill_manual(values = roycols::get_color_list(protein_data$case_control)) + xlab('sample status') + ggtitle('Interleukin-1 receptor-associated kinase 4'),
  plot_olink_expression(olink = protein_data[protein_data$pandemic == 'postpandemic', ], protein_name = 'OID20717', plot_order = c('control', 'case')) + scale_fill_manual(values = roycols::get_color_list(protein_data$case_control)) + xlab('sample status') + ggtitle('Tyrosine-protein phosphatase non-receptor type 6'),
  nrow = 3,
  ncol = 2
)
plot_grid(
  plot_olink_expression(olink = protein_data, protein_name = 'OID20477', plot_order = c('control', 'case')) + scale_fill_manual(values = roycols::get_color_list(protein_data$case_control)) + xlab('sample status'),
  plot_olink_expression(olink = protein_data, protein_name = 'OID20504', plot_order = c('control', 'case')) + scale_fill_manual(values = roycols::get_color_list(protein_data$case_control)) + xlab('sample status'),
  plot_olink_expression(olink = protein_data, protein_name = 'OID20524', plot_order = c('control', 'case')) + scale_fill_manual(values = roycols::get_color_list(protein_data$case_control)) + xlab('sample status'),
  plot_olink_expression(olink = protein_data, protein_name = 'OID20563', plot_order = c('control', 'case')) + scale_fill_manual(values = roycols::get_color_list(protein_data$case_control)) + xlab('sample status'),
  plot_olink_expression(olink = protein_data, protein_name = 'OID20577', plot_order = c('control', 'case')) + scale_fill_manual(values = roycols::get_color_list(protein_data$case_control)) + xlab('sample status'),
  plot_olink_expression(olink = protein_data, protein_name = 'OID20717', plot_order = c('control', 'case')) + scale_fill_manual(values = roycols::get_color_list(protein_data$case_control)) + xlab('sample status'),
  nrow = 3,
  ncol = 2
)
