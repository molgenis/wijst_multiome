#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_get_ncells_analysis.R
# Function: get the number of cells used in the QTL mapping analysis
############################################################################################################################

####################
# libraries        #
####################

# parse command line arguments
library(optparse)
# to make md5 checksums
library(mdfiver)

####################
# Functions        #
####################

#' Get Cell Numbers Analysis
#'
#' This function reads QTL files from a specified directory, filters them based on a given pattern, and calculates the total number of cells for each cell type.
#'
#' @param qtl_input_dir Character. The directory containing the QTL files.
#' @param qtl_prepend Character. A string to prepend to the QTL file names (default is an empty string).
#' @param qtl_append Character. A string to append to the QTL file names (default is '.covariates.txt.gz').
#' @param cell_count_column Character. The name of the column containing cell counts (default is 'CellCount').
#'
#' @return A data frame with two columns: 'cell_type' and 'cell_number', representing the cell type and the total number of cells, respectively.
#'
#' @examples
#' \dontrun{
#' qtl_input_dir <- "path/to/qtl/files"
#' result <- get_ncells_analysis(qtl_input_dir)
#' print(result)
#' }
#'
get_ncells_analysis <- function(qtl_input_dir, qtl_prepend='', qtl_append='.covariates.txt.gz', cell_count_column='CellCount') {
    # list the files in the directory
    qtl_files <- list.files(qtl_input_dir)
    # create the regex to filter the files
    qtl_file_regex <- paste('.*', qtl_append, sep = '')
    # filter the files
    qtl_files_filtered <- qtl_files[grep(qtl_file_regex, qtl_files)]
    # create a list to store the cell numbers
    cell_nr_per_celltype <- list()
    # check each file
    for (qtl_file in qtl_files_filtered) {
        # read the file
        qtl_file_table <- read.table(paste(qtl_input_dir, '/', qtl_file, sep = ''), header = T, sep = '\t', comment.char = '')
        # get the number of cells
        cell_numbers_celltype <- qtl_file_table[[cell_count_column]]
        # remove NA, those would be zero anyway
        cell_numbers_celltype <- cell_numbers_celltype[!is.na(cell_numbers_celltype)]
        # get the sum
        cell_numbers_celltype_all <- sum(cell_numbers_celltype)
        # extract the cell type from the file name
        celltype_name <- gsub(qtl_prepend, '', qtl_file)
        celltype_name <- gsub(qtl_append, '', celltype_name)
        # put in the list under this name
        cell_nr_per_celltype[[celltype_name]] <- cell_numbers_celltype_all
    }
    # turn into a dataframe
    cell_numbers_df <- data.frame('cell_type' = names(cell_nr_per_celltype), 'cell_number' = as.vector(unlist(cell_nr_per_celltype)))
    return(cell_numbers_df)
}

####################
# Main Code        #
####################

# make command line options
option_list <- list(
  make_option(c("-i", "--qtl_dir"), type="character", default=NULL,
              help="directory containing QTL input files", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL,
              help="tsv output file location", metavar="character")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# check if parameters are filled
if (is.null(opt[['qtl_dir']])) {
    stop('--qtl_dir/-i not supplied')
}
if (is.null(opt[['qtl_dir']])) {
    stop('--output_file/-o not supplied')
}

# get the numbers
cell_numbers <- get_ncells_analysis(opt[['qtl_dir']])
# get the output location
output_loc_gz <- opt[['output_file']]
# gz file ends with .gz
if (grepl('.gz$', output_loc_gz)) {
    # gzip if ends with .gz
    output_loc_gz <- gzfile(output_loc_gz)
}
write.table(cell_numbers, output_loc_gz, sep = '\t', row.names = F, col.names = T, quote = F)
# also make a checksum
mdfiver::create_md5_for_file(opt[['output_file']])