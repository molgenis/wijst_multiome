#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_differential_accessibility_merge_chunks.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(mdfiver)

####################
# Functions        #
####################

merge_chunked_topic <- function(dar_output_base_loc) {
  # check each folder, which is a cell type
  celltypes <- list.dirs(dar_output_base_loc, full.names = F, recursive = F)
  # check each cell type
  for (celltype in celltypes) {
    # now list all of the chunks
    celltype_dir <- paste(dar_output_base_loc, '/', celltype, '/', sep = '')
    chunk_dirs <- list.dirs(celltype_dir, full.names = F, recursive = F)
    # we'll save each topic
    topics_list <- list()
    # check each chunk
    for (chunk in chunk_dirs) {
      # paste together the pattern for each of the output files
      topic_files_dir <- paste(celltype_dir, '/', chunk, '/', sep = '')
      topic_files_pattern <- paste('*.tsv.gz', sep = '')
      # now check all of the files in the chunk directory
      all_files <- list.files(topic_files_dir, pattern = topic_files_pattern, full.names = F, recursive = F)
      # now filter on the ones we wan
      topic_files <- all_files[grepl(paste(celltype, '_Topic\\d+.tsv.gz$', sep = ''), all_files)]
      # check each topic file
      for (topic_file in topic_files) {
        # extract the topic from there
        topic <- gsub(paste(celltype, '_', sep = ''), '', topic_file)
        topic <- gsub('.tsv.gz', '', topic)
        # read the topic file
        topic_chunk_loc <- paste(topic_files_dir, topic_file, sep = '')
        topic_chunk <- read.table(topic_chunk_loc, header = T, sep = '\t')
        # check if this topic exists in the list
        if (!(topic %in% names(topics_list))) {
          topics_list[[topic]] <- list()
        }
        # add the topic and the chunk to the output file
        topic_chunk[['topic']] <- topic
        topic_chunk[['chunk']] <- chunk
        # add the chunk to that topic
        topics_list[[topic]][[chunk]] <- topic_chunk
      }
    }
    # now let us go over each topic
    for (topic_name in names(topics_list)) {
      # combine the chunks
      topic_combined <- do.call('rbind', topics_list[[topic_name]])
      # redo the B&H, as it was done per chunk, while it needs to be done on all data
      topic_combined[['adj.P.Val']] <- p.adjust(topic_combined[['P.Value']], method = 'BH')
      # add bonferroni correction
      topic_combined[['p.bonferroni']] <- p.adjust(topic_combined[['P.Value']], method = 'bonferroni')
      # order by the original p value
      topic_combined <- topic_combined[order(topic_combined[['P.Value']]), ]
      
      # paste together the output path
      topic_output_path <- paste(celltype_dir, '/', celltype, '_', topic_name, '.tsv.gz', sep = '')
      # write the result
      write.table(topic_combined, gzfile(topic_output_path), sep = '\t', row.names = F)
      # write md5
      mdfiver::create_md5_for_file(topic_output_path)
      
      # do the nominally significant ones as well
      topic_output_nominal_path <- paste(celltype_dir, '/', celltype, '_', topic_name, '.nominal_significant.tsv.gz', sep = '')
      topic_combined_nominal <- topic_combined[topic_combined[['P.Value']] < 0.05, ]
      write.table(topic_combined_nominal, gzfile(topic_output_nominal_path), sep = '\t', row.names = F)
      mdfiver::create_md5_for_file(topic_output_nominal_path)
    }
  }
}


####################
# Main Code        #
####################

# dar output folder
dar_output_base_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/output/topics20_otsu_imputed/'
merge_chunked_topic(dar_output_base_loc)
