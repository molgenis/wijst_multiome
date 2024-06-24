#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_gse_enrichr.R
# Function: perform gene set enrichment with REACTOME and enrichr
############################################################################################################################

####################
# libraries        #
####################

# for doing the enrichment
library(enrichR)
# for parallel processing
# library(foreach)


####################
# Functions        #
####################


#' get each gene and their LFC from a MAST output file created using the Seurat version of MAST
#' 
#' @param mast.output.loc the location of the output file
#' @param p.val.column the column to use for filtering the significant genes
#' @param p.val.cutoff the cutoff value to use for filtering the significant genes
#' @returns a table with the DE genes filtered by the given significance level
#' 
get_gene_and_lfc_seurat_mast <- function(mast.output.loc, p.val.column='p_val_adj', p.val.cutoff=0.05){
  # read the output
  de.output <- read.table(mast.output.loc, header = T, row.names = 1, sep = '\t', stringsAsFactors = F)
  # subset to what is significant
  de.output.significant <- de.output[de.output[[p.val.column]] < p.val.cutoff, ]
  return(de.output.significant)
}


#' get each gene and their LFC from a limma output file
#' 
#' @param limma.output.loc the location of the output file
#' @param p.val.column the column to use for filtering the significant genes
#' @param p.val.cutoff the cutoff value to use for filtering the significant genes
#' @returns a table with the DE genes filtered by the given significance level
#' 
get_gene_and_lfc_limma <- function(limma.output.loc, p.val.column='adj.P.Val', p.val.cutoff=0.05){
  # read the output
  de.output <- read.table(limma.output.loc, header = T, row.names = 1, sep = '\t', stringsAsFactors = F)
  # subset to what is significant
  de.output.significant <- de.output[de.output[[p.val.column]] < p.val.cutoff, ]
  return(de.output.significant)
}


#' get each gene and their LFC from a MAST output file created using the standalone version of MAST
#' 
#' @param mast.output.loc the location of the output file
#' @param p.val.column the column to use for filtering the significant genes
#' @param p.val.cutoff the cutoff value to use for filtering the significant genes
#' @param contrast the value in the model that we are interested in
#' @returns a table with the DE genes filtered by the given significance level
#' 
get_gene_and_lfc_standalone_mast <- function(mast.output.loc, p.val.column='bonferroni', p.val.cutoff=0.05, contrast='inflammation_statusNI'){
  # read the output
  de.output <- read.table(mast.output.loc, header = T, sep = '\t', stringsAsFactors = F)
  # subset to just the contrast we care about
  de.output.contrast <- de.output[de.output[['contrast']] == contrast &
                                    de.output[['component']] == 'H', ]
  # subset to just what is significant
  de.output.contrast.significant <- de.output.contrast[!is.na(de.output.contrast[[p.val.column]]) & de.output.contrast[[p.val.column]] < p.val.cutoff, ]
  # we can now safely set the rownames to be the genes, as with one contrast, there is only one entry per gene
  rownames(de.output.contrast.significant) <- de.output.contrast.significant[['primerid']]
  return(de.output.contrast.significant)
}


#' get each gene and their LFC from a DE output file
#' 
#' @param de.output.loc the location of the output file
#' @param de_method the methods used to generate the output file (mast_standalone, mast_seurat, limma)
#' @param p.val.column the column to use for filtering the significant genes
#' @param p.val.cutoff the cutoff value to use for filtering the significant genes
#' @param contrast the value in the model that we are interested in (only necessary for mast_standalone)
#' @returns a table with the DE genes filtered by the given significance level
#' 
get_de_output <- function(de.output.loc, de_method, p.val.column='bonferroni', p.val.cutoff=0.05, contrast='inflammation_statusNI'){
  # use the correct method
  de.output <- NULL
  if (de_method == 'mast_standalone') {
    de.output <- get_gene_and_lfc_standalone_mast(de.output.loc, p.val.column = p.val.column, p.val.cutoff = p.val.cutoff)
  }
  else if (de_method == 'mast_seurat') {
    de.output <- get_gene_and_lfc_seurat_mast(de.output.loc, p.val.column = p.val.column, p.val.cutoff = p.val.cutoff, contrast = contrast)
  }
  else if (de_method == 'limma'){
    de.output <- get_gene_and_lfc_limma(de.output.loc, p.val.column = p.val.column, p.val.cutoff = p.val.cutoff)
  }
  else{
    print('de_method not valid, valid options are: mast_standalone, mast_seurat, limma')
  }
  return(de.output)
}


#' use enrichr to do gene set enrichment on each DE output file in a directory
#' 
#' @param mast.output.loc the folder containing the DE output files
#' @param pathway.output.loc the folder to place the GSE results
#' @param de_method the methods used to generate the output file (mast_standalone, mast_seurat, limma)
#' @param p.val.column the column to use for filtering the significant genes
#' @param p.val.cutoff the cutoff value to use for filtering the significant genes
#' @param log.fc.column the column containing the log fold change of the gene
#' @param contrast.var the value in the model that we are interested in (only necessary for mast_standalone)
#' @param list.regex regular expression pattern to use when lising the DE output files
#' @param do_positive test just the positive LFC genes
#' @param do_negative test just the negative LFC genes
#' @param do_joint test genes regardless of LFC direction
#' @returns 0 if succesfull
#' do_enrichr(limma_lowm20_output_path, gse_enrichr_limma_lowm20_output_path, de_method = 'limma', log.fc.column = 'logFC', p.val.column = 'p.bonferroni')
do_enrichr <- function(mast.output.loc, pathway.output.loc, de_method='mast_seurat', log.fc.column='avg_log2FC', p.val.column='p_val_adj', p.val.cutoff=0.05, contrast.var = NULL, list.regex = NULL, do_positive=T, do_negative=T, do_joint=T){
  # check if the output folder exists
  if(!dir.exists(pathway.output.loc)){
    dir.create(pathway.output.loc, recursive = T)
  }
  # get the DE output files
  mast.output.files <- list.files(mast.output.loc, full.names = T, pattern = list.regex)
  # check each cell type
  for(mast.output.file in mast.output.files){
    #foreach(mast.output.file=mast.output.files) %do% {
    # get the basename
    mast.output.file.basename <- basename(mast.output.file)
    # get the comparison from that
    comparison <- sub(pattern = "(.*?)\\..*$", replacement = "\\1", x = mast.output.file.basename)
    # get the full path to the DE file
    de.output.loc <- mast.output.file
    # read the output
    #de.output <- read.table(de.output.loc, header = T, row.names = 1, sep = '\t', stringsAsFactors = F)
    # get only the significant result
    #de.output <- de.output[de.output[[p.val.column]] < p.val.cutoff, ]
    de.output <- get_de_output(de.output.loc = de.output.loc, de_method = de_method, p.val.column = p.val.column, p.val.cutoff = p.val.cutoff, contrast = contrast.var)
    # check if there are any results
    if(nrow(de.output) > 0){
      if (do_joint==T) {
        # grab the genes
        gene.symbols.sig <- rownames(de.output)
        # do the enrichment analysis
        enriched <- enrichr(genes = gene.symbols.sig, databases = c('Reactome_2016'))
        # grab the reactome result
        enriched.reactome <- enriched[['Reactome_2016']]
        # write a result if there is one
        if(nrow(enriched.reactome) > 0){
          # paste the output path together
          enriched.reactome.loc <- gzfile(paste(pathway.output.loc, comparison, '.tsv.gz', sep = ''))
          # write the result
          write.table(enriched.reactome, enriched.reactome.loc, sep = '\t', quote = T, row.names = F, col.names = T)
        }
      }
      # now check for the negative genes
      de.output.negative <- de.output[de.output[[log.fc.column]] < 0, ]
      # check if there are any genes left
      if(nrow(de.output.negative) > 0 & do_negative==T){
        # grab the genes
        gene.symbols.sig.negative <- rownames(de.output.negative)
        # do the enrichment analysis
        enriched.negative <- enrichr(genes = gene.symbols.sig.negative, databases = c('Reactome_2016'))
        # grab the reactome result
        enriched.reactome.negative <- enriched.negative[['Reactome_2016']]
        # write a result if there is one
        if(nrow(enriched.reactome.negative) > 0){
          # paste the output path together
          enriched.reactome.negative.loc <- gzfile(paste(pathway.output.loc, comparison, '_negative', '.tsv.gz', sep = ''))
          # write the result
          write.table(enriched.reactome.negative, enriched.reactome.negative.loc, sep = '\t', quote = T, row.names = F, col.names = T)
        }
      }
      # now check for the positive genes
      de.output.positive <- de.output[de.output[[log.fc.column]] > 0, ]
      # check if there are any genes left
      if(nrow(de.output.positive) > 0 & do_positive==T){
        # grab the genes
        gene.symbols.sig.positive <- rownames(de.output.positive)
        # do the enrichment analysis
        enriched.positive <- enrichr(genes = gene.symbols.sig.positive, databases = c('Reactome_2016'))
        # grab the reactome result
        enriched.reactome.positive <- enriched.positive[['Reactome_2016']]
        # write a result if there is one
        if(nrow(enriched.reactome.positive) > 0){
          # paste the output path together
          enriched.reactome.positive.loc <- gzfile(paste(pathway.output.loc, comparison, '_positive', '.tsv.gz', sep = ''))
          # write the result
          write.table(enriched.reactome.positive, enriched.reactome.positive.loc, sep = '\t', quote = T, row.names = F, col.names = T)
        }
      }
    }
    else{
      print(paste('no genes significant for', mast.output.file))
    }
  }
  return(0)
}

####################
# settings         #
####################

setEnrichrSite("Enrichr") # Human genes

####################
# Main Code        #
####################

# locations of limma output
limma_output_mo_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/differential_expression/limma_dream/output/stimulation/'
limma_output_1m_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/differential_expression/limma_dream/output/stimulation_1m/'

# location of GSE output
gse_output_mo_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/gse/enrichr/limma_dream/output/stimulation/'
gse_output_1m_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/gse/enrichr/limma_dream/output/stimulation_1m/'

# create paths if they don't exist
dir.create(gse_output_mo_loc, recursive = T, showWarnings = F)
dir.create(gse_output_1m_loc, recursive = T, showWarnings = F)

# and do enrichr
do_enrichr(limma_output_mo_loc, gse_output_mo_loc, de_method = 'limma', log.fc.column = 'logFC', p.val.column = 'p.bonferroni', do_positive = T, do_negative = F, do_joint = F)
do_enrichr(limma_output_mo_loc, gse_output_mo_loc, de_method = 'limma', log.fc.column = 'logFC', p.val.column = 'p.bonferroni', do_positive = F, do_negative = T, do_joint = F)
do_enrichr(limma_output_mo_loc, gse_output_mo_loc, de_method = 'limma', log.fc.column = 'logFC', p.val.column = 'p.bonferroni', do_positive = F, do_negative = F, do_joint = T)
do_enrichr(limma_output_1m_loc, gse_output_1m_loc, de_method = 'limma', log.fc.column = 'logFC', p.val.column = 'p.bonferroni', do_positive = T, do_negative = F, do_joint = F)
do_enrichr(limma_output_1m_loc, gse_output_1m_loc, de_method = 'limma', log.fc.column = 'logFC', p.val.column = 'p.bonferroni', do_positive = F, do_negative = T, do_joint = F)
do_enrichr(limma_output_1m_loc, gse_output_1m_loc, de_method = 'limma', log.fc.column = 'logFC', p.val.column = 'p.bonferroni', do_positive = F, do_negative = F, do_joint = T)
