#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Martijn van der Werf, Roy Oelen
# Name: imputeCelltypes.R
# Function: based on RNA celltype annotation and ATAC clusters with LSI, impute celltypes for cells with no labels in ATAC
#
############################################################################################################################


setwd("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/archr_preprocess_samples/")
library(ArchR)

addArchRThreads(threads=4)
addArchRGenome("hg38")
h5disableFileLocking()

proj <- loadArchRProject('Project_common_atac_rna_cells')

################################
# IMPUTE ONLY MISSING CELL TYPES
################################
archr_metadata <- proj@cellColData
archr_metadata$celltype[archr_metadata$celltype=="unannotated"] <- NA
archr_metadata$cell_type_missing_imputed <- NA
# go through the grouping we have for the entire object
for(group in unique(archr_metadata$Clusters)){
    # subset to get only this group
    archr_group <- archr_metadata[archr_metadata$Clusters == group,]
    best_group <- 'unknown'
    best_number <- 0
    # check against the reference column
    for(reference in unique(archr_metadata$celltype)){
      # we don't care for the NA reference, if we had all data, we wouldn't need to do this anyway
      if(is.na(reference) == F){
        # grab the number of cells in this group, with this reference
        number_of_reference_in_group <- nrow(archr_group[archr_group$celltype == reference & is.na(archr_group$celltype) == F,])
        correctpercent <- number_of_reference_in_group/ncol(archr_group)
        print(paste(group,"matches",reference,correctpercent,sep=" "))
        # update numbers if better match
        if(number_of_reference_in_group > best_number){
          best_number <- number_of_reference_in_group
          best_group <- reference
        }
      }
    }
    print(paste("setting ident:",best_group,"for group", group, sep=" "))
    # set this best identity
    archr_metadata[archr_metadata$Clusters == group,]$cell_type_missing_imputed <- ifelse(!is.na(archr_group$celltype), archr_group$celltype, best_group)
    # force cleanup
    rm(archr_group)
}

#####################################
# IMPUTE ALL CELL TYPES
#####################################
archr_metadata <- proj@cellColData
archr_metadata$celltype[archr_metadata$celltype=="unannotated"] <- NA
archr_metadata$cell_type_lowerres_imputed <- NA
# go through the grouping we have for the entire object
for(group in unique(archr_metadata$Clusters)){
    # subset to get only this group
    archr_group <- archr_metadata[archr_metadata$Clusters == group,]
    best_group <- 'unknown'
    best_number <- 0
    # check against the reference column
    for(reference in unique(archr_metadata$celltype)){
      # we don't care for the NA reference, if we had all data, we wouldn't need to do this anyway
      if(is.na(reference) == F){
        # grab the number of cells in this group, with this reference
        number_of_reference_in_group <- nrow(archr_group[archr_group$celltype == reference & is.na(archr_group$celltype) == F,])
        correctpercent <- number_of_reference_in_group/ncol(archr_group)
        print(paste(group,"matches",reference,correctpercent,sep=" "))
        # update numbers if better match
        if(number_of_reference_in_group > best_number){
          best_number <- number_of_reference_in_group
          best_group <- reference
        }
      }
    }
    print(paste("setting ident:",best_group,"for group", group, sep=" "))
    # set this best identity
    archr_metadata[archr_metadata$Clusters == group,]$cell_type_lowerres_imputed <- best_group
    # force cleanup
    rm(archr_group)
}

#############################
# PLOT LSI IMPUTED CELL TYPES 
#############################

# Plot LSI clusters for only missing cell type imputed cells
umap_celltypes_missing_imputed <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cell_type_missing_imputed", embedding = "UMAP")
plotPDF(umap_celltypes_missing_imputed, name = "UMAP_celltypes_missing_imputed.pdf", ArchRProj = proj, addDOC = FALSE, width = 7, height = 7)

# Plot LSI clusters for all cell type imputed cells
umap_clusters_celltypes_imputed <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cell_type_lowerres_imputed", embedding = "UMAP")
plotPDF(umap_clusters_celltypes_imputed, name = "UMAP_celltypes_all_imputed.pdf", ArchRProj = proj, addDOC = FALSE, width = 7, height = 7)

# Plot clusters based on LSI
umap_clusters <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(umap_clusters, name = "UMAP_clusters.pdf", ArchRProj = proj, addDOC = FALSE, width = 7, height = 7)
#####################
# Check LSI components

metadata <- read.csv('../metadata/mo_cellevel_metadata.tsv', sep='\t')
metadata$archr_barcode <- paste0(metadata$batch, "#", metadata$barcode_1)
proj_common@cellColData['barcode_archr'] <- rownames(proj_common@cellColData)

reduced_dims <- getReducedDims(proj_common)
reduced_dims <- as.data.frame(reduced_dims)
reduced_dims$stimulation <- NA
reduced_dims$stimulation[which(rownames(reduced_dims) %in% metadata$archr_barcode)] <- metadata$filtered_condition_imputed
reduced_dims$seq_depth <- proj_common@cellColData$nFrags[match(rownames(proj_common@cellColData), rownames(reduced_dims))]
reduced_dims$sample <- str_sub(rownames(reduced_dims), start=1, end=12)

p1 <- ggplot(reduced_dims[complete.cases(reduced_dims),], aes(x=LSI1, y=LSI2, color=seq_depth)) + geom_point()

# Split dataframe by sample
reduced_dims_per_sample <- split(reduced_dims, reduced_dims$sample)

p1 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(1:8)]), aes(x=LSI1, y=LSI2, color=sample)) + geom_point(size=0.1) +scale_color_manual(values=palette()) + guides(colour = guide_legend(override.aes = list(size=6)))
p2 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(9:16)]), aes(x=LSI1, y=LSI2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p3 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(17:24)]), aes(x=LSI1, y=LSI2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette()) + guides(colour = guide_legend(override.aes = list(size=6)))
p4 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(25:32)]), aes(x=LSI1, y=LSI2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p5 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(33:40)]), aes(x=LSI1, y=LSI2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p6 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(41:48)]), aes(x=LSI1, y=LSI2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p7 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(49:56)]), aes(x=LSI1, y=LSI2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p8 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(57:64)]), aes(x=LSI1, y=LSI2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p9 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(65:72)]), aes(x=LSI1, y=LSI2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p10 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(73:80)]), aes(x=LSI1, y=LSI2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))

plotList <- list(p1, p2, p3, p4, p5, p6, p7, p8,p9,p10)

ggsave(file = 'Project_common_cells/Plots/subsplot_samples_LSI1_LSI2.pdf', width=15, height=15, arrangeGrob(grobs = plotList, ncol = 2))

# Check UMAP as well
reduced_dims <- cbind(reduced_dims, proj_common@embeddings$UMAP$df)
colnames(reduced_dims)[34] <- "UMAP1"
colnames(reduced_dims)[35] <- "UMAP2"
reduced_dims_per_sample <- split(reduced_dims, reduced_dims$sample)

p1 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(1:8)]), aes(x=UMAP1, y=UMAP2, color=sample)) + geom_point(size=0.1) +scale_color_manual(values=palette()) + guides(colour = guide_legend(override.aes = list(size=6)))
p2 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(9:16)]), aes(x=UMAP1, y=UMAP2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p3 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(17:24)]), aes(x=UMAP1, y=UMAP2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette()) + guides(colour = guide_legend(override.aes = list(size=6)))
p4 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(25:32)]), aes(x=UMAP1, y=UMAP2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p5 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(33:40)]), aes(x=UMAP1, y=UMAP2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p6 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(41:48)]), aes(x=UMAP1, y=UMAP2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p7 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(49:56)]), aes(x=UMAP1, y=UMAP2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p8 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(57:64)]), aes(x=UMAP1, y=UMAP2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p9 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(65:72)]), aes(x=UMAP1, y=UMAP2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))
p10 <- ggplot(do.call('rbind', reduced_dims_per_sample[c(73:80)]), aes(x=UMAP1, y=UMAP2, color=sample)) + geom_point(size=0.1)+scale_color_manual(values=palette())+ guides(colour = guide_legend(override.aes = list(size=6)))

plotList <- list(p1, p2, p3, p4, p5, p6, p7, p8,p9,p10)

ggsave(file = 'Project_common_atac_rna_cells/Plots/subsplot_samples_UMAP1_UMAP2.pdf', width=15, height=15, arrangeGrob(grobs = plotList, ncol = 2))