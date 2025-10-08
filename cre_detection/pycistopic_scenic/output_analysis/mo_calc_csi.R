#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_calc_csi.R
# Function: calculate CSI for regulons and plot them
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(HiClimR)
library(ggplot2)
library(ggdendro)
library(patchwork)


####################
# Functions        #
####################


calc_cors <- function(x,
                      rownames_col = NULL,
                      nSplit = 1
                      ) {
  # convert data.table to data.frame
  if (inherits(x, "data.table")) {
    x <- as.data.frame(x)
  }
  
  # further convert dataframe
  if (is.data.frame(x)) {
    # use a column with the rownames if supplied
    if (!is.null(rownames_col)) {
      # warn if rownames already exist
      if (!is.null(rownames(x)) && any(rownames(x) != seq_len(nrow(x)))) {
        warning("Existing rownames will be replaced by values from 'rownames_col'.")
      }
      # extract that column
      rn <- x[[rownames_col]]
      # remove it from the original table
      x[[rownames_col]] <- NULL
      # then convert into a regular matrix
      x <- as.matrix(x)
      # and set that as the rownames of the tale
      rownames(x) <- rn
    }
    # don't bother with the rowname that is a column name otherwise
    else {
      x <- as.matrix(x)
    }
  }
  
  # if a sparse matrix convert
  if (inherits(x, "sparseMatrix")) {
    x <- as.matrix(x)
  }
  # and make sure that the data was numeric to start with
  if (!is.numeric(x)) {
    stop("Input must be numeric after conversion.")
  }
  
  # now compute the actual correlations
  cor_mat <- HiClimR::fastCor(
    t(x),
    nSplit = nSplit,
    upperTri = F
  )
  
  # add back the rownames if we had them
  if (!is.null(rownames(x))) {
    rownames(cor_mat) <- colnames(cor_mat) <- rownames(x)
  }
  
  return(cor_mat)
}


calc_csi <- function(cor_mat) {
  # extract the featurs the matrix is computed on
  features <- colnames(cor_mat)
  # get how many features these were
  n <- ncol(cor_mat)
  
  # make a list to store the CSIs
  results <- vector("list", n * n)
  # set up the counter
  idx <- 1
  # check each feature
  for (i in seq_len(n)) {
    # against each other feature
    for (j in i:n) {
      # extract the correlation of this pair
      r_ij <- cor_mat[i, j]
      # extract all the correlations for the featurs in this pair
      cor_i <- cor_mat[i, ]
      cor_j <- cor_mat[j, ]
      
      # calculate the fraction of values that are smaller or equal
      fraction_lower <- mean((cor_i <= r_ij) & (cor_j <= r_ij))
      
      # store the CSI of (i,j)
      results[[idx]] <- list(
        # with the names
        feature_1 = features[i],
        feature_2 = features[j],
        CSI = fraction_lower
      )
      # increase the counter
      idx <- idx + 1
      # when i and j are not the same, we also need to store (j,i)
      if (i != j) {
        # in the list again
        results[[idx]] <- list(
          feature_1 = features[j],
          feature_2 = features[i],
          CSI = fraction_lower
        )
        # increase the counter
        idx <- idx + 1
      }
    }
  }
  
  # if we made a list larger than needed, remove what we didn't need
  results <- results[seq_len(idx - 1)]
  # make a big dataframe of all of the lists and their values
  results <- do.call(rbind, lapply(results, as.data.frame))
  
  # set the features to be characters
  results$feature_1 <- as.character(results$feature_1)
  results$feature_2 <- as.character(results$feature_2)
  # and ensure the CSI is a numerical value
  results$CSI <- as.numeric(results$CSI)
  
  return(results)
}

csi_long_to_matrix <- function(csi_long) {
  # reshape into wide matrix
  mat <- reshape2::acast(csi_long, feature_1 ~ feature_2, value.var = "CSI")
  
  # enforce symmetry (in case only upper triangle was stored)
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  
  return(mat)
}


clust_csi <- function(csi,
                      method = "complete",
                      dist_fun = function(x) stats::dist(1 - x), 
                      k = NULL,
                      h = NULL) {
  # detect input type
  if (is.data.frame(csi)) {
    if (!all(c("feature_1","feature_2","CSI") %in% colnames(csi))) {
      stop("Long-format input must have columns: feature_1, feature_2, CSI")
    }
    csi_mat <- csi_long_to_matrix(csi)
  } else if (is.matrix(csi)) {
    csi_mat <- csi
  } else {
    stop("Input must be either a correlation-like matrix or a long-format data.frame")
  }
  
  # check symmetry
  if (!all(csi_mat == t(csi_mat))) {
    warning("Input matrix is not perfectly symmetric. Forcing symmetry.")
    csi_mat[lower.tri(csi_mat)] <- t(csi_mat)[lower.tri(csi_mat)]
  }
  
  # compute distance
  dist_mat <- dist_fun(csi_mat)
  
  # hierarchical clustering
  hc <- hclust(dist_mat, method = method)
  ordered_features <- hc$labels[hc$order]
  
  # convert to long format with ordered factors
  long <- reshape2::melt(csi_mat, varnames = c("feature_1", "feature_2"), value.name = "CSI")
  long$feature_1 <- factor(long$feature_1, levels = ordered_features)
  long$feature_2 <- factor(long$feature_2, levels = ordered_features)
  
  # optional cluster cutting
  clusters <- NULL
  if (!is.null(k)) {
    clusters <- cutree(hc, k = k)
  } else if (!is.null(h)) {
    clusters <- cutree(hc, h = h)
  }
  # put the results in a list
  result_list <- list(
    hc = hc,
    order = ordered_features,
    dist_mat = dist_mat,
    long = long, 
    clusters = clusters
  )
  
  return(result_list)
}



plot_csi_heatmap <- function(clust_result,
                             palette = c("navy", "white", "firebrick3"), 
                             cluster_palette = NULL, 
                             bar_height = 1) {
  if (!all(c("hc","long","order") %in% names(clust_result))) {
    stop("Input must be the result of clust_csi()")
  }
  
  # extract the hierchial clustering
  hc <- clust_result$hc
  # order the features based on cluster
  ordered_features <- clust_result$order
  # extract the long-format ggplot2-compatible table
  csi_long <- clust_result$long
  
  # make the dendogram 
  dendro_data <- ggdendro::dendro_data(hc, type = "rectangle")
  # add the dendotram to a ggplot, making it into a segment
  dendro_plot <- ggplot(segment(dendro_data)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0)) +
    # set expand so the two plots align
    scale_x_continuous(expand = c(0,0), limits = c(0.5, length(ordered_features) + 0.5))
  
  # grab clusters from object
  clusters <- clust_result$clusters
  
  # if clusters exist, add colored bars
  if (!is.null(clusters)) {
    # make a df to map the features to the clusters
    cluster_df <- data.frame(label = names(clusters),
                             cluster = factor(clusters))
    # add the dendogram data
    cluster_df <- merge(cluster_df, dendro_data$labels[,c("label", "x")], by =" label")
    
    # collect the palette if not supplied
    if (is.null(cluster_palette)) {
      cluster_palette <- scales::hue_pal()(length(unique(clusters)))
    }
    # make the dendogram plot, adding the clusters
    dendro_plot <- dendro_plot +
      # as a one-row tile
      geom_tile(data = cluster_df,
                aes(x = x, y = -bar_height/2, fill = cluster),
                width = 1, height = bar_height) +
      # with the correct palette
      scale_fill_manual(values = cluster_palette) +
      # no legend because the actual clusters doesn't really matter in the plot
      theme(legend.position = "none")
  }
  
  # make the heatmap
  heatmap_plot <- ggplot(csi_long, aes(x = feature_1, y = feature_2, fill = CSI)) +
    # which are just tiles
    geom_tile() +
    # set expand so the two plots align
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_gradientn(colors = palette) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.title = element_blank(),
          plot.margin = margin(0,0,0,0))
  
  # merge the two plots
  combined <- dendro_plot / heatmap_plot + plot_layout(heights = c(1, 4))
  return(combined)
}

####################
# Main code        #
####################
