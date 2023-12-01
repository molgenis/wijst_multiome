#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_plot_markers.R
# Function: merge the seurat objects
############################################################################################################################

####################
# libraries        #
####################

# required to create object
library(Seurat)
# for plots
library(ggplot2)
library(cowplot)


####################
# Functions        #
####################


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')


####################
# Main Code        #
####################

# location of where to place the objects
seurat_objects_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/seurat_preprocess_samples//objects/'
object_all_cluster_filtered_ctd_cond_loc <- paste(seurat_objects_loc, 'mo_all_souped_clus_filtered_ctd_cond_20231129.rds', sep = '')

# read object
object_all <- readRDS(object_all_cluster_filtered_ctd_cond_loc)

# plot CD4T
plot_grid(
  plot_grid(
    ggdraw() + draw_label('CD4+ T general', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CD3D')), FeaturePlot(object_all, features = ('CD3E')), FeaturePlot(object_all, features = ('CD3G')), nrow = 1
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('CD4+ T Naive', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CCR7')), FeaturePlot(object_all, features = ('SELL')), NULL, nrow = 1, ncol = 3
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('CD4+ T Memory', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('LTB')), FeaturePlot(object_all, features = ('IL7R')), FeaturePlot(object_all, features = ('S100A4')), nrow = 1
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('CD4+ T Absent', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CD8A')), FeaturePlot(object_all, features = ('CD8B')), NULL, nrow = 1, ncol = 3
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  nrow = 4
)


# plot CD8T
plot_grid(
  plot_grid(
    ggdraw() + draw_label('CD8+ T general', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CD3D')), FeaturePlot(object_all, features = ('CD3E')), FeaturePlot(object_all, features = ('CD3G')), FeaturePlot(object_all, features = ('CD8A')),
      FeaturePlot(object_all, features = ('CD8B')), FeaturePlot(object_all, features = ('GZMB')), FeaturePlot(object_all, features = ('PRF1')), NULL, nrow = 2, ncol = 4
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('CD8+ T Naive', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CCR7')), FeaturePlot(object_all, features = ('SELL')), NULL, NULL, nrow = 1, ncol = 4
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('CD8+ T Memory', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('IL7R')), FeaturePlot(object_all, features = ('S100A4')), NULL, NULL, nrow = 1, ncol = 4
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  nrow = 3,
  rel_heights = c(2, 1, 1)
)


# plot NK
plot_grid(
  plot_grid(
    ggdraw() + draw_label('NK general', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('NKG7')), FeaturePlot(object_all, features = ('GNLY')), NULL, nrow = 1, ncol = 3
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('NKdim', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('FCGR3A')), FeaturePlot(object_all, features = ('GZMB')), FeaturePlot(object_all, features = ('PRF1')),
      nrow = 1
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('NKbright', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('KLRC1')), NULL, NULL, nrow = 1, ncol = 3
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('NK Absent', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CD8A')), FeaturePlot(object_all, features = ('CD8B')), NULL, nrow = 1, ncol = 3
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('NKbright Absent', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('FCGR3A')), FeaturePlot(object_all, features = ('GZMB')), FeaturePlot(object_all, features = ('PRF1')), nrow = 1
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  nrow = 5
)


# plot monocytes
plot_grid(
  plot_grid(
    ggdraw() + draw_label('monocyte general', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CD14')), NULL, NULL, nrow = 1, ncol = 3
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('cMono', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('LYZ')), FeaturePlot(object_all, features = ('S100A9')), FeaturePlot(object_all, features = ('CSF3R')),
      nrow = 1
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('ncMono', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('FCGR3A')), FeaturePlot(object_all, features = ('LYN')), FeaturePlot(object_all, features = ('CSF1R')), 
      FeaturePlot(object_all, features = ('IFITM1')), FeaturePlot(object_all, features = ('IFITM2')), FeaturePlot(object_all, features = ('IFITM3')), 
      nrow = 2
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('cMono Absent', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('FCGR3A')), FeaturePlot(object_all, features = ('LYN')), FeaturePlot(object_all, features = ('CSF1R')), 
      FeaturePlot(object_all, features = ('IFITM1')), FeaturePlot(object_all, features = ('IFITM2')), FeaturePlot(object_all, features = ('IFITM3')), 
      nrow = 2
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('ncMono Absent', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('LYZ')), FeaturePlot(object_all, features = ('S100A9')), FeaturePlot(object_all, features = ('CSF3R')), nrow = 1
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  nrow = 5,
  rel_heights = c(1, 1, 2, 2, 1)
)


# plot B
plot_grid(
  plot_grid(
    ggdraw() + draw_label('B general', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CD79A')), nrow = 1, ncol = 1
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('B non-plasma', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('MS4A1')), nrow = 1, ncol = 1
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  nrow = 2
)


# plot DC
plot_grid(
  plot_grid(
    ggdraw() + draw_label('mDC', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CD1C')), FeaturePlot(object_all, features = ('ITGAX')),
      nrow = 1, ncol = 3
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('pDC', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CLEC4C')), NULL,
      nrow = 1, ncol = 2
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('DC absent', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CD14')), NULL,
      nrow = 1, ncol = 2
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('mDC Absent', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CLEC4C')), NULL, 
      nrow = 1, ncol = 2
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  plot_grid(
    ggdraw() + draw_label('pDC Absent', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('CD1C')), FeaturePlot(object_all, features = ('ITGAX')), 
      nrow = 1, ncol = 2
    ),
    rel_heights=c(0.1, 1), nrow = 2
  ),
  nrow = 5
)


# plot megakaryocytes
plot_grid(
  plot_grid(
    ggdraw() + draw_label('megakaryocyte general', fontface='bold'),
    plot_grid(
      FeaturePlot(object_all, features = ('GP9')), FeaturePlot(object_all, features = ('ITGA2B')), 
      FeaturePlot(object_all, features = ('PF4')), FeaturePlot(object_all, features = ('PPBP')), 
      nrow = 2, ncol = 2
    ),
    rel_heights=c(0.1, 1), nrow = 2
  )
)
