#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_figure5_plots.R
# Function: create dummy plots to explain the mechanism in figure 5
#
############################################################################################################################

####################
# libraries        #
####################

library(ggplot2)
library(cowplot)


####################
# Functions        #
####################


####################
# Main code        #
####################

# the ERAP levels in non-risk
nr_erap_df <- data.frame(
  x = factor(c(rep('NF-KB\n100%', times = 20), rep('50/50', times = 10), rep('CREL\n100%', times = 10)), levels = c('NF-KB\n100%', '50/50', 'CREL\n100%')), 
  erap2 = c(seq(1.5, 1.6, length.out = 20), 
            seq(0.75, 0.85, length.out = 10), 
            seq(0, 0.1, length.out = 10))
)
nr_erap_df_p <- ggplot(
  data = nr_erap_df, 
  mapping = aes(x = x, y = erap2)
)  + geom_jitter() + 
  ylim(c(0,2.2)) +
  xlab('') + ylab('ERAP2') +
  theme(axis.text.y=element_blank(), axis.ticks.y = element_blank()) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# the ERAP levels in risk
ra_erap_df <- data.frame(
  x = factor(c(rep('NF-KB\n100%', times = 10), rep('50/50', times = 10), rep('CREL\n100%', times = 20)), levels = c('NF-KB\n100%', '50/50', 'CREL\n100%')), 
  erap2 = c(seq(2, 2.1, length.out = 10), 
            seq(1, 1.1, length.out = 10), 
            seq(0, 0.1, length.out = 20))
)
ra_erap_df_p <- ggplot(
  data = ra_erap_df, 
  mapping = aes(x = x, y = erap2)
)  + geom_jitter() + 
  ylim(c(0,2.2)) +
  xlab('') + ylab('ERAP2') +
  theme(axis.text.y=element_blank(), axis.ticks.y = element_blank()) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))

# the line going down for nfkb in non-risk
nr_nfkb_df <- data.frame(
  x = c(0, 0.5, 1), 
  nfkb = c(1, 0.5, 0)
)
nr_nfkb_df_p <- ggplot(
  data = nr_nfkb_df, mapping = aes(x = x, y = nfkb)
) + geom_area(fill = "skyblue", alpha = 0.8) + 
  ylim(c(0,1.2)) +
  xlab('') + ylab('NF-KB') + 
  theme(axis.text.x=element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank()) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# the line going up for C-REL in non-risk
nr_crel_df <- data.frame(
  x = c(0, 0.5, 1), 
  crel = c(0, 0.5, 1)
)
nr_crel_df_p <- ggplot(
  data = nr_crel_df, mapping = aes(x = x, y = crel)
) + geom_area(fill = "#EE4B2B", alpha = 0.8) + 
  xlab('') + ylab('CREL') + 
  ylim(c(0,1.2)) + 
  theme(axis.text.x=element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank()) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))

# the line going dow for the nfkb in risk
ra_nfkb_df <- data.frame(
  x = c(0, 0.5, 1), 
  nfkb = c(.8, 0.4, 0)
)
ra_nfkb_df_p <- ggplot(
  data = ra_nfkb_df, mapping = aes(x = x, y = nfkb)
) + geom_area(fill = "skyblue", alpha = 0.8) + 
  xlab('') + ylab('NF-KB') + 
  ylim(c(0,1.2)) + 
  theme(axis.text.x=element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank()) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# the line going up for C-REL in risk
ra_crel_df <- data.frame(
  x = c(0, 0.5, 1), 
  crel = c(0, 0.6, 1.2)
)
ra_crel_df_p <- ggplot(
  data = ra_crel_df, mapping = aes(x = x, y = crel)
) + geom_area(fill = "#EE4B2B", alpha = 0.8) + 
  xlab('') + ylab('CREL') + 
  ylim(c(0,1.2)) + 
  theme(axis.text.x=element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank()) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# all together
p_nr_vs_ra <- plot_grid(
  plot_grid(
    ggdraw() + draw_label("Non-risk allele (T)", fontface = "bold", size = 16), 
    nr_erap_df_p, 
    nr_nfkb_df_p, 
    nr_crel_df_p, 
    nrow = 4, 
    ncol = 1, 
    rel_heights = c(1,3,3,3)
  ), 
  plot_grid(
    ggdraw() + draw_label("Risk allele (G)", fontface = "bold", size = 16), 
    ra_erap_df_p, 
    ra_nfkb_df_p, 
    ra_crel_df_p, 
    nrow = 4, 
    ncol = 1, 
    rel_heights = c(1,3,3,3)
  ), 
  nrow = 1, ncol = 2
)
# show the plot
p_nr_vs_ra
# save the plot
ggsave(paste0('~/multiome/plots/mo_figure_5_nfkb_rel.pdf'), plot = p_nr_vs_ra, width = 5, height = 5)

