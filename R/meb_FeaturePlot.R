
#' @title meb_FeaturePlot
#'
#' @description Create Seurat FeaturePlot with consistent and attractive formatting
#'
#'
#' @param seurat_obj Seurat object
#' @param features genes to plot
#' @param pt_size point size - default 1.25
#' @param reduction
#' @param scale_data
#' @param order
#'
#' @return a cowplot gtable with ggplot objects
#' @export


meb_FeaturePlot <- function(seurat_obj, features, pt_size = 1.25, reduction = "umap", scale_data = F, order = T) {
  require(Seurat)
  require(cowplot)
  require(tidyverse)
  require(MetBrewer)
  require(dplyr)

  #set x and y axis labels
  xlab <- paste(toupper(reduction), "1")
  ylab <- paste(toupper(reduction), "2")

  cols <- met.brewer("Troy",n=15,type="continuous")

  slot <- ifelse(scale_data == T, "scale.data", "data")

  plots <- Seurat::FeaturePlot(seurat_obj,
                               features = features,
                               pt.size = pt_size,
                               combine = F,
                               reduction = reduction,
                               slot = slot,
                               order = order)
  plots <- lapply(plots, function(x) x +
                    xlab(xlab) +
                    ylab(ylab) +
                    scale_colour_gradientn(colours = cols) +
                    theme(axis.ticks = element_blank(),
                          axis.text = element_blank(),
                          axis.title = element_text(size = 20),
                          axis.line = element_line(size = 1)))
  return_plots <- cowplot::plot_grid(plotlist = plots, ncol = 1)
  return(return_plots)
}


