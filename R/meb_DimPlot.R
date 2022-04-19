
#' @title meb_DimPlot
#'
#' @description Create Seurat DimPlot with consistent and attractive formatting
#'
#' @param seurat_obj Seurat object
#' @param pt_size point size - default 1.25
#' @param reduction
#' @param group_by metadata to group by - default "orig.ident"
#' @param cols color palette
#' @param dims_plot
#' @param split_by
#' @param label add labels to plot
#'
#' @return a ggplot object
#' @export

meb_DimPlot <- function(seurat_obj, pt_size = 1.25, reduction = "umap", group_by = "orig.ident", cols = NULL, dims_plot = 1:2, split_by = NULL, label = F) {

  require(Seurat)
  require(cowplot)
  require(tidyverse)
  require(MetBrewer)
  require(dplyr)

  #set x and y axis labels
  xlab <- paste(toupper(reduction), "1")
  ylab <- paste(toupper(reduction), "2")

  #set color map: tol for datasets with 12 groups, alphabet1 for datasets with more than 12 groups (but less than 22)
  n_groups <- dplyr::select(seurat_obj@meta.data, c(group_by)) %>% unique() %>% count()
  cols = met.brewer("Renoir", n=n_groups[,])

  #produce plots
  if (label == T){
    plot <- Seurat::DimPlot(seurat_obj,
                            group.by = group_by,
                            pt.size = pt_size,
                            cols = cols,
                            label = T,
                            reduction = reduction) +
      xlab(xlab) +
      ylab(ylab) +
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_line(size = 1.5),
            axis.title = element_text(size = 12)) +
      NoLegend()
  } else {
    plot <- Seurat::DimPlot(seurat_obj,
                            group.by = group_by,
                            pt.size = pt_size,
                            cols = cols,
                            label = F) +
      xlab(xlab) +
      ylab(ylab) +
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.line = element_line(size = 1.5),
            axis.title = element_text(size = 12))
  }

  return(plot)
}

