
#' @title meb_QC_plots
#'
#' @description Takes a Seurat object and produces plots for basic QC metrics: RNA counts, feature counts, % mito, log10 genes per UMI
#' @param seurat_obj Seurat object
#'
#' @return grid of ggplot objects of QC metrics
#' @export
#'


meb_QC_plots <- function(seurat_obj) {

  require(Seurat)
  require(tidyverse)
  require(magrittr)
  require(data.table)
  require(gridExtra)

  #add QC metadata to Seurat objects
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

  #collect metadata
  metadata <- as.data.table(seurat_obj@meta.data)
  vars <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "log10GenesPerUMI")
  metadata <- metadata[, ..vars]

  #produce plots
  plot.list <- metadata[,2:5] %>% imap( ~ {

    give.n <- function(x) {
      return(c(y=median(x) + max (x)/10, label = round(median(x), 2)))
    }

    colors.ls <- setNames(
      c('darkolivegreen3', 'chocolate1', 'cornflowerblue', 'darkgoldenrod1'),
      c("nCount_RNA", "nFeature_RNA", "percent.mt", "log10GenesperUMI")
    )

    plot <- ggplot(data = metadata, aes(x = orig.ident, y = .x)) +
      geom_violin(trim = FALSE, fill = colors.ls[.y])+
      ggtitle(label = .y) + ylab(label = .y) +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank()
      ) +
      theme(
        plot.title = element_text(size = 12, hjust = 0.75),
        axis.line = element_line(colour = "black", size = 0.35),
        axis.title.x = element_blank(),
        legend.position = 'none'
      ) +
      stat_summary(fun = median, geom = "point", col="black") +
      stat_summary(fun.data = give.n,
                   geom = "text",
                   col = "black")


  })

  return_plots <- grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]], ncol = 2)

  return(return_plots)
}
