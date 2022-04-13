
#' @title meb_rename_ensemblids
#'
#' @description Translate zebrafish ensembl ids in a matrix to gene names
#'
#' @param data a scRNAseq data matrix with counts and zebrafish ensembl ids
#' @return a data matrix with ensembl ids translated into zebrafish gene names
#' @export
#'

meb_rename_ensemblids <- function(data) {

  #check to make sure these are the only dependencies
  require(org.Dr.eg.db)
  require(tidyverse)

  #wrangle matrix to get ensemblids for translation
  message('Creating data frame...')
  data <- as.data.frame(t(data))
  colnames(data) <- data[1,]
  data <- data[-1,]
  ensemblids <- rownames(data)

  #map ensemblids to zfish genes, add to matrix
  message('Mapping ENSEMBL IDs to fish genes...')
  symbols <- mapIds(org.Dr.eg.db, ensemblids, keytype="ENSEMBL", column="SYMBOL", multiVals = "asNA")
  data$geneid <- symbols
  data <- data[!is.na(data$geneid),]
  data<- data[!duplicated(data$geneid),]
  rownames(data) <- data$geneid

  message('Done!')
  return(data)
}
