
#' @title meb_translateEnsIDs_drerio
#'
#' @description Translate a vector of zebrafish EnsemblIDs to gene symbols. Replaces un-mapped IDs with "un-translated".
#'
#' @param ensemblIDs a vector of danio rerio EnsemblIDs
#' @return a matched vector danio rerio gene symbols
#' @export
#'

meb_translateEnsIDs_drerio <- function(ensemblIDs = NULL) {

  require(bioRmart)

  #make sure input is a vector
  if(!is.vector(ensemblIDs)) {
    message('Input is not a vector')
    return(ensemblIDs)
  }

  #get zfish mart
  ensembl <- useEnsembl(biomart = "genes")
  ensembl.con <- useMart("ensembl", dataset = "drerio_gene_ensembl")

  message('Translating...')

  #translate
  translated <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ensemblIDs,
                      mart = ensembl.con)

  #restructure to data frame for merging
  ensemblIDs <- ensemblIDs %>% as.data.frame()
  colnames(ensemblIDs) <- c("ensembl_gene_id")
  joined <- dplyr::left_join(ensemblIDs, translated, by = "ensembl_gene_id", all.y=TRUE)
  new_names <- joined$external_gene_name
  new_names[is.na(new_names)] = c("un-mapped")

  message('Done!')
  return(new_names)
}

