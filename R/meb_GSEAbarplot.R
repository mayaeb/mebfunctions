
#' @title meb_GSEAbarplot
#'
#' @description Create consistently formatted bar plots from GSEA analysis output. Adapted from a similar function written by Miranda Hunter.
#'
#' @param result output from fgsea analysis
#' @param n_terms number of terms to plot
#' @param metric which metric to plot - options are "NES" or "pval"
#' @param fill option to fill bar. if F, bars will all be black
#' @param change_labels reformat pathway/component labels - removes "GOCC" or "GOBP" header, removes dashes, changes to lowercase
#' @param pval_fill fills bars based on p value
#'
#' @return a ggplot object
#' @export


meb_GSEAbarplot <- function(result, n_terms, metric = "NES", fill, change_labels = TRUE, pval_fill = TRUE) {

  require(pals)
  require(dplyr)
  require(MetBrewer)



  #set color scheme for filled bars
  if (fill == T) {
    cols_pval <- pals::viridis(n_terms)
    cols_NES <- pals::magma(n_terms)
  } else if (fill == F) {
    cols_pval <- "black"
    cols_NES <- "black"
  }

  #remove annotations and change to lower calse
  if (change_labels == TRUE) {
    result$pathway <- result$pathway %>%
      gsub(pattern = "GOBP_", replacement = "", x = ., ignore.case = F) %>%
      gsub(pattern = "GOCC_", replacement = "", x = ., ignore.case = F) %>%
      gsub(pattern = "_", replacement = " ", x = .) %>%
      gsub(pattern = "plus", replacement = "+", x = .) %>%
      tolower()
  }

  #plot p-value bar
  if (metric == "pval") {
    data <- result %>% filter(NES>0) %>% arrange(pval) %>% slice(., 1:n_terms)
    if (pval_fill) {
      data$pathway <- factor(data$pathway, levels = rev(data$pathway))
      plot <- ggplot(data, aes(pval, pathway, fill = pval)) +
        geom_bar(stat = "identity") +
        scale_fill_gradientn(colours = met.brewer("Tam",n=1000, type = "continuous")) +
        theme_minimal() +
        labs(x = "p-value",
             y = "GO Term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(plot)
    } else {
      data$pathway <- factor(data$pathway, levels = rev(data$pathway))
      plot <- ggplot(data, aes(x = pval, y = pathway)) +
        geom_bar(stat = "identity", fill = cols_pval) +
        theme_minimal() +
        labs(x = "p-value",
             y = "GO term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(plot)

    }
  } else if (metric == "NES") {
    data <- result %>% arrange(-NES) %>% slice(., 1:n_terms)
    if (pval_fill) {
      data$pathway <- factor(data$pathway, levels = rev(data$pathway))
      plot <- ggplot(data, aes(NES, pathway, fill = pval)) +
        geom_bar(stat = "identity") +
        scale_fill_gradientn(colours = met.brewer("Tam",n=1000, type = "continuous")) +
        theme_minimal() +
        labs(x = "Normalized Enrichment Score",
             y = "GO Term") +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(plot)
    } else {
      data$pathway <- factor(data$pathway, levels = rev(data$pathway))
      plot <- ggplot(data, aes(NES, pathway, fill = pval)) +
        geom_bar(stat = "identity", fill = cols_NES) +
        theme_minimal() +
        labs(x = "Normalized Enrichment Score",
             y = "GO Term") +
        theme_minimal() +
        theme(axis.title.x = element_text(size = 13, color = "black"),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.text.x = element_text(size = 12, color = "black"))
      return(plot)
    }

  } else {
    stop("Metric must be NES or p-value")
  }
}
