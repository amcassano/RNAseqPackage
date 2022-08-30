#' Filter and Rank Significant Differentially Expressed Genes
#'
#' @description This function \itemize{
#' \item keeps only differentially expressed genes that are significant
#' \item adds an indicator of if genes in the experimental group are upregulated or downregulated compared to baseline
#' \item and calculates a change metric which will combine direction, Log2 Fold Change, and adjusted Pvalue}
#'
#' @param result dataframe contains \itemize{
#' \item gene annotations
#' \item log2 FC and adjusted p value}
#' @param padj_cutoff double, cutoff level for adjusted P vals to keep, defaults to 0.05
#' @param l2fc_cutoff double, cutoff level for log2 fold change, defaults to 0.8 (abs value)
#'
#' @return ranked, a sorted and filtered data frame containing \itemize{
#' \item annotations
#' \item change direction
#' \item fold change
#' \item adjusted p value
#' \item change metric (log2 fold change * - log10 padj)}
#' @export
#'
#' @examples
signif_deg <- function(result, padj_cutoff = 0.05, l2fc_cutoff = 0.8) {

  # create new df, removing any samples with NA as the value for FC or Pvals
  ranked <- tidyr::drop_na(result, c(log2FoldChange, padj))
  ranked <- tidyr::drop_na(result, c(log2FoldChange, padj))

  # keep only results that have an adjusted pvalue of less than 0.05 and FC more than 0.8
  ranked <- ranked %>%
    dplyr::filter(padj <= padj_cutoff) %>%
    dplyr::filter(abs(log2FoldChange) >= l2fc_cutoff) %>%
    dplyr::filter(MGI_Symbol != "")

  # add the change direction and change metric
  ranked <- ranked %>%
    dplyr::mutate(change_dir = ifelse(ranked$log2FoldChange < 0, "Down", "Up")) %>%
    dplyr::mutate(change_metric = ranked$log2FoldChange * -log10(ranked$padj))

  # keep only the relevant columns & arrange in order of change metric
  rownames(ranked) <- ranked$GeneID
  ranked <- ranked %>%
    dplyr::select(MGI_Symbol,
                  MGI_Desc,
                  GeneType,
                  change_dir,
                  log2FoldChange,
                  padj,
                  change_metric) %>%
    dplyr::arrange(desc(abs(change_metric)))

  # rename the columns
  ranked <- ranked %>%
    dplyr::rename("Change Direction" = "change_dir") %>%
    dplyr::rename("Adj. P Value" = "padj") %>%
    dplyr::rename("Log2 Fold Change" = "log2FoldChange") %>%
    dplyr::rename("Change Metric" = "change_metric")

  return(ranked)
}


#' Export Gene List
#'
#' @param tab data frame, contains ranked DEGs, the output of signif_deg
#' @param filename string, name to save the exported file as
#' @param direction string ("up", "down", or NULL), which genes to export
#' @param saveas string, extention for filename, defaults to ".csv"
#'
#' @return CSV file exported to CWD
#' @export
#'
#' @examples
export_genelist <- function(tab, filename, direction = NULL, saveas = ".csv") {

  # if direction is "up", export only those genes with a + change metric
  if(is.null(direction)){
    write.csv(tab, file = BiocGenerics::paste(filename, saveas, sep = ""))
  }
  else if(stringi::stri_cmp_equiv(direction, "up", strength = 1)){
    upreg <- tab %>%
      dplyr::filter(`Change Metric` > 0) %>%
      dplyr::select(-`Change Direction`)

    utils::write.csv(upreg, file = BiocGenerics::paste(filename, "_upregulated", saveas, sep = ""))
  }
  # if direction is "down", export only genes with - change metric
  else if (stringi::stri_cmp_equiv(direction, "down", strength = 1)) {
    downreg <- tab %>%
      dplyr::filter(`Change Metric` < 0) %>%
      dplyr::select(-`Change Direction`)

    utils::write.csv(downreg, file = BiocGenerics::paste(filename, "_downregulated", saveas, sep = ""))
  }
  # if no direction is given export all genes
  else {
    utils::write.csv(tab, file = paste(filename, saveas, sep = ""))
  }
}
