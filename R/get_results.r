#' Get DE results
#'
#' Get the results of a deseq object
#'
#' @param contrast_str string, defaults to "Condition", the comparison of interest for pairwise comparison
#' @param numerator string, name of the conditon of interest for comparison
#' @param denominator string, name of baseline condition for comparison
#' @param tidy_result boolean, defaults to TRUE, dictates if the results are returned as a dataframe or obj
#' @param deseq_obj data
#' @param shrinkage boolean, defaults to TRUE, dictates if lfc shrinking formula is used
#'
#' @return DESeq Results object, or a dataframe
#' @export
get_DEresults <- function(contrast_str = "Condition", numerator, denominator, tidy_result = TRUE, deseq_obj, shrinkage = TRUE){
  con <- c(contrast_str, numerator, denominator)

  if(shrinkage){
    results <- DESeq2::lfcShrink(dds = deseq_obj,
                                contrast = con,
                                type = "ashr")

    if(tidy_result){
      results <- as.data.frame(results)
      results <- tibble::rownames_to_column(results, var = "GeneID")
    }
  }
  else{
    results <- DESeq2::results(deseq_obj, con, tidy = tidy_result)

    if(tidy_result){
    results <- dplyr::rename(results, "GeneID" = "row")
    }

  }

  return(results)

}