#' Rename Columns
#'
#' @description  rename columns, this allows for consistent naming and easier access to data
#'
#' @param colstr string, column to be renamed
#'
#' @return string, edited column name
#' @export
#'
#' @examples
rename_columns <- function(colstr){

  colstr <- BiocGenerics::paste(colstr)

  if (stringi::stri_cmp_equiv(colstr, "ensembl_gene_id", strength = 1)) {
    return("GeneID")
  }

  if(stringi::stri_cmp_equiv(colstr, "mgi_symbol", strength = 1)){
    return("MGI_Symbol")
  }

  else if(stringi::stri_cmp_equiv(colstr, "mgi_description", strength = 1)){
    return("MGI_Desc")
  }

  else if(stringi::stri_cmp_equiv(colstr, "gene_biotype", strength = 1)){
    return("GeneType")
  }

  else if(stringi::stri_cmp_equiv(colstr, "entrezgene_id", strength = 1)){
    return("EntrezID")
  }

  else if(stringi::stri_cmp_equiv(colstr, "go_id", strength = 1)){
    return("GO_ID")
  }

  else{
    return(colstr)
  }
}


#' Annotate Biomart
#'
#' add human readable gene IDs, descriptions, and gene biotypes to the data
#'
#' @param res DESeqResults object, must have ensembl gene ID in 'GeneID' column
#' @param gm biomart genemap, contains GeneID and MGI_symbol columns
#'
#' @return results object now with annotations
#' @export
#'
#' @examples
annotate_biomart <- function(res, gm) {

  # join the annotations with the results
  res <- res %>% dplyr::left_join(gm, by = c("GeneID" = "ensembl_gene_id"))

  # rename the columns
  colList <- colnames(res)
  newCols <- purr::map(colList, rename_columns)
  colnames(res) <- NULL
  colnames(res) <- newCols

  res <- res %>%
    dplyr::relocate(GeneID, MGI_Symbol)

  # return annotated results table
  return(res)
}
