#' Sample TFs data
#' 
#' This data set can directly be generated use the function scGRN_getTF and data set sample_interation
#' 
#' @format a data.table with 100 rows and 7 variables
#' \describe{
#'  \item{gene}{gene id}
#'  \item{promoter}{promoter region for target gene}
#'  \item{enhancer}{linked enhancer for a target gene}
#'  \item{promoter_TF}{possible transcription factors that would bind the promoter, its element type is list}
#'  \item{enhancer_TF}{possible transcription factors that bind to the enhancer, ts element type is list}
#' }
"sample_TFs"