#' Atlas1006 Output Data
#'
#' Output from running Atlas1006 data through compare_DAA_methods()
#'
#' @format ## `atlas1006_output`
#' A data frame with 115 rows and 5 columns:
#' \describe{
#'   \item{taxon}{Taxon name}
#'   \item{aldex2, ancombc, deseq2}{Whether each method found a significant DAA}
#'   \item{score}{Number of methods that found the taxon to be significant}
#'   ...
#' }
#' @source library(microbiome)
"atlas1006_output"
