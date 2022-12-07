#' Atlas1006 Output Data
#'
#' Output from running Atlas1006 data through compare_DAA_methods()
#'
#' @format ## `atlas1006_output`
#' A data frame with 115 rows and 6 columns:
#' \describe{
#'   \item{taxon}{Taxon name}
#'   \item{aldex2, ancombc, deseq2}{Benjamini-Hochberg adjusted p-value
#'   resulting from the respective DAA analysis}
#'   \item{fisher_p}{Resulting p-value from the Fisher method meta-analysis of the aldex2,
#'   ancombc, and deseq2 adjusted p-values}
#'   \item{rawcount}{Number of methods that found the taxon to be significant}
#' }
#' @source library(microbiome)
"atlas1006_output"
