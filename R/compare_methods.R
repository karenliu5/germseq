#' Compare methods of microbial differential abundance analysis (DAA)
#'
#' `compare_DAA_methods()` applies three popular methods of DAA to a phyloseq object
#' and returns a table with all three results. Note that this only works with R4.0.5.
#'
#' @param ps a phyloseq object containing the data that is to be analyzed.
#' @param group a string name of a binary variable to compare.
#' @param prev.thr the prevalence threshold, between (0, 1]. Default value is 0.1
#'
#' @returns A data frame containing whether each method found a statistically
#' significant difference between the groups for each taxon.
#'
#' @examples
#' library(microbiome)
#' data(atlas1006)
#' atlas1006 <- phyloseq::subset_samples(atlas1006, (bmi_group == "lean" | bmi_group == "obese"))
#' compare_DAA_methods(ps = atlas1006, group = "bmi_group", prev.thr = 0.1)
#'
#' @export
#' @importFrom magrittr "%>%"

compare_DAA_methods <- function(ps, group, prev.thr = 0.1){
  # Check prev.thr
  if(prev.thr >= 1 | prev.thr < 0){
    stop("Error: prev.thr must be a number between [0, 1).")
  }

  # Check for Phyloseq object

  # Check if group is valid
  if(! assertthat::is.string(group)){ # Group is a string
    stop("Error: group must be a string.")
  }

  col_names <- colnames(phyloseq::sample_data(ps))
  if(! (group %in% col_names)){
    stop("Error: group could not be found in phyloseq object.")
  } else {
    ind <- which(group == col_names)
  }

  #Automatically remove NAs
  num_NA <- sum(is.na(phyloseq::sample_data(ps)[,ind]))
  if(num_NA > 0){
    var_values <- phyloseq::sample_data(ps)[[group]]
    ps <- phyloseq::prune_samples(!is.na(var_values), ps)
    warning("Warning: Rows with missing values were removed")
  }

  # Check for two groups
  if(nrow(unique(phyloseq::sample_data(ps)[,ind])) != 2){
    stop("Error: Make sure group takes only two values.")
  }

  # Preprocess Data
  prevdf <- microbiome::prevalence(ps)

  mask <- as.logical(prevdf > prev.thr)
  keepTaxa <- names(prevdf[mask])

  ps_filt <- phyloseq::prune_taxa(keepTaxa, ps)

  # DESeq2
  form <- stats::as.formula(paste("~", group))
  deseq2_format <- phyloseq::phyloseq_to_deseq2(ps_filt, form)
  deseq2_format <- DESeq2::DESeq(deseq2_format, test = "Wald", fit="parametric")
  deseq2_res <- DESeq2::results(deseq2_format)

  # ANCOMBC
  ancombc_form <- ANCOMBC::ancombc(
    phyloseq = ps_filt,
    formula = group,
    p_adj_method = "fdr",
    zero_cut = 1.0,
    lib_cut = 0,
    group = group,
    struc_zero = TRUE,
    neg_lb = TRUE,
    tol = 1e-5,
    max_iter = 100,
    conserve = TRUE,
    alpha = 0.05,
    global = TRUE
  )

  # ALDEx2
  aldex2_da <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps_filt)),
                             phyloseq::sample_data(ps_filt)[[group]],
                             test = "t", effect = TRUE,
                             denom = "iqlr")

  summ <- dplyr::full_join(data.frame(taxon = row.names(aldex2_da), aldex2 = (aldex2_da$wi.eBH < 0.05)),
                           data.frame(taxon = row.names(ancombc_form$res$diff_abn), ancombc = ancombc_form$res$diff_abn[,1]), by = "taxon") %>%
    dplyr::full_join(data.frame(taxon = row.names(deseq2_res), deseq2 = deseq2_res$padj < 0.05), by="taxon")

  return(summ)
}



#' Compare methods of microbial differential abundance analysis (DAA)
#'
#' `compare_DAA_methods_2()` applies three popular methods of DAA to a phyloseq object
#' and returns a table with all three results. Updated for R4.2.2.
#'
#' @param ps a phyloseq object containing the data that is to be analyzed.
#' @param group a string name of a binary variable to compare.
#' @param prev.thr the prevalence threshold, between (0, 1]. Default value is 0.1
#'
#' @returns A data frame containing whether each method found a statistically
#' significant difference between the groups for each taxon.
#'
#' @examples
#' library(microbiome)
#' data(atlas1006)
#' atlas1006 <- phyloseq::subset_samples(atlas1006, (bmi_group == "lean" | bmi_group == "obese"))
#' compare_DAA_methods_2(ps = atlas1006, group = "bmi_group", prev.thr = 0.1)
#'
#' @export
#' @importFrom magrittr "%>%"

compare_DAA_methods_2 <- function(ps, group, prev.thr = 0.1){
  # Check prev.thr
  if(prev.thr >= 1 | prev.thr < 0){
    stop("Error: prev.thr must be a number between [0, 1).")
  }

  # Check for Phyloseq object

  # Check if group is valid
  if(! assertthat::is.string(group)){ # Group is a string
    stop("Error: group must be a string.")
  }

  col_names <- colnames(phyloseq::sample_data(ps))
  if(! (group %in% col_names)){
    stop("Error: group could not be found in phyloseq object.")
  } else {
    ind <- which(group == col_names)
  }

  #Automatically remove NAs
  num_NA <- sum(is.na(phyloseq::sample_data(ps)[,ind]))
  if(num_NA > 0){
    var_values <- phyloseq::sample_data(ps)[[group]]
    ps <- phyloseq::prune_samples(!is.na(var_values), ps)
    warning("Warning: Rows with missing values were removed")
  }

  # Check for two groups
  if(nrow(unique(phyloseq::sample_data(ps)[,ind])) != 2){
    stop("Error: Make sure group takes only two values.")
  }

  # Preprocess Data
  prevdf <- microbiome::prevalence(ps)

  mask <- as.logical(prevdf > prev.thr)
  keepTaxa <- names(prevdf[mask])

  ps_filt <- phyloseq::prune_taxa(keepTaxa, ps)

  # DESeq2
  form <- stats::as.formula(paste("~", group))
  deseq2_format <- phyloseq::phyloseq_to_deseq2(ps_filt, form)
  deseq2_format <- DESeq2::DESeq(deseq2_format, test = "Wald", fit="parametric")
  deseq2_res <- DESeq2::results(deseq2_format)

  # ANCOMBC
  ancombc_form <- ANCOMBC::ancombc(
    data = ps_filt,
    formula = group,
    p_adj_method = "fdr",
    prv_cut = 0,
    lib_cut = 0,
    group = group,
    struc_zero = TRUE,
    neg_lb = TRUE,
    tol = 1e-5,
    max_iter = 100,
    conserve = TRUE,
    alpha = 0.05,
    global = TRUE
  )

  # ALDEx2
  aldex2_da <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps_filt)),
                             phyloseq::sample_data(ps_filt)[[group]],
                             test = "t", effect = TRUE,
                             denom = "iqlr")

  summ <- dplyr::full_join(data.frame(taxon = row.names(aldex2_da), aldex2 = (aldex2_da$wi.eBH < 0.05)),
                           data.frame(taxon = ancombc_form$res$diff_abn$taxon, ancombc = ancombc_form$res$diff_abn[,3]), by = "taxon") %>%
    dplyr::full_join(data.frame(taxon = row.names(deseq2_res), deseq2 = deseq2_res$padj < 0.05), by="taxon")

  return(summ)
}


