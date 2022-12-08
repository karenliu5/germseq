#' Compare methods of microbial differential abundance analysis (DAA) (R4.0.5)
#'
#' `compare_DAA_methods()` applies three popular methods of DAA to a phyloseq object
#' and returns a table with all three results. Note that this only works with R4.0.5.
#'
#' @param ps a phyloseq object containing the data that is to be analyzed.
#' @param group a string name of a binary variable to compare.
#' @param prevThr the prevalence threshold, between (0, 1]. Default value is 0.1
#'
#' @returns A data frame containing booleans indicating whether each method found
#' a statistically significant difference between the groups for each taxon after
#' adjusting for multiple hypotheses using FDR.
#'
#' @examples
#' library(microbiome)
#' data(atlas1006)
#' atlas1006 <- phyloseq::subset_samples(atlas1006, (bmi_group == "lean" | bmi_group == "obese"))
#' compare_DAA_methods(ps = atlas1006, group = "bmi_group", prevThr = 0.1)
#'
#'@references
#'Gloor, G., Fernandes, A., Macklaim, J., Albert, A., Links, M., Quinn, T., Wu, J.R.,
#'Wong, R.G., and B. Lieng (2013). ALDEx2: Analysis Of Differential Abundance Taking
#'Sample Variation Into Account. R package version 1.30.0.
#'\href{https://bioconductor.org/packages/release/bioc/html/ALDEx2.html}{Link}.
#'
#'Huang, L. (2020). ANCOMBC: Microbiome differential abudance and correlation
#'analyses with bias correction. R package version 2.0.1.
#'\href{https://www.bioconductor.org/packages/release/bioc/html/ANCOMBC.html}{Link}.
#'
#'Lahti, L. and S. Shetty (2019). microbiome: Microbiome Analytics. R package version
#'1.20.0. \href{https://www.bioconductor.org/packages/release/bioc/html/microbiome.html}{Link}.
#'
#'Lahti, L., Shetty, S., Ernst, F.M. et al. (2021). Orchestrating Microbiome Analysis
#'with Bioconductor. \href{microbiome.github.io/oma/}{Link}.
#'
#'Love, M., Huber, W., and S. Anders (2014). DESeq2: Differential gene expression
#'analysis based on the negative binomial distribution. R package version 1.38.0.
#'\href{https://bioconductor.org/packages/release/bioc/html/DESeq2.html}{Link}.
#'
#'McMurdle, P.J. and S. Holmes (2013). phyloseq: Handling and analysis of high-throughput
#'microbiome census data. R package version 1.42.0.
#'\href{https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html}{Link}.
#'
#'Olberding, N. (2019). Introduction to the Statistical Analysis of Microbiome Data in R.
#'\href{https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/}{Link}.
#'
#'Winkler, A.M. (2012). The logic of the Fisher method to combine P-values.
#'\href{https://brainder.org/2012/05/11/the-logic-of-the-fisher-method-to-combine-p-values/}{Link}.
#'
#'@export
#' @importFrom magrittr "%>%"
compare_DAA_methods <- function(ps, group, prevThr = 0.1){

  # Check if inputs are valid

  # Check if prevalence threshold, prevThr, is between 0 and 1.
  if (prevThr >= 1 | prevThr < 0){
    stop("Error: prevThr must be a number between [0, 1).")
  } else {;}

  # Check if group is a string
  if (! assertthat::is.string(group)){ # Group is a string
    stop("Error: group must be a string.")
  } else {;}

  # Check if group is a valid variable to pass
  col_names <- colnames(phyloseq::sample_data(ps))
  if (! (group %in% col_names)){
    stop("Error: group could not be found in phyloseq object.")
  } else {
    ind <- which(group == col_names)
  }

  # Automatically remove NAs
  num_NA <- sum(is.na(phyloseq::sample_data(ps)[ ,ind]))
  if (num_NA > 0){
    var_values <- phyloseq::sample_data(ps)[[group]]
    ps <- phyloseq::prune_samples(!is.na(var_values), ps)
    warning("Warning: Rows with missing values were removed")
  } else {;}

  # Check for two groups
  if (nrow(unique(phyloseq::sample_data(ps)[ ,ind])) != 2){
    stop("Error: Make sure group takes only two values.")
  } else {;}



  # Apply prevalence filter to data

  prevdf <- microbiome::prevalence(ps) # Get prevalence

  mask <- as.logical(prevdf > prevThr)
  keepTaxa <- names(prevdf[mask]) # Identify taxa with prevalence higher than threshold

  ps_filt <- phyloseq::prune_taxa(keepTaxa, ps) # Filter out taxa with prevalence lower than threshold


  # Run analyses on data

  message("Now running DESeq2")  # DESeq2

  form <- stats::as.formula(paste("~", group))

  deseq2_format <- phyloseq::phyloseq_to_deseq2(ps_filt, form) %>% suppressMessages()
  deseq2_format <- DESeq2::DESeq(deseq2_format, test = "Wald", fit="parametric") %>% suppressMessages()
  deseq2_out <- DESeq2::results(deseq2_format)


  message("Now running ANCOMBC")  # ANCOMBC

  ancom_taxon <- NULL
  ancom_q <- NULL

  if(R.version$major < 4 | R.version$minor < 2.0){ # Check for right version of R
    ancombc_out <- ANCOMBC::ancombc(
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
    ) %>% suppressMessages() %>% suppressWarnings()
    ancom_taxon <- row.names(ancombc_out$res$q_val)
    ancom_q <- ancombc_out$res$q_val[,1]
  } else {
    ancombc_out <- ANCOMBC::ancombc(
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
    ) %>% suppressMessages() %>% suppressWarnings()
    ancom_taxon <- ancombc_out$res$q_val[,1]
    ancom_q <- c(ancombc_out$res$q_val[,3])
  }

  message("Now running ALDEx2") # ALDEx2

  aldex2_out <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps_filt)),
                             phyloseq::sample_data(ps_filt)[[group]],
                             test = "t", effect = TRUE,
                             denom = "iqlr") %>% suppressMessages() %>% suppressWarnings()


  # Creating final data.frame of results

  message("Summarizing results")

  # Extract Benjamini-Hochberg adjusted p-values
  summ <- dplyr::full_join(data.frame(taxon = row.names(aldex2_out),
                                      aldex2 = (aldex2_out$wi.eBH)),
                           data.frame(taxon = ancom_taxon,
                                      ancombc = ancom_q),
                           by = "taxon") %>%
    dplyr::full_join(data.frame(taxon = row.names(deseq2_out),
                                deseq2 = deseq2_out$padj),
                     by="taxon")

  # Apply fisher method for meta-analysis of p-values
  meta_fisher_out <- apply(summ[,2:4], 1, metap::sumlog)
  meta_fish_p <- unlist(do.call(rbind, meta_fisher_out)[,3])
  summ <- cbind(summ, fisher_p = meta_fish_p)

  # Provide total number of DAA methods which found significance
  summ <- cbind(summ, rawcount = ((summ$aldex2 < 0.05)
                                             + (summ$ancombc < 0.05)
                                             + (summ$deseq2 < 0.05)))

  return(summ)

}



# [END]
