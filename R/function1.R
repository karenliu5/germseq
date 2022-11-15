### Gotta add the info LMAO
### Analysis functionality
### First, the group shit -> but additionally, what is the effect of DIVERSITY
### Then, survival function -> what is the effect of diversity (alpha diversity)
### /diff. concentrations of survival

# Takes in a phyloseq object, only takes a binary argument
DifferencesFunction <- function(ps, group, prev.thr = 0.1){
  # Check prev.thr
  if(prev.thr >= 1 | prev.thr < 0){
    stop("Error: prev.thr must be a number between [0, 1).")
  }

  # Check for Phyloseq object

  # Check if group is valid
  if(! assertthat::is.string(group)){ # Group is a string
    stop("Error: group must be a string.")
  }

  sample_meta <- sample_data(ps)
  meta_names <- colnames(sample_meta)

  if(! (group %in% meta_names)){
    stop("Error: group could not be found in phyloseq object.")
  } else {
    ind <- which(group == meta_names)
  }

  #Automatically remove NAs
  num_NA <- sum(is.na(sample_data(ps)[,ind]))
  if(num_NA > 0){
    var_values <- sample_data(ps)[[group]]
    ps <- prune_samples(!is.na(var_values), ps)
    warning("Warning: Rows with missing values were removed")
    return(ps)
  }

  #return()
  # if(num_NA > 0){
  #   ps <- subset_samples(ps, stats::complete.cases(sample_data(ps)))
  #    warning("Warning: Rows with missing values in sample data were removed.")
  # }

  # Check for two groups
  if(nrow(unique(sample_data(ps)[,ind])) != 2){
    stop("Error: Make sure group takes only two values.")
  }

  # Preprocess Data
  prevdf <- microbiome::prevalence(ps)

  mask <- as.logical(prevdf > prev.thr)
  keepTaxa <- names(prevdf[mask])

  ps_filt <- prune_taxa(keepTaxa, ps)

  # DESeq2
  form <- as.formula(paste("~", group))
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

  return(aldex2_da)
  # diagdds <- phyloseq_to_deseq2(atlas1006_genus, ~ bmi_group)
  # diagdds <- DESeq2::DESeq(diagdds, test = "Wald", fit="parametric")
  #
  # res <- DESeq2::results(diagdds)
  # df <- as.data.frame(res)

}

# library(microbiome)
# data(atlas1006)
# sample_data(atlas1006)$bmi_group
# atlas1006 <- subset_samples(atlas1006, (bmi_group == "lean" | bmi_group == "obese"))
# atlas1006 <- subset_samples(atlas1006, stats::complete.cases(sample_data(atlas1006)))
# DifferencesFunction(atlas1006, "bmi_group")
#
# mask <- as.logical(microbiome::prevalence(atlas1006) > 0.1)
# sum(mask)
# names(microbiome::prevalence(atlas1006)[mask])
#
# atlas1006 <- subset_samples(atlas1006, ! is.na(bmi_group))
# group <- as.formula("~ bmi_group")
# sample_data(atlas1006)$bmi_group
# atlas1006 <- subset_samples(atlas1006, ! is.na(sample_data(atlas1006)[,7]))
#
#
# deseq2_format <- phyloseq::phyloseq_to_deseq2(atlas1006, group)
#
# ALDEx2::aldex(data.frame(phyloseq::otu_table(atlas1006)),
#               phyloseq::sample_data(atlas1006)[["bmi_group"]],
#               test = "t", effect = TRUE,
#               denom = "iqlr")
#
# phyloseq::sample_data(atlas1006)$bmi_group
#
# group <- "bmi_group"
# sample_data(atlas1006)[[group]]
#
# atlas1006
# nrow(otu_table(atlas1006))
# as.vector(as.list(phyloseq::sample_data(atlas1006)[,7])[1])
# phyloseq::sample_data(atlas1006)$bmi_group
# nrow(data.frame(phyloseq::otu_table(atlas1006)))
# # ind <- which("bmi_group" == colnames(sample_data(atlas1006)))
# #
# "bmi_grou" %in% colnames(sample_data(atlas1006))
#
# nrow(unique(sample_data(atlas1006)[,ind]))
#
# library(microbiome)
# data(atlas1006)
#
# typeof(atlas1006)
# length(unique(sample_data(atlas1006)$bmi_group))
#

# atlas1006 <- subset_samples(atlas1006, (bmi_group == "lean" | bmi_group == "obese"))
#
# prevdf <- apply(X = otu_table(atlas1006),
#                  MARGIN = ifelse(taxa_are_rows(atlas1006), yes = 1, no = 2),
#                  FUN = function(x){sum(x > 0)})
# prevdf <- data.frame(Prevalence = prevdf,
#                      TotalAbundance = taxa_sums(atlas1006),
#                      tax_table(atlas1006))
#
# prevThres <- 0.10 * nsamples(atlas1006)
# keepTaxa <- rownames(prevdf)[(prevdf$Prevalence >= prevThres)]
# ps2 <- prune_taxa(keepTaxa, atlas1006)
#
# length(get_taxa_unique(ps2, taxonomic.rank="Genus"))
# atlas1006_genus <- phyloseq::tax_glom(ps2, taxrank="Genus", NArm= TRUE)
# atlas1006
#
# diagdds <- phyloseq_to_deseq2(atlas1006_genus, ~ bmi_group)
# diagdds <- DESeq2::DESeq(diagdds, test = "Wald", fit="parametric")
#
# res <- DESeq2::results(diagdds)
# df <- as.data.frame(res)
#
# sample_data(atlas1006_genus)
# tax_table(atlas1006_genus)
#
# out = ANCOMBC::ancombc(
#   phyloseq = atlas1006_genus,
#   formula = "bmi_group",
#   p_adj_method = "fdr",
#   zero_cut = 0.9,
#   lib_cut = 0,
#   group = "bmi_group",
#   struc_zero = TRUE,
#   neg_lb = TRUE,
#   tol = 1e-5,
#   max_iter = 100,
#   conserve = TRUE,
#   alpha = 0.05,
#   global = TRUE
# )
#
# res2 <- out$res
#
# res2$diff_abn



