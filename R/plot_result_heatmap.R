#' Visualize resulting p-values from differential abundance analysis (DAA) methods and meta-analysis using Fisher's method
#'
#' `plot_result_heatmap()` visualizes the -log(p-values) resulting from three
#' different methods of metagenomic DAA (DESeq2, ANCOM-BC, and ALDEx2), and their meta-analysis
#' using Fisher's method.
#'
#'  @param daa_output a data frame containing the results of the three methods, for format refer
#' to output from germseq::compare_DAA_methods() and germseq::atlas1006_output
#'
#' @returns An interactive heatmap containing the proportion of Taxa found significantly different
#' by the methods.
#'
#' @examples
#' plot_result_heatmap(germseq::atlas1006_output)
#'
#' @references
#' Dewey, M. (2022). metap: Meta-Analysis of Significance Values. R package version 1.8.
#' \href{https://cran.r-project.org/web/packages/metap/index.html}{Link}.
#'
#' Wickham, H. (2022). ggplot2: Create Elegant Data
#' Visualisations Using the Grammar of Graphics. R package version 3.4.0.
#' \href{https://cran.r-project.org/web/packages/ggplot2/index.html}{Link}.
#'
#' @export
#' @importFrom magrittr "%>%"
#'
plot_result_heatmap <- function(daa_output){

  # Mutate output to get -log p-values
  daa_output <- daa_output %>%
    dplyr::mutate(
      deseq2 = -log10(deseq2),
      ancombc = -log10(ancombc),
      aldex2 = -log10(aldex2),
      fisher_p = -log10(fisher_p)
    )

  # Change format of dataframe for input into graph
  heatmap_input <- daa_output %>%
    tidyr::pivot_longer(c(deseq2, ancombc, aldex2, fisher_p),
                        names_to = "method",
                        values_to = "pval")

  # Add additional text for interactive plot
  heatmap_input <- heatmap_input %>%
    dplyr::mutate(text = paste0("Taxon: ", taxon, "\n",
                         "Method: ", method, "\n",
                         "p-value: ", -(pval^10), "\n",
                         "-log10(p): ", round(pval, 4)))

  # Load palette
  pal <- wesanderson::wes_palette("Zissou1", 50, type = "continuous")

  # Create heatmap
  p <- heatmap_input %>%
    ggplot2::ggplot(mapping = aes(x = taxon, y = method, fill = pval, text=text)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = c("black", pal),
                                  name = "-log10(p)",
                                  limit = c(-log10(0.5), min(35, max(heatmap_input$pval))),
                                  oob=scales::squish) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::ggtitle("Resulting -log10(p) by Method")

  # Make heatmap interactive
  p <- plotly::ggplotly(p, tooltip="text")

  return(p)
}

# [END]
