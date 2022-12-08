#' Heatmap -log10(p)
#' @export
#' @importFrom magrittr "%>%"
plot_result_heatmap <- function(daa_output){
  daa_output <- daa_output %>%
    dplyr::mutate(
      deseq2 = -log10(deseq2),
      ancombc = -log10(ancombc),
      aldex2 = -log10(aldex2),
      fisher_p = -log10(fisher_p)
    )

  heatmap_input <- daa_output %>%
    tidyr::pivot_longer(c(deseq2, ancombc, aldex2, fisher_p),
                        names_to = "method",
                        values_to = "-log10(p)")

  pal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")
}
