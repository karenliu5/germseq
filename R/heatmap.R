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
                        values_to = "pval")

  heatmap_input <- heatmap_input %>%
    dplyr::mutate(text = paste0("Taxon: ", taxon, "\n",
                         "Method: ", method, "\n",
                         "p-value: ", -(pval^10), "\n",
                         "-log10(p): ", round(pval, 4)))

  pal <- wesanderson::wes_palette("Zissou1", 5, type = "continuous")

  p <- heatmap_input %>% ggplot2::ggplot(aes(x = taxon, y = method, fill = pval)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colours = c("#000000", pal),
                                  name = "-log10(p)",
                                  limit = c(-log10(0.5),
                                            min(30, max(heatmap_input$pval))),
                                  oob=scales::squish) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::ggtitle("-log10(p) for Each Method")

  p <- plotly::ggplotly(p, tooltip="text")

  return(p)
}
