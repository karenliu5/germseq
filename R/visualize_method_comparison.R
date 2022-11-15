#' Visualize performances of different methods of microbial differential abundance analysis (DAA)
#'
#' `visualize_performances_barchart()` visualizes the proportion of significant taxon
#' identified by three popular current methods of microbial differential abundance analysis,
#' DESeq2, ANCOM-BC, and ALDEx2.
#'
#' @param daa_output a data frame containing the results of the three methods, for format refer
#' to output from germseq::compare_DAA_methods()
#'
#' @returns A stacked barchart containing the proportion of Taxa found significantly different
#' by the methods.
#'
#' @examples
#' library(germseq)
#' data("atlas1006_output")
#' visualize_performances_barchart(atlas1006_output)
#'
#' @export
visualize_performances_barchart <- function(daa_output) {

  N <- nrow(daa_output)

  # Modify data.frame for suitable input into ggplot
  methods <- c(rep("aldex2", 2), rep("ancombc", 2), rep("deseq2", 2))
  significance <- rep(c("Significant", "Not Significant"), 3)
  value <- c(sum(daa_output$aldex2), N - sum(daa_output$aldex2),
             sum(daa_output$ancombc), N - sum(daa_output$ancombc),
             sum(daa_output$deseq2), N- sum(daa_output$deseq2))

  graph_input <- data.frame(methods, significance, value)


  # Create stacked barchart
  pl <- ggplot2::ggplot(graph_input, aes(fill=significance, y=value, x=methods)) +
    ggplot2::geom_bar(position="fill", stat="identity") +
    ggplot2::ggtitle("Proportion of Taxa found Significantly Different \nbetween Conditions by Method") +
    ggplot2::xlab("Methods") +
    ggplot2::ylab("Proportion")

  return(pl)
}


#' Visualize overlap between differential abundance analysis (DAA) method results
#'
#' `visualize_overlap_piechart()` visualizes the overlap between taxon found significant
#' identified by three popular current methods of microbial differential abundance analysis,
#' DESeq2, ANCOM-BC, and ALDEx2.
#'
#' @param daa_output a data frame containing the results of the three methods, for format refer
#' to output from germseq::compare_DAA_methods()
#'
#' @returns A pie chart containing the proportion of Taxa found significant by different
#' numbers of methods.
#'
#' @examples
#' library(germseq)
#' data("atlas1006_output")
#' visualize_overlap_piechart(atlas1006_output)
#'
#' @export
visualize_overlap_piechart <- function(daa_output) {

  # Modify data.frame for suitable input into ggplot
  score_summ <- table(daa_output$score)

  graph_input <- data.frame(method =  c("no methods", "one method",
                                        "two methods", "all methods"),
                            value = c(as.numeric(score_summ)))

  # Create pie chart
  pie <- ggplot2::ggplot(graph_input, aes(fill=method, y=value, x="")) +
    ggplot2::geom_bar(width = 1, stat="identity") +
    ggplot2::ggtitle("Proportion of Taxon found Significant \nby Number of Methods") +
    ggplot2::ylab("Proportion") +
    coord_polar("y", start = 0)

  return(pie)
}
