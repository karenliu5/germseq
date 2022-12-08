#' Visualize performances of different methods of microbial differential abundance analysis (DAA)
#'
#' `visualize_performances()` visualizes the proportion of significant taxon
#' identified by three popular current methods of microbial differential abundance analysis,
#' DESeq2, ANCOM-BC, and ALDEx2.
#'
#' @param daa_output a data frame containing the results of the three methods, for format refer
#' to output from germseq::compare_DAA_methods() and germseq::atlas1006_output
#'
#' @returns A stacked barchart containing the proportion of Taxa found significantly different
#' by the methods.
#'
#' @examples
#' visualize_performances(germseq::atlas1006_output)
#'
#' @references
#' Wickham, H. (2022). ggplot2: Create Elegant Data
#'Visualisations Using the Grammar of Graphics. R package version 3.4.0.
#'\href{https://cran.r-project.org/web/packages/ggplot2/index.html}{Link}.
#'
#' @export
visualize_performances <- function(daa_output) {

  graph_input <- get_graph_input1(daa_output)

  # Create stacked barchart
  pl <- ggplot2::ggplot(graph_input, aes(fill=significance, y=value, x=methods)) +
    ggplot2::geom_bar(position="fill", stat="identity") +
    ggplot2::scale_fill_manual(values = wesanderson::wes_palette("Royal1")) +
    ggplot2::ggtitle("Proportion of Taxa found Significantly Different \nbetween Conditions by Method") +
    ggplot2::xlab("Methods") +
    ggplot2::ylab("Proportion")

  return(pl)
}

#'
#'Private Helper Function
#'
#' Process data into a suitable format for visualization with
#' visualize_performances()
#'
#' @param daa_ouput A dataframe to visualize.
get_graph_input1 <- function(daa_output){
  N <- nrow(daa_output)

  # Modify data.frame for suitable input into ggplot
  methods <- c(rep("aldex2", 2), rep("ancombc", 2), rep("deseq2", 2))
  significance <- rep(c("Significant", "Not Significant"), 3)
  value <- c(sum(daa_output$aldex2 < 0.05), N - sum(daa_output$aldex2 < 0.05),
             sum(daa_output$ancombc < 0.05), N - sum(daa_output$ancombc < 0.05),
             sum(daa_output$deseq2 < 0.05), N- sum(daa_output$deseq2 < 0.05))

  graph_input <- data.frame(methods, significance, value)

  return(graph_input)
}





#' Visualize overlap between differential abundance analysis (DAA) method results
#'
#' `visualize_overlap()` visualizes the overlap between taxon found significant
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
#' visualize_overlap(germseq::atlas1006_output)
#'
#' @export
visualize_overlap <- function(daa_output) {

  # Modify data.frame for suitable input into ggplot
  graph_input <- get_graph_input2(daa_output)

  # Create pie chart
  pie <- ggplot2::ggplot(graph_input, aes(fill=method, y=value, x="")) +
    ggplot2::geom_bar(width = 1, stat="identity") +
    ggplot2::scale_fill_manual(values = wesanderson::wes_palette("Darjeeling1")) +
    ggplot2::ggtitle("Proportion of Taxon found Significant \nby Number of Methods") +
    ggplot2::ylab("Proportion") +
    ggplot2::coord_polar("y", start = 0)

  return(pie)
}

#'
#'Private Helper Function
#'
#' Process data into a suitable format for visualization with
#' visualize_overlap
#'
#' @param daa_ouput A dataframe to visualize.
get_graph_input2 <- function(daa_output){
  score_summ <- table(daa_output$rawcount)

  graph_input <- data.frame(method =  c("no methods", "one method",
                         "two methods", "all methods"),
             value = c(as.numeric(score_summ)))

  return(graph_input)
}


# [END]
