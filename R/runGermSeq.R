#' Launch Shiny app for germseq
#'
#' A function that launches the Shiny app for germseq.
#' The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return Launches a Shiny page, no return value.
#'
#' @examples
#' \dontrun{
#' germseq::runTestingPackage()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' Silva, A. (2022) TestingPackage: An Example R Package For
#' BCB410H. Unpublished. https://github.com/anjalisilva/TestingPackage
#'
#' @export
#' @importFrom shiny runApp

runGermSeq <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "germseq")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}

#[END]
