library(shiny)

ui <- fluidPage(

  # Title
  titlePanel("Meta-Analysis of DAA using Fisher's Method"),

  # Add sidebar layout containing input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      br(),

      # Example Datasets ----
      tags$b("Example Datasets"),
      tags$p("Adapted from the atlas1006 dataset from the R package, microbiome."),
      tags$p("Lahti, L., Salojarvi, J., Salonen, A., Scheffer, M., and W. de Vos. (2014). Tipping elements in the human intestinal ecosystem. Nature Comm 5(4344). doi: https://doi.org/10.1038/ncomms5344."),
      uiOutput("tab1"),
      uiOutput("tab2"),
      uiOutput("tab3"),

      br(),

      # Inputs ----
      fileInput(inputId = "otu",
                label = "OTU table",
                accept = c(".csv")),
      fileInput(inputId = "tax",
                label = "Tax Table",
                accept = c(".csv")),
      fileInput(inputId = "samp",
                label = "Sample Data",
                accept = c(".csv")),
      textInput(inputId = "group",
                label = "variable of interest"),
      numericInput(inputId = "prevThr",
                   label = "Prevalence Threshold",
                   value = 0.1,
                   min = 0,
                   max = 1,
                   step = 0.05),

      br(),

      actionButton(inputId = "button2",
                   label = "Run"),

      width = 3

    ),

    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Results",
                           br(),
                           DT::dataTableOutput("tab")),
                  tabPanel("Heat map",
                           br(),
                           plotly::plotlyOutput("heatmap")),
                  tabPanel("Raw plot for Pie",
                           br(),
                           plotly::plotlyOutput("pie")),
                  tabPanel("Bar chart",
                           br(),
                           plotly::plotlyOutput("bar"))
                  ),

      width=8
    )
  )
)

server <- function(input, output) {

  otuInput <- eventReactive(eventExpr = input$button2, {
    if (! is.null(input$otu))
      as.matrix(read.csv(input$otu$datapath,
                         sep = ",",
                         header = TRUE,
                         row.names = 1,
                         check.names = FALSE))
  })

  taxInput <- eventReactive(eventExpr = input$button2, {
    if (! is.null(input$tax))
      as.matrix(read.csv(input$tax$datapath,
                         sep = ",",
                         header = TRUE,
                         row.names = 1,
                         check.names = FALSE))
  })

  sampInput <- eventReactive(eventExpr = input$button2, {
    if (! is.null(input$samp))
      read.csv(input$samp$datapath,
                         sep = ",",
                         header = TRUE,
                         row.names = 1)
  })

  physeq <- eventReactive(eventExpr = input$button2, {
    phyloseq::phyloseq(
      phyloseq::otu_table(otuInput(), taxa_are_rows = TRUE),
      phyloseq::tax_table(taxInput()),
      phyloseq::sample_data(sampInput()))
  })

  groupInput <- eventReactive(eventExpr = input$button2, {
    if(! is.null(input$group))
      input$group
  })

  thrInput <- eventReactive(eventExpr = input$button2, {
    if(! is.null(input$prevThr))
      input$prevThr
  })

  startDAA <- eventReactive(eventExpr = input$button2, {
    germseq::compare_DAA_methods(
      ps = physeq(),
      group = groupInput(),
      prevThr = thrInput())
  })

  output$tab <- DT::renderDataTable({
    if(! is.null(startDAA))
      startDAA()
  })

  output$heatmap <- plotly::renderPlotly({
    if(! is.null(startDAA))
      germseq::plot_result_heatmap(daa_output = startDAA()) %>% plotly::layout(height = 600)
  })

  output$pie <- plotly::renderPlotly({
    if(! is.null(startDAA))
      germseq::visualize_overlap(daa_output = startDAA())
  })

  output$bar <- plotly::renderPlotly({
    if(! is.null(startDAA))
      germseq::visualize_performances(daa_output = startDAA())
  })

  url1 <- a("Example OTU table", href="https://raw.githubusercontent.com/karenliu5/germseq/master/inst/extdata/atlas1006_otu.csv")
  output$tab1 <- renderUI({
    tagList("Download:", url1)
  })

  url2 <- a("Example tax table", href="https://raw.githubusercontent.com/karenliu5/germseq/master/inst/extdata/atlas1006_samp.csv")
  output$tab2 <- renderUI({
    tagList("Download:", url2)
  })

  url3 <- a("Example sample data", href="https://raw.githubusercontent.com/karenliu5/germseq/master/inst/extdata/atlas1006_tax.csv")
  output$tab3 <- renderUI({
    tagList("Download:", url3)
  })
}

shinyApp(ui, server)
