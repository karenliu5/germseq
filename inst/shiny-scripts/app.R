
library(shiny)

ui <- fluidPage(
  titlePanel("Meta-Analysis of DAA using Fisher's Method"),

  sidebarLayout(

    sidebarPanel(
      fileInput(inputId = "otu",
                label = "OTU Table",
                accept = c(".csv")),
      fileInput(inputId = "tax",
                label = "Tax Table",
                accept = c(".csv")),
      fileInput(inputId = "samp",
                label = "Sample Table",
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
                   label = "Run")

    ),

    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Heat map",
                           plotOutput("heatmap")),
                  tabPanel("Raw plot for Pie",
                           plotOutput("pie")),
                  tabPanel("Bar chart",
                           plotOutput("bar"))
                  )
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
      phyloseq::otu_table(otuInput),
      phyloseq::tax_table(taxInput),
      phyloseq::sample_data(sampInput))
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
      ps = physeq,
      group = groupInput,
      prevThr = thrInput)
  })
}

shinyApp(ui, server)
