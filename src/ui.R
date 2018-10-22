#########################################################################################################################################
### Design the dashboad, Make the header
#########################################################################################################################################
header <- dashboardHeader(
  title = "DEkupl viewer"
  # dropdownMenu(
  #   type = "messages",
  #   messageItem(
  #     from = "Support",
  #     message = "The new server is ready.",
  #     icon = icon("life-ring"),
  #   )
  # )
)

#########################################################################################################################################
### Make sidebar Menu
#########################################################################################################################################
sidebar <- dashboardSidebar(
  disable = TRUE,            # It's disabled
  sidebarMenu(
    menuItem("DEkupl Annotation", tabName = "annot", icon = icon("file-archive-o"))
  )
)

#########################################################################################################################################
### Make body
#########################################################################################################################################
body <- dashboardBody(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  tabItems(
    tabItem(tabName = "annot",
      ######### Filters columns
      column(width = 3,
        fluidRow(
          selectInput("preset", "Filter Presets", choices=as.character(list('-'))),
          hr()
        ),
        fluidRow(
          div(style="margin-bottom:10px;", actionButton("reset", "Reset All", icon = icon("glyphicon glyphicon-repeat", lib = "glyphicon")))
        ),
        fluidRow(
          box(
            width = 12, title = "Binary options", status = "primary", solidHeader = TRUE, collapsible = TRUE,
            radioButtons("isMapped",    label = "Is mapped",                        choices = list("yes" = "TRUE", "no" = "FALSE", "both" = "NA"),  selected = "NA", inline = T),
            radioButtons("geneIsDiff",  label = "Gene is differentially expressed", choices = list("yes" = "TRUE", "no" = "FALSE", "both" = "NA"),  selected = "NA", inline = T),
            radioButtons("hasGene",     label = "Affected to a gene",               choices = list("yes" = "TRUE", "no" = "FALSE", "both" = "NA"),  selected = "NA", inline = T),
            radioButtons("hasASGene",   label = "Affected to an anti-sens gene",    choices = list("yes" = "TRUE", "no" = "FALSE", "both" = "NA"),  selected = "NA", inline = T),
            radioButtons("isExonic",    label = "Is exonic",                        choices = list("yes" = "TRUE", "no" = "FALSE", "both" = "NA"),  selected = "NA", inline = T),
            radioButtons("isIntronic",  label = "Is intronic",                      choices = list("yes" = "TRUE", "no" = "FALSE", "both" = "NA"),  selected = "NA", inline = T)
          )
        ),
        ## pvalue filter
        fluidRow(
          box(
            width = 12, title = "Pvalue", status = "primary", solidHeader = TRUE, collapsible = TRUE, "Pvalue distribution", br(),
            withSpinner(plotOutput('pvaluePlot', height = "50px")),
            sliderInput(
              inputId = "pvalue",
              label = "values",
              min = 0,
              max = 1,
              step = 0.00001,
              value = c(0,0.1)
            ),
            div(style="display:inline-block; width: 49%", numericInput("minPvalue", "min", value= 0, step = 0.001)),
            div(style="display:inline-block; width: 49%", numericInput("maxPvalue", "max", value= 0.1, step = 0.001)),
            checkboxInput("pvalueKeepNA", "Keep NA", value = TRUE)
          )
        ),
        ## du_pvalue filter
        fluidRow(
          box(
            width = 12, title = "DU Pvalue", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "Du Pvalue distribution", br(),
            withSpinner(plotOutput('duPvaluePlot', height = "50px")),
            sliderInput(
              inputId = "duPvalue",
              label = "values",
              min = 0,
              max = 1,
              step = 0.00001,
              value = c(0,0.1)
            ),
            div(style="display:inline-block; width: 49%", numericInput("minDuPvalue", "min", value= 0, step = 0.001)),
            div(style="display:inline-block; width: 49%", numericInput("maxDuPvalue", "max", value= 0.1, step = 0.001)),
            checkboxInput("duPvalueKeepNA", "Keep NA", value = TRUE)
          )
        ),
        ## clipped_3p filter
        fluidRow(
          box(
            width = 12, title = "3p clipped number", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "Clipped 3p distribution", br(),
            withSpinner(plotOutput('clipped3pPlot', height = "50px")),
            sliderInput(
              inputId = "clipped3p",
              label = "values",
              min = 0,
              max = 200,
              value = c(0, 200)
            ),
            div(style="display:inline-block; width: 49%", numericInput("minClipped3p", "min", value= 0)),
            div(style="display:inline-block; width: 49%", numericInput("maxClipped3p", "max", value= 200)),
            checkboxInput("clipped3pKeepNA", "Keep NA", value = TRUE)
          )
        ),
        ## nb_splice filter
        fluidRow(
          box(
            width = 12, title = "Number of splices", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "Number of Splices", br(),
            withSpinner(plotOutput('nbSplicePlot', height = "50px")),
            sliderInput(
              inputId = "nbSplice",
              label = "values",
              min = 0,
              max = 10,
              value = c(0, 10)
            ),
            div(style="display:inline-block; width: 49%", numericInput("minNbSplice", "min", value= 0)),
            div(style="display:inline-block; width: 49%", numericInput("maxNbSplice", "max", value= 10)),
            checkboxInput("nbSpliceKeepNA", "Keep NA", value = TRUE)
          )
        ),
        ## nb_snv filter
        fluidRow(
          box(
            width = 12, title = "Number of SNV", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "Number of SNV", br(),
            withSpinner(plotOutput('nbSnvPlot', height = "50px")),
            sliderInput(
              inputId = "nbSnv",
              label = "values",
              min = 0,
              max = 50,
              value = c(0, 50)
            ),
            div(style="display:inline-block; width: 49%", numericInput("minNbSnv", "min", value= 0)),
            div(style="display:inline-block; width: 49%", numericInput("maxNbSnv", "max", value= 50)),
            checkboxInput("nbSnvKeepNA", "Keep NA", value = TRUE)
          )
        ),
        ## nb_hit filter
        fluidRow(
          box(
            width = 12, title = "Number of Hits", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "Number of Hits", br(),
            withSpinner(plotOutput('nbHitPlot', height = "50px")),
            sliderInput(
              inputId = "nbHit",
              label = "values",
              min = 0,
              max = 300,
              value = c(0, 300)
            ),
            div(style="display:inline-block; width: 49%", numericInput("minNbHit", "min", value= 0)),
            div(style="display:inline-block; width: 49%", numericInput("maxNbHit", "max", value= 300)),
            checkboxInput("nbHitKeepNA", "Keep NA", value = TRUE)
          )
        ),
        ## contig length filter
        fluidRow(
          box(
            width = 12, title = "Contig length", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "Contig length", br(),
            withSpinner(plotOutput('contigSizePlot', height = "50px")),
            sliderInput(
              inputId = "contigSize",
              label = "values",
              min = 1,
              max = 1200,
              step = 50,
              value = c(1, 1200)
            ),
            div(style="display:inline-block; width: 49%", numericInput("minContigSize", "min", value= 1)),
            div(style="display:inline-block; width: 49%", numericInput("maxContigSize", "max", value= 1200)),
            checkboxInput("contigSizeKeepNA", "Keep NA", value = TRUE)
          )
        ),
        ## custom filter
        fluidRow(
          box(
            width = 12, title = "Customized filter", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "Free R code. exemple: grepl('(AATAAA|ATTAAA|AGTAAA|TATAAA).*AAAAA$', contig)", br(),
            textInput("customFilter", "Input text")
          )
        )
      ),
      # Table column
      column(width = 9,
        tabsetPanel(
          tabPanel(
            title = "Table",
            fluidRow(
              box(
                title = "contigs", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                strong(textOutput("datatableSelectedItems"), align = "center"),
                div(style="position:relative;margin-bottom:10px;", actionButton("showCols", "Choose cols", icon = icon("glyphicon glyphicon-list", lib = "glyphicon")), downloadButton(outputId = "downloadTable", label = "Download")),
                withSpinner(DT::dataTableOutput("table"))
              )
            )
          ),
          tabPanel(
            title = "Heatmap",
            fluidRow(
              box(
                title = "contigs", status = "primary", solidHeader = TRUE, collapsible = FALSE, width = 12, height = 800,
                strong(textOutput("heatmapSelectedItems"), align = "center"),
                # downloadButton(outputId = "downloadHeatmap", label = "Download"),
                withSpinner(plotOutput("heatmap"))
              )
            )
          ),
          tabPanel(
            title = "PCA",
            fluidRow(
              box(
                title = "contigs", status = "primary", solidHeader = TRUE, collapsible = FALSE, width = 12, height = 800,
                strong(textOutput("pcaSelectedItems"), align = "center"),
                div(style="position:relative;left:90%", downloadButton(outputId = "downloadPCA", label = "Download")),
                withSpinner(plotOutput("pca"))
              )
            )
          ),
          tabPanel(
            title = "Volcano plot",
            fluidRow(
              box(
                title = "contigs", status = "primary", solidHeader = TRUE, collapsible = FALSE, width = 12, height = 800,
                strong(textOutput("volcanoSelectedItems"), align = "center"),
                div(style="position:relative;left:90%", downloadButton(outputId = "downloadVolcano", label = "Download")),
                withSpinner(plotOutput("volcano"))
              )
            )
          )
        )
      )
    )
  )
)

#########################################################################################################################################
# define ui app
#########################################################################################################################################
ui <- dashboardPage(header, sidebar, body, skin="black")