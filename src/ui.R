#########################################################################################################################################
### Design the dashboad, Make the header
#########################################################################################################################################
header <- dashboardHeader(
  title = "DEkupl Annot dashboard"
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
        # fluidRow(
        #   selectInput("preset", "Filter Presets", choices=as.character(list('-', 'splice', 'splice_DU', 'polyA', 'polyA_DU', 'antisense'))),  # choices=colnames(df)
        #   hr()
        # ),
        # fluidRow(
        #   actionButton("reset", "Get All"),
        #   helpText("Get all items, including NA values")
        # ),
        fluidRow(
          # helpText("Allows to Get all items, including NA values"),
          checkboxInput("switchToFilterMode", "Disable/Activate filters", value = TRUE)
        ),
        ## du_pvalue filter
        fluidRow(
          box(
            width = 12, title = "Pvalue", status = "primary", solidHeader = TRUE, collapsible = TRUE, "some description", br(),
            sliderInput(
              inputId = "pvalue",
              # label = textOutput("max_pvalue"),
              label = "value",
              min = 0,
              max = 1,
              step = 0.00001,
              value = c(0,0.1)
            )
          )
        ),
        ## clipped_3p filter
        fluidRow(
          box(
            width = 12, title = "3p clipped number", status = "primary", solidHeader = TRUE, collapsible = TRUE, "some description", br(),
            sliderInput(
              inputId = "clipped3p",
              label = textOutput("max_clipped3p"),
              min = 0,
              max = 200,
              value = c(1, 200)
            )
          )
        ),
        ## nb_splice filter
        fluidRow(
          box(
            width = 12, title = "Number of splices", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "some description", br(),
            sliderInput(
              inputId = "nbSplice",
              label = textOutput("max_splice"),
              min = 0,
              max = 10,
              value = c(0, 10)
            )
          )
        ),
        ## nb_snv filter
        fluidRow(
          box(
            width = 12, title = "Number of SNV", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "some description", br(),
            sliderInput(
              inputId = "nbSnv",
              label = textOutput("max_snv"),
              min = 0,
              max = 50,
              value = c(0, 50)
            )
          )
        ),
        ## nb_hit filter
        fluidRow(
          box(
            width = 12, title = "Number of Hits", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "some description", br(),
            sliderInput(
              inputId = "nbHit",
              label = textOutput("max_hit"),
              min = 0,
              max = 300,
              value = c(0, 300)
            )
          )
        ),
        ## contig length filter
        fluidRow(
          box(
            width = 12, title = "Contig length", status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "some description", br(),
            sliderInput(
              inputId = "contigSize",
              label = textOutput("max_contig"),
              min = 1,
              max = 1200,
              step = 50,
              value = c(30, 1200)
            )
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
                title = "contigs", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12, strong(textOutput("datatableSelectedItems"), align = "center"), DT::dataTableOutput("table")
              )
            )
          ),
          tabPanel(
            title = "Heatmap",
            fluidRow(
              box(
                title = "contigs", status = "primary", solidHeader = TRUE, collapsible = FALSE, width = 12, height = 800, strong(textOutput("heatmapSelectedItems"), align = "center"), plotOutput("heatmap")
              )
            )
          ),
          tabPanel(
            title = "PCA",
            fluidRow(
              box(
                title = "contigs", status = "primary", solidHeader = TRUE, collapsible = FALSE, width = 12, height = 800, strong(textOutput("pcaSelectedItems"), align = "center"), plotOutput("pca")
              )
            )
          ),
          tabPanel(
            title = "Volcano plot",
            fluidRow(
              box(
                title = "contigs", status = "primary", solidHeader = TRUE, collapsible = FALSE, width = 12, height = 800, strong(textOutput("volcanoSelectedItems"), align = "center"), plotOutput("volcano")
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