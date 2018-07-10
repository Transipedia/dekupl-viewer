#########################################################################################################################################
### Design the dashboad, Make the header
#########################################################################################################################################
header <- dashboardHeader(
  title = "DEkupl Annot dashboard",
  dropdownMenu(
    type = "messages",
    messageItem(
      from = "Support",
      message = "The new server is ready.",
      icon = icon("life-ring"),
    )
  )
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
        ## du_pvalue filter
        fluidRow(
          box(
            width = 12, title = "DU Pvalue", status = "info", solidHeader = TRUE, collapsible = TRUE, "some description", br(),
            sliderInput(
              inputId = "duPvalue",
              label="value",
              min = 0,
              max = 1,
              step = 0.00001,
              value = 0.01
            )
          )
        ),
        ## nb_splice filter
        fluidRow(
          box(
            width = 12, title = "Number of splices", status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "some description", br(),
            sliderInput(
              inputId = "nbSplice",
              label="value",
              min = 0,
              max = 10,
              value = c(0, 10)
            )
          )
        ),
        ## clipped_3p filter
        fluidRow(
          box(
            width = 12, title = "3p clipped number", status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "some description", br(),
            sliderInput(
              inputId = "clipped3p",
              label="value",
              min = 0,
              max = 100,
              value = c(1, 100)
            )
          )
        ),
        ## nb_snv filter
        fluidRow(
          box(
            width = 12, title = "Number of SNV", status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "some description", br(),
            sliderInput(
              inputId = "nbSnv",
              label="value",
              min = 0,
              max = 30,
              value = c(0, 30)
            )
          )
        ),
        ## nb_hit filter
        fluidRow(
          box(
            width = 12, title = "Number of Hits", status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "some description", br(),
            sliderInput(
              inputId = "nbHit",
              label="value",
              min = 0,
              max = 30,
              value = c(0, 30)
            )
          )
        ),
        ## contig length filter
        fluidRow(
          box(
            width = 12, title = "Contigs lenght", status = "info", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE, "some description", br(),
            sliderInput(
              inputId = "contigSize",
              label="value",
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
                title = "contigs", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 12, DT::dataTableOutput("table")
              )
            )
          ),
          tabPanel(
            title = "Heatmap",
            fluidRow(
              box(
                title = "contigs", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 12, uiOutput("heatmap")
              )
            )
          ),
          tabPanel(
            title = "PCA",
            fluidRow(
              box(
                title = "contigs", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 12, uiOutput("pca")
              )
            )
          ),
          tabPanel(
            title = "Volcano graph",
            fluidRow(
              box(
                title = "contigs", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 12, uiOutput("volcano")
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