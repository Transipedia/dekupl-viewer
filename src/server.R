# Define server

server <- function(input, output, session) {
  set.seed(122)
  
  #########################################################################################################################################
  ###  load Data from input files
  #########################################################################################################################################

  loadContigsfile <- function() {
    d <- read.table(opt$contig, header=TRUE)
    print ('contigs file dim :'); print(dim(d))
    return (d) 
  }

  loadSamplesCdtfile <- function() {
    d <- read.table(opt$sample, header=TRUE)
    print ('sample conditions file dim :'); print(dim(d))
    return (d) 
  }

  contigsDF <- loadContigsfile()
  samplesCdtDF <- loadSamplesCdtfile()

  #########################################################################################################################################
  ###  Usual functions
  #########################################################################################################################################
  # TODO: use max value in slider UI
  maxOf <- function(col) {
    return (max(data[col], na.rm = TRUE))
  }

  # return selected line in the table
  selectedRow <- function() {
    return(dataTableFilters()[input$table_rows_selected,])
  }

  # return selected line infos
  selectedRowInfos <- function() {
    # ellipsis for contig
    # contig = ifelse(nchar(as.character(selectedRow()$contig)) <= 75, as.character(selectedRow()$contig), paste(substr(selectedRow()$contig, 1, 75), '...', sep=''))
    return (paste("Gene: ", selectedRow()$gene_symbol, " | Interval: ", selectedRow()$chromosome, ":", selectedRow()$start, "-", selectedRow()$end, " | CIGAR: ", selectedRow()$cigar, " | Strand: ", selectedRow()$gene_strand, " | biotype: ", selectedRow()$gene_biotype ))
  }

  # create new dataframe from samplesCdts and add sample differential expressions (from selectedRow)
  completeSampleCdts <- function() {
    newCdt = cbind(samplesCdtDF, data.frame(diff = 0))
    for(i in 1:nrow(newCdt)) {
      newCdt[i,'diff'] = selectedRow()[as.character(newCdt[i, 'sample'])]
    }
    return (newCdt)
  }

  ### filter dataTable depending input filters
  dataTableFilters <- reactive({
    s <- subset(
      contigsDF,
      du_pvalue <= input$duPvalue &
      nb_splice >= input$nbSplice[[1]] &
      nb_splice <= input$nbSplice[[2]] &
      clipped_3p >= input$clipped3p[[1]] &
      clipped_3p <= input$clipped3p[[2]] &
      nb_snv >= input$nbSnv[[1]] &
      nb_snv <= input$nbSnv[[2]] &
      nb_hit >= input$nbHit[[1]] &
      nb_hit <= input$nbHit[[2]] &
      contig_size >= input$contigSize[[1]] &
      contig_size <= input$contigSize[[2]]
    )
    # get cols we want to display
    samples <- do.call(paste, c(as.list(samplesCdtDF['sample']), sep = "")) # keep samples to be able to compute diff in boxplot
    s <- s[,c(samples, "contig", "tag", "contig_size", "chromosome", "start", "end", "gene_id", "gene_symbol", "gene_strand", "gene_biotype", "exonic", "intronic", "gene_is_diff", "cigar", "is_mapped", "pvalue", "du_pvalue", "du_stat", "meanA", "meanB", "log2FC", "nb_insertion", "nb_deletion", "nb_splice", "nb_snv", "clipped_3p", "clipped_5p", "is_clipped_3p", "is_clipped_5p", "query_cover", "alignment_identity", "nb_hit", "nb_mismatch", "strand", "as_gene_id", "as_gene_symbol", "as_gene_strand", "as_gene_biotype", "upstream_gene_id", "upstream_gene_strand", "upstream_gene_symbol", "upstream_gene_dist", "downstream_gene_id", "downstream_gene_strand", "downstream_gene_symbol", "downstream_gene_dist")]
    return(s)
  })

  #########################################################################################################################################
  ###  Outputs
  #########################################################################################################################################
  
  ######################################  Table page 

  # Set an event observable (on click) on table row, open a modal
  observeEvent(input$table_rows_selected, {
    showModal(
      modalDialog(
        title = selectedRowInfos(),
        size = 'l',
        helpText("Contig: ", selectedRow()$contig),  # TODO: this must be multi lines
        helpText("Tag: ", selectedRow()$tag),
        tabsetPanel(
          tabPanel(
            title = "Contig vizualiser",
            fluidRow(
              plotOutput(outputId="boxplotInModal")
            ),
            downloadButton(outputId = "downloadboxplot", label = "Download")
          )
        )
      )
    )
  })

  ## Draw boxplot in modal when click in a table row
  output$boxplotInModal <- renderPlot({
    boxplot()
  })

  boxplot <- function() {
    ggplot(data = completeSampleCdts(), aes(x = condition, y = diff, fill = condition)) +
      geom_boxplot() +
      # geom_point()
      geom_jitter(position = position_jitter(width = .2)) # disperse les points sur la largeur
  }

  ## download boxplot in modal
  output$downloadboxplot <- downloadHandler(
    filename = function() {
      paste(as.character(selectedRow()$gene_symbol), "_", as.character(selectedRow()$chromosome), ":", as.character(selectedRow()$start), "-", as.character(selectedRow()$end), '.png', sep='')
    },
    content = function(file) {
      ggsave(file, plot = boxplot(), device = "png", width=10, height=6)
    }
  )

  ###  Build the main table dataframe
  output$table <- DT::renderDataTable(
    {
      DT::datatable(
        dataTableFilters(),
        rownames = FALSE,
        escape = FALSE,
        selection = 'single',
        options = list(
          pageLength = 18,
          columnDefs = list(list(
            visible=FALSE,
            targets=c(0 : (nrow(samplesCdtDF['sample']) + 1))   # to hide samples and contig + tag cols
          ))
        )
      ) %>%
        formatRound(c("pvalue", "meanA", "meanB", "log2FC", "query_cover", "alignment_identity"), 3) %>%
        formatRound(c("du_pvalue", "du_stat"), 10) %>%
        formatStyle(c('chromosome', 'gene_symbol', 'gene_biotype'), fontWeight = 'bold') %>% 
        # formatStyle('is_mapped', fontWeight = 'bold', color = styleEqual(c('true', 'false'), c('green', 'red')))  %>%
        formatStyle('pvalue', target = 'row', cursor = 'pointer')
    },
    server = FALSE
  )

  ######################################  Heatmap page 

  output$heatmap <- renderPlot(
    {
      data = dataTableFilters()
      if (nrow(data) > 0) {
        samples <- do.call(paste, c(as.list(samplesCdtDF['sample']), sep = ""))
        mat = as.matrix(data[, c(samples)])
        base_mean = rowMeans(mat)
        mat_scaled = t(apply(mat, 1, scale))

        formatSample = gsub("s\\d+_", "", colnames(mat))
        ha = HeatmapAnnotation(df = data.frame(samples = formatSample))

        Heatmap(
          mat_scaled, name = "expression", km = 5, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
          top_annotation = ha, top_annotation_height = unit(4, "mm"), 
          show_row_names = FALSE, show_column_names = FALSE
        ) +
        Heatmap(base_mean, name = "base mean", show_row_names = FALSE, width = unit(5, "mm")) +
        Heatmap(
          data$contig_size, name = "contig size", col = colorRamp2(c(0, 1300), c("white", "orange")),
          heatmap_legend_param = list(
            at = c(0, 200, 500, 800, 1000, 1300), 
            labels = c("0b", "200b", "500b", "800b", "1Kb", "1.3kb")
          ),
          width = unit(5, "mm")
        ) +
        Heatmap(data$gene_biotype, name = "biotype", width = unit(5, "mm"))
      } else {
        print('no data')
      }
    },
    height=700, bg="transparent"
  )

  output$selectedItems <- renderText({ 
    paste("Selected items: ", nrow(dataTableFilters()))
  })

}
