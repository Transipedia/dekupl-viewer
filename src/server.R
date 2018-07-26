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

  outputSelectedItems <- function() {
    d <- renderText({ 
      paste("Selected items: ", nrow(dataTableFilters()))
    })
    return (d)
  }

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
        conditions <- samplesCdtDF[,'condition']
        mat = as.matrix(data[, c(samples)])
        base_mean = rowMeans(mat)
        mat_scaled = t(apply(mat, 1, scale))
        formatSample = gsub("s\\d+_", "", colnames(mat))

        haConditions = HeatmapAnnotation(df = data.frame(conditions = conditions), show_annotation_name = TRUE)
        haSamples = HeatmapAnnotation(df = data.frame(samples = formatSample), show_annotation_name = TRUE)

        Heatmap(
          mat_scaled, name = "expression", km = 5, col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
          top_annotation = haConditions,
          bottom_annotation = haSamples,
          cluster_columns = TRUE,
          # column_title_gp = gpar(),
          show_row_names = FALSE,
          show_column_names = TRUE
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

  output$heatmapSelectedItems <- outputSelectedItems()

  ######################################  PCA page

  output$pca <- renderPlot(
    {
      data = dataTableFilters()
      if (nrow(data) > 0) {
        samples <- do.call(paste, c(as.list(samplesCdtDF['sample']), sep = ""))
        mat = t(as.matrix(data[, c(samples)]))

        res.pca <- prcomp(mat, scale = TRUE)
        # Graphique des individus.
        fviz_pca_ind(
          res.pca,
          col.ind = "cos2", # Colore by cos2 ((qualité de représentation))
          gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
          repel = TRUE
        )
        # Graphique des variables. Coloration en fonction de la contribution des variables. (Mais pas de header sur la mat transposée)
        # fviz_pca_var(res.pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
        # Biplot des individus et des variables
        # fviz_pca_biplot(res.pca, repel = TRUE, col.var = "#2E9FDF", col.ind = "#696969")
      } else {
        print('no data')
      }
    },
    height=700, bg="transparent"
  )
  
  output$pcaSelectedItems <- outputSelectedItems()
  
  ######################################  Volcano page
  # Genes that are highly dysregulated are farther to the left and right sides, while highly significant changes appear higher on the plot.
  output$volcano <- renderPlot(
    {
      data = dataTableFilters()
      if (nrow(data) > 0) {
        with(data, plot(log2FC, -log10(pvalue), pch=20, xlim=c(-2.5, 2)))

        # Add colored points: red if du_pvalue (padj ?) <0.05, orange of log2FC>1, green if both)
        with(subset(data, du_pvalue<.05 ), points(log2FC, -log10(pvalue), pch=20, col="red"))
        with(subset(data, abs(log2FC)>1), points(log2FC, -log10(pvalue), pch=20, col="orange"))
        with(subset(data, du_pvalue<.05 & abs(log2FC)>1), points(log2FC, -log10(pvalue), pch=20, col="green"))

        # Label points with the textxy function from the calibrate plot
        library(calibrate)
        with(subset(data, du_pvalue<.05 & abs(log2FC)>1), textxy(log2FC, -log10(pvalue), labs=gene_symbol, cex=.8))
        legend(x="topleft",
          legend=c("du_pvalue < 0.05", "log2FC > 1", "both"),
          col=c("red", "orange", "green"), box.lty=0, pch=20
        )
      } else {
        print('no data')
      }
    },
    height=700, bg="transparent"
  )

  output$volcanoSelectedItems <- outputSelectedItems()



}
