# Define server

server <- function(input, output, session) {
  set.seed(122)
  
  #########################################################################################################################################
  ###  load Data from input files
  #########################################################################################################################################

  loadSamplesCdtfile <- function() {
    d <- read.table(opt$sample, header=TRUE, sep="\t")
    for (sample in d$sample) {
      if (grepl("(^[0-9])|-", sample)) {
        print ("ERROR: unexpected name of samples. Must not start with a number or contains a '-' character")
        stop()
      }
    }
    print ('sample conditions file dim :'); print(dim(d))
    return (d)
  }
  loadContigsfile <- function() {
    d <- read.table(opt$contig, header=TRUE, sep="\t")
    print ('contigs file dim :'); print(dim(d))
    return (d) 
  }

  loadFiltersPresetfile <- function() {
    d <- read.table(opt$filters, header=TRUE, sep="\t")
    print ('filters preset file dim :'); print(dim(d))
    return (d) 
  }

  contigsDF <- loadContigsfile()
  samplesCdtDF <- loadSamplesCdtfile()
  filtersPresetDF <- loadFiltersPresetfile()

  samples <- do.call(paste, c(as.list(samplesCdtDF['sample']), sep = ""))
  excludedCols <- c(samples, "contig", "tag", "nb_merged_kmers")
  currentCols = setdiff(names(contigsDF), excludedCols)

  minPvalues = min(contigsDF$pvalue, na.rm = TRUE);              maxPvalues = max(contigsDF$pvalue, na.rm = TRUE)
  minNbSplices = min(contigsDF$nb_splice, na.rm = TRUE);         maxNbSplices = max(contigsDF$nb_splice, na.rm = TRUE)
  minClipped3ps = min(contigsDF$clipped_3p, na.rm = TRUE);       maxClipped3ps = max(contigsDF$clipped_3p, na.rm = TRUE)
  minNbSnvs = min(contigsDF$nb_snv, na.rm = TRUE);               maxNbSnvs = max(contigsDF$nb_snv, na.rm = TRUE)
  minNbHits = min(contigsDF$nb_hit, na.rm = TRUE);               maxNbHits = max(contigsDF$nb_hit, na.rm = TRUE)
  minContigSizes = min(contigsDF$contig_size, na.rm = TRUE);     maxContigSizes = max(contigsDF$contig_size, na.rm = TRUE)
  minDuPvalues = min(contigsDF$du_pvalue, na.rm = TRUE);         maxDuPvalues = max(contigsDF$du_pvalue, na.rm = TRUE); 

  # update UI inputs
  updateSliderInput(session, "pvalue", value = NULL, min = minPvalues, max = maxPvalues)
  updateSliderInput(session, "nbSplice", value = NULL, min = minNbSplices, max = maxNbSplices)
  updateSliderInput(session, "clipped3p", value = NULL, min = minClipped3ps, max = maxClipped3ps)
  updateSliderInput(session, "nbSnv", value = NULL, min = minNbSnvs, max = maxNbSnvs)
  updateSliderInput(session, "nbHit", value = NULL, min = minNbHits, max = maxNbHits)
  updateSliderInput(session, "contigSize", value = NULL, min = minContigSizes, max = maxContigSizes)

  updateSelectInput(session, "preset", choices = c('-', as.vector(rbind(as.character(filtersPresetDF$event)))))

  #########################################################################################################################################
  ###  Usual functions
  #########################################################################################################################################

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
  
  getColsIndex <- function(cols) {
    l <- c()
    for (i in cols) {
      l <- c(l, which(colnames(contigsDF)==i) - 1)
    }
    return (l)
  }

  hidedColsIndexes <- function(cols) {
    if(!is.null(cols)) {
      return (getColsIndex(setdiff(names(contigsDF), cols)))
    }
    return (getColsIndex(excludedCols))
  }

  ######################################  filter dataTable depending input filters
  dataTableFilters <- reactive({
    # from classic filters
    inputPvalueMin <- input$pvalue[[1]];         inputPvalueMax <- input$pvalue[[2]]
    inputNbSpliceMin <- input$nbSplice[[1]];     inputNbSpliceMax <- input$nbSplice[[2]]
    inputClipped3pMin <- input$clipped3p[[1]];   inputClipped3pMax <- input$clipped3p[[2]]
    inputNbSnvMin <- input$nbSnv[[1]];           inputNbSnvMax <- input$nbSnv[[2]]
    inputNbHitMin <- input$nbHit[[1]];           inputNbHitMax <- input$nbHit[[2]]
    inputContigSizeMin <- input$contigSize[[1]]; inputContigSizeMax <- input$contigSize[[2]]
    customizedFilter <- if (input$customFilter != "") input$customFilter else TRUE
    # these values are defined in preset filters file. (transipedia.tsv)
    duPvalueMin <- minDuPvalues; duPvalueMax <- maxDuPvalues
    geneIsDiff <- c(FALSE, TRUE, NA)
    geneSymbol <- c(FALSE, TRUE)
    asGeneID <- c(FALSE, TRUE)
    isMapped <- c(FALSE, TRUE, NA)
    isExonic <- c(FALSE, TRUE, NA)
    isIntronic <- c(FALSE, TRUE, NA)
    customizedPresetFilter <- TRUE # this is R code defined preset-filters file, column other (transipedia.tsv)

    # from filters preset
    if (input$preset != '-') {
      filter <- filtersPresetDF[which(filtersPresetDF$event == input$preset),]

      inputPvalueMin <- minPvalues
      inputPvalueMax <- maxPvalues
      inputNbSpliceMin <- if (!is.na(filter$nb_splice_min)) filter$nb_splice_min else minNbSplices
      inputNbSpliceMax <- if (!is.na(filter$nb_splice_max)) filter$nb_splice_max else maxNbSplices
      inputClipped3pMin <- if (!is.na(filter$clipped_3p_min)) filter$clipped_3p_min else minClipped3ps
      inputClipped3pMax <- if (!is.na(filter$clipped_3p_max)) filter$clipped_3p_max else maxClipped3ps
      inputNbSnvMin <- if (!is.na(filter$nb_snv_min)) filter$nb_snv_min else minNbSnvs
      inputNbSnvMax <- if (!is.na(filter$nb_snv_max)) filter$nb_snv_max else maxNbSnvs
      inputNbHitMin <- if (!is.na(filter$nb_hit_min)) filter$nb_hit_min else minNbHits
      inputNbHitMax <- if (!is.na(filter$nb_hit_max)) filter$nb_hit_max else maxNbHits
      inputContigSizeMin <- if (!is.na(filter$contig_size_min)) filter$contig_size_min else minContigSizes
      inputContigSizeMax <- if (!is.na(filter$contig_size_max)) filter$contig_size_max else maxContigSizes

      duPvalueMin <- if (!is.na(filter$du_pvalue_min)) filter$du_pvalue_min else minDuPvalues
      duPvalueMax <- if (!is.na(filter$du_pvalue_max)) filter$du_pvalue_max else maxDuPvalues
      geneIsDiff <- if (!is.na(filter$gene_is_diff)) c(filter$gene_is_diff) else c(FALSE, TRUE, NA)
      geneSymbol <- if (is.na(filter$gene_symbol)) c(FALSE, TRUE) else if (filter$gene_symbol == 'unknown') c(FALSE) else c(TRUE)
      asGeneID <- if (is.na(filter$as_gene_id)) c(FALSE, TRUE) else if (filter$as_gene_id == 'unknown') c(FALSE) else c(TRUE)
      isMapped <- if (!is.na(filter$is_mapped)) c(filter$is_mapped) else c(FALSE, TRUE, NA)
      isExonic <- if (!is.na(filter$exonic)) c(filter$exonic) else c(FALSE, TRUE, NA)
      isIntronic <- if (!is.na(filter$intronic)) c(filter$intronic) else c(FALSE, TRUE, NA)
      customizedPresetFilter <- if (!is.na(filter$other)) c(filter$other) else TRUE

      cat(sprintf(
        "\n\ngene is diff: %s\ngene symbol: %s\nas gene id: %s\nis mapped: %s\nexonic: %s\nintronic: %s\nother: %s\n",
        filter$gene_is_diff, filter$gene_symbol, filter$as_gene_id, filter$is_mapped, filter$exonic, filter$intronic, filter$other
      ))

      # update UI slider inputs
      updateSliderInput(session, "pvalue", value = c(inputPvalueMin, inputPvalueMax))
      updateSliderInput(session, "nbSplice", value = c(inputNbSpliceMin, inputNbSpliceMax))
      updateSliderInput(session, "clipped3p", value = c(inputClipped3pMin, inputClipped3pMax))
      updateSliderInput(session, "nbSnv", value = c(inputNbSnvMin, inputNbSnvMax))
      updateSliderInput(session, "nbHit", value = c(inputNbHitMin, inputNbHitMax))
      updateSliderInput(session, "contigSize", value = c(inputContigSizeMin, inputContigSizeMax))
      updateTextInput(session, "customFilter", value = filter$other)
    }
    if (input$switchToFilterMode) {
      cat(sprintf(
        "\npvalue: [%f, %f]\nnb splice: [%f, %f]\nclipped 3p: [%f, %f]\nnb snv: [%f, %f]\nnb hit: [%f, %f]\ncontig size: [%f, %f]\ndu pvalue: [%f, %f]\ncustom filter: %s\n",
        inputPvalueMin, inputPvalueMax, inputNbSpliceMin, inputNbSpliceMax, inputClipped3pMin, inputClipped3pMax, inputNbSnvMin, inputNbSnvMax, inputNbHitMin, inputNbHitMax, inputContigSizeMin, inputContigSizeMax, duPvalueMin, duPvalueMax, input$customFilter
      ))
      
      # update UI numeric input
      updateNumericInput(session, "minPvalue", value = inputPvalueMin);         updateNumericInput(session, "maxPvalue", value = inputPvalueMax)
      updateNumericInput(session, "minNbSplice", value = inputNbSpliceMin);     updateNumericInput(session, "maxNbSplice", value = inputNbSpliceMax)
      updateNumericInput(session, "minClipped3p", value = inputClipped3pMin);   updateNumericInput(session, "maxClipped3p", value = inputClipped3pMax)
      updateNumericInput(session, "minNbSnv", value = inputNbSnvMin);           updateNumericInput(session, "maxNbSnv", value = inputNbSnvMax)
      updateNumericInput(session, "minNbHit", value = inputNbHitMin);           updateNumericInput(session, "maxNbHit", value = inputNbHitMax)
      updateNumericInput(session, "minContigSize", value = inputContigSizeMin); updateNumericInput(session, "maxContigSize", value = inputContigSizeMax)
      
      # get subset with filters
      s <- subset(
        contigsDF,
        pvalue >= inputPvalueMin & pvalue <= inputPvalueMax &
        nb_splice >= inputNbSpliceMin & nb_splice <= inputNbSpliceMax &
        clipped_3p >= inputClipped3pMin & clipped_3p <= inputClipped3pMax &
        nb_snv >= inputNbSnvMin & nb_snv <= inputNbSnvMax &
        nb_hit >= inputNbHitMin & nb_hit <= inputNbHitMax &
        contig_size >= inputContigSizeMin & contig_size <= inputContigSizeMax &
        du_pvalue >= duPvalueMin & du_pvalue <= duPvalueMax &
        is.element(gene_is_diff, geneIsDiff) &
        is.element(!is.na(gene_symbol), geneSymbol) &
        is.element(!is.na(as_gene_id), asGeneID) &
        is.element(is_mapped, isMapped) &
        is.element(exonic, isExonic) &
        is.element(intronic, isIntronic) &
        eval(parse(text = customizedFilter)) & 
        eval(parse(text = customizedPresetFilter))
      )
    } else {
      s <- contigsDF
    }
    # if(!is.null(input$select_cols)) {
    #   currentCols = input$select_cols
    # }
    return(s)
  })
  
  ####################################### observe only numeric inputs to update sliders 
  observe({
    updateSliderInput(session, "pvalue", value = c(input$minPvalue, input$maxPvalue))
    updateSliderInput(session, "clipped3p", value = c(input$minClipped3p, input$maxClipped3p))
    updateSliderInput(session, "nbSplice", value = c(input$minNbSplice, input$maxNbSplice))
    updateSliderInput(session, "nbSnv", value = c(input$minNbSnv, input$maxNbSnv))
    updateSliderInput(session, "nbHit", value = c(input$minNbHit, input$maxNbHit))
    updateSliderInput(session, "contigSize", value = c(input$minContigSize, input$maxContigSize))
  })

  #########################################################################################################################################
  ###  Main page Outputs
  #########################################################################################################################################

  ####################################### Choose cols
  observeEvent(input$showCols, {
    showModal(modalDialog(
      title = 'columns',
      size = 'm',
      checkboxGroupInput(inputId = "select_cols", 
        label = "Select cols", 
        choices = names(contigsDF),
        selected = currentCols
      )
    ))
  })

  ####################################### selected items
  outputSelectedItems <- function() {
    d <- renderText({
      paste("Selected items: ", nrow(dataTableFilters()), " / ", nrow(contigsDF))
    })
    return (d)
  }

  ######################################  On Table line selection 
  # Set an event observable (on click) on table row, open a modal
  observeEvent(input$table_rows_selected, {
    showModal(
      modalDialog(
        title = selectedRowInfos(),
        size = 'l',
        tabsetPanel(
          tabPanel(
            title = "Contig boxplot",
            fluidRow(
              plotOutput(outputId="boxplotInModal")
            ),
            downloadButton(outputId = "downloadboxplot", label = "Download")
          ),
          tabPanel(
            title = "Contig Details",
            fluidRow(
              DT::dataTableOutput(outputId="contigDetails")
            )
          )
        )
      )
    )
  })
  ## Contigs details
  output$contigDetails <- DT::renderDataTable({
     DT::datatable(
      t(selectedRow()),
      colnames = c('', ''),
      selection = 'single',
      options = list(
        searching = FALSE,
        paging = FALSE
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

  ##################################### Build the main table dataframe
  output$table <- DT::renderDataTable({
    DT::datatable(
      dataTableFilters(),
      rownames = FALSE,
      escape = FALSE,
      selection = 'single',
      options = list(
        pageLength = 15,
        columnDefs = list(list(
          visible=FALSE,
          targets=hidedColsIndexes(input$select_cols)
        ))
      )
    ) %>%
    formatRound(c("pvalue", "meanA", "meanB", "log2FC", "query_cover", "alignment_identity"), 3) %>%
    formatRound(c("du_pvalue", "du_stat"), 10) %>%
    formatStyle(c('chromosome', 'gene_symbol', 'gene_biotype'), fontWeight = 'bold') %>% 
    # formatStyle('is_mapped', fontWeight = 'bold', color = styleEqual(c('true', 'false'), c('green', 'red')))  %>%
    formatStyle('pvalue', target = 'row', cursor = 'pointer')
  }, server = TRUE )
  
  output$datatableSelectedItems <- outputSelectedItems()

  ######################################  Heatmap page 

  heatmap <- function() {
    data = dataTableFilters()
    if (nrow(data) > 0) {
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
  }
    
  output$heatmap <- renderPlot({
    heatmap()
  }, height=700, bg="transparent")

  ## download Heatmap
  output$downloadHeatmap <- downloadHandler(
    filename = "heatmap.png",
    content = function(file) {
      #preparing the dataset from an external function
      # corr_data <- dataTableFilters() 
      #colors of the heatmap
      col<- colorRampPalette(c("blue", "white", "red"))(20)
      png(file)
      # heatmap(x = cor(x =as.matrix(corr_data)) , col = col, symm = TRUE)
      heatmap()
      dev.off()
    }
  )

  output$heatmapSelectedItems <- outputSelectedItems()

  ######################################  PCA page

  pca <- function() {
    data = dataTableFilters()
    if (nrow(data) > 0) {
      
      conditions <- samplesCdtDF[,'condition']
      mat = t(as.matrix(data[, c(samples)]))

      res.pca <- prcomp(mat, scale = TRUE)
      # Graphique des individus.
      fviz_pca_ind(
        res.pca,
        col.ind = conditions,
        repel = TRUE
      )
    } else {
      print('no data')
    }
  }

  output$pca <- renderPlot({
    pca()
  }, height=700, bg="transparent")
  

  ## download pca
  output$downloadPCA <- downloadHandler(
    filename = "pca.png",
    content = function(file) {
      ggsave(file, plot = pca(), device = "png", width=10, height=6)
    }
  )

  output$pcaSelectedItems <- outputSelectedItems()
  
  ######################################  Volcano page
  # Genes that are highly dysregulated are farther to the left and right sides, while highly significant changes appear higher on the plot.
  volcano <- function() {
    data = dataTableFilters()
    if (nrow(data) > 0) {
      with(data, plot(log2FC, -log10(pvalue), pch=20, xlim=c(-2.5, 2)))

      # Add colored points: red if pvalue (padj ?) <0.05, orange of log2FC>1, green if both)
      with(subset(data, pvalue<.05 ), points(log2FC, -log10(pvalue), pch=20, col="red"))
      with(subset(data, abs(log2FC)>1), points(log2FC, -log10(pvalue), pch=20, col="orange"))
      with(subset(data, pvalue<.05 & abs(log2FC)>1), points(log2FC, -log10(pvalue), pch=20, col="green"))

      # Label points with the textxy function from the calibrate plot
      library(calibrate)
      with(subset(data, pvalue<.05 & abs(log2FC)>1), textxy(log2FC, -log10(pvalue), labs=gene_symbol, cex=.8))
      legend(x="topleft",
        legend=c("pvalue < 0.05", "log2FC > 1", "both"),
        col=c("red", "orange", "green"), box.lty=0, pch=20
      )
    } else {
      print('no data')
    }
  }
  
  output$volcano <- renderPlot({
    volcano()
  }, height=700, bg="transparent")

  ## download pca
  output$downloadVolcano <- downloadHandler(
    filename = "volcano.png",
    content = function(file) {
      png(file)
      volcano()
      dev.off()
    }
  )

  output$volcanoSelectedItems <- outputSelectedItems()

#########################################################################################################################################
###  Other Outputs
#########################################################################################################################################

  minPlotTheme <- theme(
    axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background=element_rect(fill = "transparent",colour = NA),
    plot.background=element_rect(fill = "transparent",colour = NA),
    axis.line = element_blank(),
    panel.border=element_blank(),plot.margin=margin(t = 0, unit='cm')
  )
  plot_color <- "#608fc3"
  plot_line_color <- "#26c487"
  
  ######################################### pvalue Plot 
  output$pvaluePlot <- renderPlot({
    ggplot(contigsDF, aes(pvalue)) +
    geom_density(fill = plot_color, color = NA) +
    xlim(minPvalues, maxPvalues) +
    geom_vline(xintercept = input$pvalue[[1]], color = plot_line_color) +
    geom_vline(xintercept = input$pvalue[[2]], color = plot_line_color) +
    minPlotTheme
  })

  ######################################### clipped 3p Plot 
  output$clipped3pPlot <- renderPlot({
    ggplot(contigsDF, aes(clipped_3p)) +
    geom_density(fill = plot_color, color = NA) +
    xlim(minClipped3ps, maxClipped3ps) +
    geom_vline(xintercept = input$clipped3p[[1]], color = plot_line_color) +
    geom_vline(xintercept = input$clipped3p[[2]], color = plot_line_color) +
    minPlotTheme
  })

  ######################################### nbSplice Plot 
  output$nbSplicePlot <- renderPlot({
    ggplot(contigsDF, aes(nb_splice)) +
    geom_density(fill = plot_color, color = NA) +
    xlim(minNbSplices, maxNbSplices) +
    geom_vline(xintercept = input$nbSplice[[1]], color = plot_line_color) +
    geom_vline(xintercept = input$nbSplice[[2]], color = plot_line_color) +
    minPlotTheme
  })

  ######################################### nbSnv Plot 
  output$nbSnvPlot <- renderPlot({
    ggplot(contigsDF, aes(nb_snv)) +
    geom_density(fill = plot_color, color = NA) +
    xlim(minNbSnvs, maxNbSnvs) +
    geom_vline(xintercept = input$nbSnv[[1]], color = plot_line_color) +
    geom_vline(xintercept = input$nbSnv[[2]], color = plot_line_color) +
    minPlotTheme
  })

  ######################################### nbHit Plot 
  output$nbHitPlot <- renderPlot({
    ggplot(contigsDF, aes(nb_hit)) +
    geom_density(fill = plot_color, color = NA) +
    xlim(minNbHits, maxNbHits) +
    geom_vline(xintercept = input$nbHit[[1]], color = plot_line_color) +
    geom_vline(xintercept = input$nbHit[[2]], color = plot_line_color) +
    minPlotTheme
  })

  ######################################### contigSize Plot 
  output$contigSizePlot <- renderPlot({
    ggplot(contigsDF, aes(contig_size)) +
    geom_density(fill = plot_color, color = NA) +
    xlim(minContigSizes, maxContigSizes) +
    geom_vline(xintercept = input$contigSize[[1]], color = plot_line_color) +
    geom_vline(xintercept = input$contigSize[[2]], color = plot_line_color) +
    minPlotTheme
  })

}
