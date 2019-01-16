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
  minDuPvalues = min(contigsDF$du_pvalue, na.rm = TRUE);         maxDuPvalues = max(contigsDF$du_pvalue, na.rm = TRUE); 
  minNbSplices = min(contigsDF$nb_splice, na.rm = TRUE);         maxNbSplices = max(contigsDF$nb_splice, na.rm = TRUE)
  minClipped3ps = min(contigsDF$clipped_3p, na.rm = TRUE);       maxClipped3ps = max(contigsDF$clipped_3p, na.rm = TRUE)
  minNbSnvs = min(contigsDF$nb_snv, na.rm = TRUE);               maxNbSnvs = max(contigsDF$nb_snv, na.rm = TRUE)
  minNbHits = min(contigsDF$nb_hit, na.rm = TRUE);               maxNbHits = max(contigsDF$nb_hit, na.rm = TRUE)
  minContigSizes = min(contigsDF$contig_size, na.rm = TRUE);     maxContigSizes = max(contigsDF$contig_size, na.rm = TRUE)
  
  # update UI inputs
  updateSliderInput(session, "pvalue", value = NULL, min = minPvalues, max = maxPvalues)
  updateSliderInput(session, "duPvalue", value = c(minDuPvalues, maxDuPvalues), min = minDuPvalues, max = maxDuPvalues)
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

  getFilterValues <- function() {
    f <- list(
      "inputPvalueMin"      = input$pvalue[[1]],
      "inputPvalueMax"      = input$pvalue[[2]],
      "inputDuPvalueMin"    = input$duPvalue[[1]],
      "inputDuPvalueMax"    = input$duPvalue[[2]],
      "inputNbSpliceMin"    = input$nbSplice[[1]],
      "inputNbSpliceMax"    = input$nbSplice[[2]],
      "inputClipped3pMin"   = input$clipped3p[[1]],
      "inputClipped3pMax"   = input$clipped3p[[2]],
      "inputNbSnvMin"       = input$nbSnv[[1]],
      "inputNbSnvMax"       = input$nbSnv[[2]],
      "inputNbHitMin"       = input$nbHit[[1]],
      "inputNbHitMax"       = input$nbHit[[2]],
      "inputContigSizeMin"  = input$contigSize[[1]],
      "inputContigSizeMax"  = input$contigSize[[2]],
      "customizedFilter"    = if (input$customFilter != "") input$customFilter else TRUE,
      "isMapped"            = if(input$isMapped != "NA") as.logical(input$isMapped) else c(FALSE, TRUE, NA),
      "geneIsDiff"          = if(input$geneIsDiff != "NA") as.logical(input$geneIsDiff) else c(FALSE, TRUE, NA),
      "hasGene"             = if(input$hasGene != "NA") as.logical(input$hasGene) else c(FALSE, TRUE),
      "hasASGene"           = if(input$hasASGene != "NA") as.logical(input$hasASGene) else c(FALSE, TRUE),
      "isExonic"            = if(input$isExonic != "NA") as.logical(input$isExonic) else c(FALSE, TRUE, NA),
      "isIntronic"          = if(input$isIntronic != "NA") input$isIntronic else c(FALSE, TRUE, NA)
    )
    return(f)
  }

  filterToString <- function() {
    f = getFilterValues()
    # s <- sprintf(
    #   "#pvalue: [%f, %f]\n#nb splice: [%f, %f]\n#is intronic: %s\n",
    #   f$inputPvalueMin, f$inputPvalueMax, f$inputNbSpliceMin, f$inputNbSpliceMax, input$isIntronic
    # )
    
    s <- sprintf(
      "#pvalue: [%f, %f]\n#nb splice: [%f, %f]\n#clipped 3p: [%f, %f]\n#nb snv: [%f, %f]\n#nb hit: [%f, %f]\n#contig size: [%f, %f]\n#du pvalue: [%f, %f]\n#is mapped: %s\n#gene is diff: %s\n#has gene: %s\n#has AS gene: %s\n#is exonic: %s\n#is intronic: %s\n#custom filter: %s",
      f$inputPvalueMin, f$inputPvalueMax, f$inputNbSpliceMin, f$inputNbSpliceMax, f$inputClipped3pMin, f$inputClipped3pMax, f$inputNbSnvMin, f$inputNbSnvMax, f$inputNbHitMin, f$inputNbHitMax, f$inputContigSizeMin, f$inputContigSizeMax, f$inputDuPvalueMin, f$inputDuPvalueMax, input$isMapped, input$geneIsDiff, input$hasGene, input$hasASGene, input$isExonic, input$isIntronic, f$customizedFilter
    )

    return(s)
  }

  ######################################  filter dataTable depending input filters
  dataTableFilters <- reactive({

    filterValues <- getFilterValues()

    cat(filterToString())
    
    # get subset with filters
    s <- subset(
      contigsDF,
      ((pvalue      >= filterValues$inputPvalueMin & pvalue           <= filterValues$inputPvalueMax)     | (input$pvalueKeepNA & is.na(pvalue))) &
      ((du_pvalue   >= filterValues$inputDuPvalueMin & du_pvalue      <= filterValues$inputDuPvalueMax)   | (input$duPvalueKeepNA & is.na(du_pvalue))) &
      ((nb_splice   >= filterValues$inputNbSpliceMin & nb_splice      <= filterValues$inputNbSpliceMax)   | (input$nbSpliceKeepNA & is.na(nb_splice))) &
      ((clipped_3p  >= filterValues$inputClipped3pMin & clipped_3p    <= filterValues$inputClipped3pMax)  | (input$clipped3pKeepNA & is.na(clipped_3p))) &
      ((nb_snv      >= filterValues$inputNbSnvMin & nb_snv            <= filterValues$inputNbSnvMax)      | (input$nbSnvKeepNA & is.na(nb_snv))) &
      ((nb_hit      >= filterValues$inputNbHitMin & nb_hit            <= filterValues$inputNbHitMax)      | (input$nbHitKeepNA & is.na(nb_hit))) &
      ((contig_size >= filterValues$inputContigSizeMin & contig_size  <= filterValues$inputContigSizeMax) | (input$contigSizeKeepNA & is.na(contig_size))) &
      is.element(gene_is_diff,        filterValues$geneIsDiff) &
      is.element(!is.na(gene_id),     filterValues$hasGene) &
      is.element(!is.na(as_gene_id),  filterValues$hasASGene) &
      is.element(is_mapped,           filterValues$isMapped) &
      is.element(exonic,              filterValues$isExonic) &
      is.element(intronic,            filterValues$isIntronic) &
      eval(parse(text =               filterValues$customizedFilter))
    )

    return(s)
  })
  
  ####################################### observe only profil filters, get values from filters preset file
  observe({
    if (input$preset != '-') {
      filterPreset <- filtersPresetDF[which(filtersPresetDF$event == input$preset),]

      inputPvalueMin      <- minPvalues
      inputPvalueMax      <- maxPvalues
      inputIsMapped       <- filterPreset$is_mapped
      inputGeneIsDiff     <- filterPreset$gene_is_diff
      inputHasGene        <- filterPreset$has_gene
      inputHasASGene      <- filterPreset$has_as_gene
      inputIsIntronic     <- filterPreset$intronic
      inputIsExonic       <- filterPreset$exonic
      inputDuPvalueMin    <- if (!is.na(filterPreset$du_pvalue_min)) filterPreset$du_pvalue_min else minDuPvalues
      inputDuPvalueMax    <- if (!is.na(filterPreset$du_pvalue_max)) filterPreset$du_pvalue_max else maxDuPvalues
      inputNbSpliceMin    <- if (!is.na(filterPreset$nb_splice_min)) filterPreset$nb_splice_min else minNbSplices
      inputNbSpliceMax    <- if (!is.na(filterPreset$nb_splice_max)) filterPreset$nb_splice_max else maxNbSplices
      inputClipped3pMin   <- if (!is.na(filterPreset$clipped_3p_min)) filterPreset$clipped_3p_min else minClipped3ps
      inputClipped3pMax   <- if (!is.na(filterPreset$clipped_3p_max)) filterPreset$clipped_3p_max else maxClipped3ps
      inputNbSnvMin       <- if (!is.na(filterPreset$nb_snv_min)) filterPreset$nb_snv_min else minNbSnvs
      inputNbSnvMax       <- if (!is.na(filterPreset$nb_snv_max)) filterPreset$nb_snv_max else maxNbSnvs
      inputNbHitMin       <- if (!is.na(filterPreset$nb_hit_min)) filterPreset$nb_hit_min else minNbHits
      inputNbHitMax       <- if (!is.na(filterPreset$nb_hit_max)) filterPreset$nb_hit_max else maxNbHits
      inputContigSizeMin  <- if (!is.na(filterPreset$contig_size_min)) filterPreset$contig_size_min else minContigSizes
      inputContigSizeMax  <- if (!is.na(filterPreset$contig_size_max)) filterPreset$contig_size_max else maxContigSizes

      
      cat(sprintf(
        "\n\nProfil filter: %s\ngene is diff: %s\ngene symbol: %s\nas gene id: %s\nis mapped: %s\nexonic: %s\nintronic: %s\nother: %s\n",
        input$preset, filterPreset$gene_is_diff, filterPreset$gene_symbol, filterPreset$as_gene_id, filterPreset$is_mapped, filterPreset$exonic, filterPreset$intronic, filterPreset$other
      ))

      # update UI slider inputs
      updateSliderInput(session, "pvalue", value = c(inputPvalueMin, inputPvalueMax))
      updateSliderInput(session, "duPvalue", value = c(inputDuPvalueMin, inputDuPvalueMax))
      updateSliderInput(session, "nbSplice", value = c(inputNbSpliceMin, inputNbSpliceMax))
      updateSliderInput(session, "clipped3p", value = c(inputClipped3pMin, inputClipped3pMax))
      updateSliderInput(session, "nbSnv", value = c(inputNbSnvMin, inputNbSnvMax))
      updateSliderInput(session, "nbHit", value = c(inputNbHitMin, inputNbHitMax))
      updateSliderInput(session, "contigSize", value = c(inputContigSizeMin, inputContigSizeMax))
      updateTextInput(session,   "customFilter", value = filterPreset$other)

      # Binary operations
      updateRadioButtons(session, "isIntronic", selected = toString(inputIsIntronic))
      updateRadioButtons(session, "isExonic",   selected = toString(inputIsExonic))
      updateRadioButtons(session, "hasASGene",  selected = toString(inputHasASGene))
      updateRadioButtons(session, "hasGene",    selected = toString(inputHasGene))
      updateRadioButtons(session, "geneIsDiff", selected = toString(inputGeneIsDiff))
      updateRadioButtons(session, "isMapped",   selected = toString(inputIsMapped))
    }
  })

  ####################################### observe only numeric inputs to update sliders 
  observe({
    updateSliderInput(session, "pvalue", value = c(input$minPvalue, input$maxPvalue))
    updateSliderInput(session, "duPvalue", value = c(input$minDuPvalues, input$maxDuPvalues))
    updateSliderInput(session, "clipped3p", value = c(input$minClipped3p, input$maxClipped3p))
    updateSliderInput(session, "nbSplice", value = c(input$minNbSplice, input$maxNbSplice))
    updateSliderInput(session, "nbSnv", value = c(input$minNbSnv, input$maxNbSnv))
    updateSliderInput(session, "nbHit", value = c(input$minNbHit, input$maxNbHit))
    updateSliderInput(session, "contigSize", value = c(input$minContigSize, input$maxContigSize))
  })

  ####################################### reset
  observeEvent(input$reset, {
    updateSliderInput(session, "pvalue", value = c(minPvalues, maxPvalues));             updateCheckboxInput(session, "pvalueKeepNA", value = TRUE)
    updateSliderInput(session, "duPvalue", value = c(minDuPvalues, maxDuPvalues));        updateCheckboxInput(session, "duPvalueKeepNA", value = TRUE)
    updateSliderInput(session, "clipped3p", value = c(minClipped3ps, maxClipped3ps));     updateCheckboxInput(session, "clipped3pKeepNA", value = TRUE)
    updateSliderInput(session, "nbSplice", value = c(minNbSplices, maxNbSplices));        updateCheckboxInput(session, "nbSpliceKeepNA", value = TRUE)
    updateSliderInput(session, "nbSnv", value = c(minNbSnvs, maxNbSnvs));                 updateCheckboxInput(session, "nbSnvKeepNA", value = TRUE)
    updateSliderInput(session, "nbHit", value = c(minNbHits, maxNbHits));                 updateCheckboxInput(session, "nbHitKeepNA", value = TRUE)
    updateSliderInput(session, "contigSize", value = c(minContigSizes, maxContigSizes));  updateCheckboxInput(session, "contigSizeKeepNA", value = TRUE)
    updateTextInput(session, "customFilter", value = "")
    updateSelectInput(session, "preset", selected = "-")

    # Reset binary filters
    updateRadioButtons(session, "isMapped", selected = "NA")
    updateRadioButtons(session, "geneIsDiff", selected = "NA")
    updateRadioButtons(session, "hasGene", selected = "NA")
    updateRadioButtons(session, "hasASGene", selected = "NA")
    updateRadioButtons(session, "isIntronic", selected = "NA")
    updateRadioButtons(session, "isExonic", selected = "NA")
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

  ####################################### Download table
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("selected_contigs",".tsv", sep = "")
    },
    content = function(file) {
      con <- file(file, 'w')
      writeLines(filterToString(), con = file)
      write.table(dataTableFilters(), file = file, append=TRUE, sep = "\t", col.names = TRUE, row.names = FALSE)
    }
  )

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
        )),
        dom = 'Bfrtip'
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
      colnames(mat_scaled) = formatSample
      # set conditions color. Caution: must not have more than 5 conditions (2 usually).
      # this is to transform condition to a named vector
      uniqCdts = unique(conditions)
      colorCdts = c("black", "green", "blue", "red", "yellow")
      colorCdtsVec = list(matrix(colorCdts[1: length(uniqCdts)], nrow = 1, ncol = length(uniqCdts), byrow = TRUE, dimnames = list(c(), uniqCdts)))
      colorCdtsVec = unlist(as.data.frame(colorCdtsVec)[1,])
      names(colorCdtsVec) = uniqCdts
      
      haConditions = HeatmapAnnotation(df = data.frame(conditions = conditions), col = list(conditions = colorCdtsVec), show_annotation_name = TRUE)

      Heatmap(
        mat_scaled, name = "expression", km = 5, col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        top_annotation = haConditions,
        cluster_columns = TRUE,
        # column_title_gp = gpar(),
        show_row_names = FALSE,
        show_column_names = TRUE
      ) +
      Heatmap(base_mean, name = "base mean", col = colorRamp2(c(0, 1000), c("white", "grey")), show_row_names = FALSE, width = unit(5, "mm")) +
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
    filename = "heatmap.pdf",
    content = function(file) {
      pdf(file=file, paper="a4r", width=30)
      h = heatmap()
      HM <- draw(h)
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

      # Add colored points: red if pvalue (padj ?) <0.05, orange of log2FC>1, gray if both)
      with(subset(data, pvalue<.05 ), points(log2FC, -log10(pvalue), pch=20, col="red"))
      with(subset(data, abs(log2FC)>1), points(log2FC, -log10(pvalue), pch=20, col="orange"))
      with(subset(data, pvalue<.05 & abs(log2FC)>1), points(log2FC, -log10(pvalue), pch=20, col="gray"))

      # Label points with the textxy function from the calibrate plot
      library(calibrate)
      with(subset(data, pvalue<.05 & abs(log2FC)>1), textxy(log2FC, -log10(pvalue), labs=gene_symbol, cex=.8))
      legend(x="topleft",
        legend=c("pvalue < 0.05", "log2FC > 1", "both"),
        col=c("red", "orange", "gray"), box.lty=0, pch=20
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

  ######################################### Du pvalue Plot 
  output$duPvaluePlot <- renderPlot({
    ggplot(contigsDF, aes(du_pvalue)) +
    geom_density(fill = plot_color, color = NA) +
    xlim(minDuPvalues, maxDuPvalues) +
    geom_vline(xintercept = input$duPvalue[[1]], color = plot_line_color) +
    geom_vline(xintercept = input$duPvalue[[2]], color = plot_line_color) +
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
