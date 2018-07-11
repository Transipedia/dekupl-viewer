#########################################################################################################################################
###  Define args options
#########################################################################################################################################

library(optparse)

option_list = list(
  make_option(
    c("-c", "--contig"), type="character", default=NULL, help="contigs differential expression file (DiffContigsInfos.tsv)"
  ),
  make_option(
    c("-s", "--sample"), type="character", default=NULL, help="sample conditions file (sample_conditions_full.tsv)"
  ),
  make_option(
    c("-t", "--test"), action="store_true", default=NULL, help="run app with a test dataset"
  )
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$test) & (is.null(opt$contig) | is.null(opt$sample))) {
  print_help(opt_parser)
  stop("-c and -s arguments must be supplied (input files).\n", call.=FALSE)
}

if (!is.null(opt$test)) {
  opt$contig <- '../toy/DiffContigsInfos.tsv'
  opt$sample <- '../toy/sample_conditions_full.tsv'
  print("running app with test dataset")
}

#########################################################################################################################################
###  check input files
#########################################################################################################################################

checkInputfiles <- function() {
  if (!file.exists(opt$contig)) {
    stop("file not found ", opt$contig, ".\n", call.=FALSE)
  }
  if (!file.exists(opt$sample)) {
    stop("file not found ", opt$sample, ".\n", call.=FALSE)
  }
}

checkInputfiles()

#########################################################################################################################################
###  run the app
#########################################################################################################################################

library(shiny)
library(DT)
library(ggplot2)
library(shinydashboard)
library(ComplexHeatmap)
library(circlize)

app <- shinyAppDir("./")
runApp(app, 8080, FALSE, '0.0.0.0')
