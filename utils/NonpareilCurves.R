#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = F)
can_json <- FALSE

# Load stuff
if (!suppressPackageStartupMessages(
      requireNamespace("optparse", quietly = TRUE)))
  stop("R package 'optparse' is required.")
if (!suppressPackageStartupMessages(
      requireNamespace("Nonpareil", quietly = TRUE)))
  stop("R package 'Nonpareil' is required.")
if (suppressPackageStartupMessages(
      requireNamespace("jsonlite", quietly = TRUE)))
  can_json <- TRUE

optopt <- list(
  list(
    opt_str = "--labels", metavar = "SAMPLE1,SAMPLE2,...",
    help = "Labels to be used for the samples"
  ),
  list(
    opt_str = "--json", metavar = "FILE",
    help = paste0("Output file, Nonpareil curve processed data in JSON format",
             ifelse(can_json, "",
             "\n                (not available, requires jsonlite)"))
  ),
  list(
    opt_str = "--tsv", metavar = "FILE",
    help = "Output file, Nonpareil curve summaries in TSV format"
  ),
  list(
    opt_str = "--csv", metavar = "FILE",
    help = "Output file, Nonpareil curve summaries in CSV format"
  ),
  list(
    opt_str = "--pdf", metavar = "FILE",
    help = "Output file, Nonpareil curve plot in PDF format"
  ),
  list(
    opt_str = "--pdf-size", metavar = "WIDTH,HEIGHT",
    help = "Width and height of the PDF file (if --pdf)"
  ),
  list(
    opt_str = "--xlim", metavar = "FROM,TO",
    help = "Minimum and maximum points in the X-axis (in bp, if --pdf)"
  ),
  list(
    opt_str = "--col", metavar = "COL1,COL2,...",
    help = "Colors to be used for the Nonpareil curves (if --pdf)"
  ),
  list(
    opt_str = "--no-observed", action = "store_true",
    help = "Do not plot observed curves (if --pdf)"
  ),
  list(
    opt_str = "--no-model", action = "store_true",
    help = "Do not plot model curves (if --pdf)"
  )
)

opt <- optparse::parse_args(
  optparse::OptionParser(
    option_list = lapply(optopt, function(x) do.call(optparse::make_option, x)),
    description = list(),
    usage = 'usage: %prog [options] sample1.npo [sample2.npo ...]'
  ),
  positional_arguments = TRUE
)

if (length(opt[["args"]]) == 0)
  stop("You must pass at least one Nonpareil (npo) file, see -h")
opt.val <- opt[["options"]]
opt.val$files <- opt[["args"]]
opt <- opt.val

# Pad options with FALSE to simplify the code below (so I can avoid `is.null()`)
for (i in sapply(optopt, function(x) sub("^--", "", x[["opt_str"]])))
  opt[[paste0("do.", i)]] <- !is.null(opt[[i]])

# Do the work

##
# Function by @fgvieira via https://github.com/lmrodriguezr/nonpareil/issues/63,
# as modified and available in https://github.com/MultiQC/MultiQC with some
# minor additional modifications
export_curve <- function(object){
  # Extract variables
  n <- names(attributes(object))[c(1:12, 21:29)]
  x <- sapply(n, function(v) attr(object, v))
  names(x) <- n
  # Extract vectors
  n <- names(attributes(object))[13:20]
  y <- lapply(n, function(v) attr(object, v))
  names(y) <- n
  curve_json <- c(x, y)

  # Add model
  if (object$has.model) {
    x_min   <- 1e3
    x_max   <- signif(tail(attr(object, "x.adj"), n = 1) * 10, 1)
    x.model <- exp(seq(log(x_min), log(x_max), length.out = 1e3))
    y.model <- predict(object, lr = x.model)
    curve_json <- append(curve_json, list(x.model = x.model))
    curve_json <- append(curve_json, list(y.model = y.model))
  }

  curve_json
}

##
# Function by @fgvieira via https://github.com/lmrodriguezr/nonpareil/issues/63,
# also available in https://github.com/MultiQC/MultiQC
export_set <- function(object){
  y <- lapply(object$np.curves, "export_curve")
  names(y) <- sapply(object$np.curves, function(n) n$label)
  jsonlite::prettify(jsonlite::toJSON(y, auto_unbox = TRUE))
}

# Read curves
cat("Parsing Redundancy Files\n")
args <- list(opt[["files"]], plot = FALSE)
if (opt[["do.col"]])    args$col    <- strsplit(opt[["col"]], ",")[[1]]
if (opt[["do.labels"]]) args$labels <- strsplit(opt[["labels"]], ",")[[1]]
nps <- do.call(Nonpareil::Nonpareil.set, args)
if (opt[["do.json"]] && !can_json)
  stop("Requested JSON output but missing 'jsonlite' package")

# Write summaries
cat("Generating Summaries\n")
sum <- summary(nps)
if (opt[["do.json"]])
  writeLines(export_set(nps), con = opt[["json"]])
if (opt[["do.csv"]])
  write.csv(sum, file = opt[["csv"]])
if (opt[["do.tsv"]])
  write.table(sum, file = opt[["tsv"]], quote = FALSE, sep = "\t")

# Plot curves
if (opt[["do.pdf"]]) {
  cat("Plotting Curves\n")
  args <- list(opt[["pdf"]])
  if (opt[["do.pdf-size"]]) {
    wh <- as.numeric(strsplit(opt[["pdf-size"]], ",")[[1]])
    args$width  <- wh[1]
    args$height <- wh[2]
  }
  do.call(pdf, args)
  args <- list(nps)
  if (opt[["do.no-observed"]]) args$plot.observed <- FALSE
  if (opt[["do.no-model"]]) args$plot.model <- FALSE
  if (opt[["do.xlim"]])
    args$xlim <- as.numeric(strsplit(opt[["xlim"]], ",")[[1]])
  do.call(plot, args)
  t <- dev.off()
}

