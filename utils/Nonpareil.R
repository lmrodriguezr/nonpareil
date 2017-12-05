# Classes
setClass("Nonpareil.Curve",
  ### A single Nonpareil curve. This object can be produced by `Nonpareil.curve`
  ### and supports S4 methods `plot`, `summary`, `print`, and `predict`
  representation(
  # Dataset info
  file='character',      ##<< Input .npo file
  label='character',     ##<< Name of the dataset
  col='character',       ##<< Color of the dataset
  # Input .npo metadata
  L='numeric',           ##<< Read length
  AL='numeric',          ##<< Adjusted read length (same as L for alignment)
  R='numeric',           ##<< Number of reads
  LR='numeric',          ##<< Effective sequencing effort used
  overlap='numeric',     ##<< Minimum read overlap
  ksize='numeric',       ##<< K-mer size (for kmer kernel only)
  log.sample='numeric',  ##<< Multiplier of the log-sampling (or zero if linear)
  kernel='character',    ##<< Read-comparison kernel
  version='character',   ##<< Nonpareil version used
  # Input .npo data
  x.obs='numeric',       ##<< Rarefied sequencing effort
  x.adj='numeric',       ##<< Adjusted rarefied sequencing effort
  y.red='numeric',       ##<< Rarefied redundancy (observed)
  y.cov='numeric',       ##<< Rarefied coverage (corrected)
  y.sd= 'numeric',       ##<< Standard deviation of rarefied coverage
  y.p25='numeric',       ##<< Percentile 25 (1st quartile) of rarefied coverage
  y.p50='numeric',       ##<< Percentile 50 (median) of rarefied coverage
  y.p75='numeric',       ##<< Percentile 75 (3rd quartile) of rarefied coverage
  # Estimated coverage
  kappa='numeric',       ##<< Dataset redundancy
  C='numeric',           ##<< Dataset coverage
  consistent='logical',  ##<< Is the data sufficient for accurate estimation?
  # Projected coverage
  star='numeric',        ##<< Coverage considered 'nearly complete'
  has.model='logical',   ##<< Was the model successfully estimated?
  warning='character',   ##<< Warnings generated on consistency or model fitting
  LRstar='numeric',      ##<< Projected seq. effort for nearly complete coverage
  modelR='numeric',      ##<< Pearson's R for the estimated model
  diversity='numeric',   ##<< Dataset Nd index of sequence diversity
  model='list',          ##<< Fitted sigmoidal model
  # Call
  call='call'            ##<< Call producing this object
  ), package='Nonpareil'
  );
setClass("Nonpareil.Set",
  ### Collection of `Nonpareil.Curve` objects. This object can be produced by
  ### `Nonpareil.curve.batch` and supports S4 methods `plot`, `summary`, and
  ### `print`
  representation(
  np.curves='list',      ##<< List of `Nonpareil.Curve` objects
  call='call'            ##<< Call producing this object
  ), package='Nonpareil'
  );

# S3 Methods
setMethod("$", "Nonpareil.Curve", function(x, name) attr(x, name))
setMethod("$<-", "Nonpareil.Curve",
      function(x, name, value) { attr(x, name) <- value ; x })
setMethod("$", "Nonpareil.Set", function(x, name) attr(x, name))
setMethod("$<-", "Nonpareil.Set",
      function(x, name, value) { attr(x, name) <- value ; x })

plot.Nonpareil.Set <- function(
      ### Plot a `Nonpareil.Set` object
      x,
      ### `Nonpareil.Set` object to plot
      col=NA,
      ### Color of the curves (vector). If passed, it overrides the colors set
      ### in the `Nonpareil.Curve` objects. Values are recycled
      labels=NA,
      ### Labels of the curves (vector). If passed, it overrides the labels set
      ### in the `Nonpareil.Curve` objects. Values are recycled
      main="Nonpareil Curves",
      ### Title of the plot
      legend.opts=list(),
      ### Any additional parameters passed to `Nonpareil.legend`. If FALSE, the
      ### legend is not displayed
      ...
      ### Any additional parameters passed to `plot.Nonpareil.Curve`
      ){
  if(!inherits(x, "Nonpareil.Set"))
    stop("'x' must inherit from class `Nonpareil.Set`")

  # Plots
  new <- TRUE;
  col <- rep(col, length.out=length(x$np.curves))
  labels <- rep(labels, length.out=length(x$np.curves))
  for(i in 1:length(x$np.curves)){
    if(!is.na(col[i])) x$np.curves[[i]]$col <- col[i]
    if(!is.na(labels[i])) x$np.curves[[i]]$label <- labels[i]
    plot(x$np.curves[[i]], new=new, main=ifelse(new, main, ""), ...)
    new <- FALSE
  }

  # Legend
  if(inherits(legend.opts, "list")){
    legend.opts[["np"]] <- x
    if(is.null(legend.opts[["x"]])) legend.opts[["x"]] <- "bottomright"
    do.call(Nonpareil.legend, legend.opts)
  }

  # Return
  invisible(x)
  ### Returns invisibly a `Nonpareil.Set` object (same as `x` input)
}
plot.Nonpareil.Curve <- function(
      ### Plot a `Nonpareil.Curve` object
      x,
      ### `Nonpareil.Curve` object to plot
      col=NA,
      ### Color of the curve. If passed, it overrides the colors set in the
      ### `Nonpareil.Curve` object
      add=FALSE,
      ### If TRUE, it attempts to use a previous (active) canvas to plot the
      ### curve
      new=!add,
      ### Inverse of `add`
      plot.observed=TRUE,
      ### Indicates if the observed (rarefied) coverage is to be plotted
      plot.model=TRUE,
      ### Indicates if the fitted model is to be plotted
      plot.dispersion=FALSE,
      ### Indicates if (and how) dispersion of the replicates should be plotted.
      ### Supported values are:
      ### FALSE: no dispersion is plotted (default),
      ### 'sd': one standard deviation around the mean,
      ### 'ci95': 95% confidence interval,
      ### 'ci90': 90% confidence interval,
      ### 'ci50': 50% confidence interval,
      ### 'iq': Inter-quartile range
      plot.diversity=TRUE,
      ### If TRUE, the diversity estimate is plotted as a small arrow below the
      ### Nonpareil curve
      xlim=c(1e3,10e12),
      ### Limits of the sequencing effort (X-axis)
      ylim=c(1e-6,1),
      ### Limits of the coverage (Y-axis)
      main=paste("Nonpareil Curve for", x$label),
      ### Title of the plot
      xlab="Sequencing effort (bp)",
      ### Label of the X-axis
      ylab="Estimated Average Coverage",
      ### Label of the Y-axis
      curve.lwd=2,
      ### Line width of the rarefied coverage
      curve.alpha=0.4,
      ### Alpha value (from 0 to 1) of the rarefied coverage
      model.lwd=1,
      ### Line width of the model
      model.alpha=1,
      ### Alpha value (from 0 to 1) of the model
      log="x",
      ### Axis to plot in logarithmic scale. Supported values are:
      ### 'x': sequencing effort (default),
      ### 'y': coverage,
      ### 'xy': both logarithmic, or
      ### '': both linear
      arrow.length=0.05,
      ### If plot.diversity=TRUE, it determines the length of the arrow to
      ### display the divsersity (as a fraction of the ylim range).
      arrow.head=arrow.length,
      ### If plot.diversity=TRUE, it determines the length of the arrow head to
      ### display the diversity index (in inches).
      ...
      ### Additional graphical parameters
      ){
  if(!inherits(x, "Nonpareil.Curve"))
    stop("'x' must inherit from class `Nonpareil.Curve`")

  # Create empty canvas
  if(new){
    plot(1, type="n", xlim=xlim, ylim=ylim, bty="l",
          xlab=xlab, ylab=ylab, main=main, xaxs="i", yaxs="i", log=log, ...)
    abline(h=c(1, x$star/100), lty=2, col="red")
    abline(v=10^seq(0,15,by=3), lty=2, col="gray80")
  }

  # Dispersion
  if(plot.dispersion!=FALSE){
    if(plot.dispersion == "sd"){
      err.y <- c(x$y.cov+x$y.sd, rev(x$y.cov-x$y.sd))
    }else if(plot.dispersion == "ci95"){
      err.y <- c(x$y.cov+x$y.sd*1.9, rev(x$y.cov-x$y.sd*1.9))
    }else if(plot.dispersion == "ci90"){
      err.y <- c(x$y.cov+x$y.sd*1.64, rev(x$y.cov-x$y.sd*1.64))
    }else if(plot.dispersion == "ci50"){
      err.y <- c(x$y.cov+x$y.sd*.67, rev(x$y.cov-x$y.sd*.67))
    }else if(plot.dispersion == "iq"){
      err.y <- c(x$y.p25, rev(x$y.p75))
    }
    polygon(c(x$x.adj, rev(x$x.adj)), border=NA,
      ifelse(err.y<=ylim[1]*0.1, ylim[1]*0.1, err.y), col=Nonpareil.col(x, .2))
  }

  # Rarefied coverage
  if(plot.observed){
    lines(x$x.adj, x$y.cov, col=Nonpareil.col(x, curve.alpha), lwd=curve.lwd);
  }

  # Model
  if(x$has.model & plot.model){
    model.lty <- ifelse(plot.observed, 2, 1)
    model.x   <- exp(seq(log(xlim[1]), log(xlim[2]), length.out=1e3));
    model.y   <- predict(x, lr=model.x);
    lines(model.x, model.y, col=Nonpareil.col(x, model.alpha), lty=model.lty,
          lwd=model.lwd)
    if(!plot.observed){
      points(x$LR, predict(x), col=Nonpareil.col(x, 1.0), pch=21, bg="white")
    }
  }

  if(x$has.model & plot.diversity & x$diversity>0){
    arrows(x0=exp(x$diversity), length=arrow.head,
          y1=ifelse(log=='y' | log=='xy' | log=='yx',
            ylim[1]*(ylim[2]/ylim[1])**arrow.length, ylim[1] + diff(ylim)*arrow.length),
          y0=ylim[1], col=Nonpareil.col(x, model.alpha));
  }

  # Return
  invisible(x)
  ### Retuns invisibly a `Nonpareil.Curve` object (same as `x` input)
}
summary.Nonpareil.Set <- function(
      ### Returns a summary of the Nonpareil.Set results
      object,
      ### `Nonpareil.Set` object
      ...
      ### Additional parameters ignored
      ){
  if(!inherits(object, "Nonpareil.Set"))
    stop("'object' must inherit from class `Nonpareil.Set`")
  y <- rbind(sapply(object$np.curves, "summary"))
  colnames(y) <- sapply(object$np.curves, function(n) n$label)
  t(y)
}
summary.Nonpareil.Curve <- function(
      ### Returns a summary of the Nonpareil.Curve results
      object,
      ### `Nonpareil.Curve` object
      ...
      ### Additional parameters ignored
      ){
  if(!inherits(object, "Nonpareil.Curve"))
    stop("'object' must inherit from class `Nonpareil.Curve`")
  n <- c("kappa","C","LR","modelR","LRstar","diversity")
  y <- sapply(n, function(v) attr(object,v))
  names(y) <- n
  y
}
print.Nonpareil.Set <- function(
      ### Prints and returns invisibly a summary of the `Nonpareil.Set` results
      x,
      ### `Nonpareil.Set` object
      ...
      ### Additional parameters ignored
      ){
  if(!inherits(x, "Nonpareil.Set"))
    stop("'x' must inherit from class `Nonpareil.Set`")
  y <- summary(x)
  cat("===[ Nonpareil.Set ]===================================\n")
  cat("Collection of", length(x$np.curves), "Nonpareil curves.\n")
  print(y)
  cat("-------------------------------------------------------\n")
  cat("call:",as.character(x$call),"\n")
  cat("-------------------------------------------------------\n")
  invisible(y)
}
print.Nonpareil.Curve <- function(
      ### Prints and returns invisibly a summary of the `Nonpareil.Curve`
      ### results
      x,
      ### `Nonpareil.Set` object
      ...
      ### Additional parameters ignored
      ){
  if(!inherits(x, "Nonpareil.Curve"))
    stop("'x' must inherit from class `Nonpareil.Curve`")
  y <- summary(x)
  yp <- cbind(y)
  colnames(yp) <- x$label
  cat("===[ Nonpareil.Curve ]=================================\n")
  print(yp)
  cat("-------------------------------------------------------\n")
  cat("call:",as.character(x$call),"\n")
  cat("-------------------------------------------------------\n")
  invisible(y)
}
predict.Nonpareil.Curve <- function(
      ### Predict the coverage for a given sequencing effort
      object,
      ### `Nonpareil.Curve` object
      lr=object$LR,
      ### Sequencing effort for the prediction (in bp)
      ...
      ### Additional parameters ignored
      ){
  if(!inherits(object, "Nonpareil.Curve"))
    stop("'object' must inherit from class `Nonpareil.Curve`")
  if(!object$has.model)
    stop("'object' must be a Nonpareil Curve with a fitted model")
  predict(object$model, list(x=lr))
}

# Ancillary functions
Nonpareil.read_metadata <- function(
      ### Read the metadata headers
      x
      ### `Nonpareil.Curve` object
      ){
  # Load key-values and defaults
  meta_data <- gsub('^# @', "", grep("^# @", readLines(x$file), value=TRUE))
  keys <- gsub(': .*', "", meta_data)
  vals <- gsub('.*: ', "", meta_data)
  x$kernel <- "alignment"
  x$log.sample <- 0

  # Set metadata
  if("ksize" %in% keys && vals[keys=="ksize"]>0) x$kernel <- "kmer";
  if("divide" %in% keys) x$log.sample <- as.numeric(vals[keys=="divide"]);
  if("logsampling" %in% keys) x$log.sample <- as.numeric(vals[keys=="logsampling"]);
  x$version   <- as.numeric(vals[keys=="version"])
  x$L         <- as.numeric(vals[keys=="L"])
  x$R         <- as.numeric(vals[keys=="R"])
  if(x$kernel=="kmer"){
    x$overlap <- 50
    x$ksize   <- as.numeric(vals[keys=="ksize"])
    x$AL      <- as.numeric(vals[keys=="AL"])
  }else{
    x$overlap <- as.numeric(vals[keys=="overlap"])
    x$AL      <- x$L
  }
  x$LR <- exp(log(x$R) + log(x$L));

  invisible(x)
}
Nonpareil.read_data <- function(
      ### Read the data tables and extract direct estimates
      x,
      ### `Nonpareil.Curve` object
      correction.factor
      ### Logical; see `Nonpareil.curve` for details
      ){
  # Read input
  a <- read.table(x$file, sep="\t", header=FALSE)
  x$x.obs <- a[,1]
  x$y.red <- a[,2]
  x$kappa <- tail(x$y.red, n=1)

  # Estimate coverage
  cor.f   <- 1.0;
  if(correction.factor) cor.f <- Nonpareil.coverage_factor(x)
  for(i in 2:6) a[, i] <- a[, i]^cor.f;
  x$y.cov <- a[, 2]
  x$y.sd  <- a[, 3]
  x$y.p25 <- a[, 4]
  x$y.p50 <- a[, 5]
  x$y.p75 <- a[, 6]
  x$C <- tail(x$y.cov, n=1)

  # Adjust sequencing effort
  x$x.adj <- exp(
    max(log(x$x.obs)) + (x$C^0.27)*(log(x$x.obs) - max(log(x$x.obs)))
  )
  # Obsolete corrections {
  #   x$x.adj <- exp(log(x$x.adj)*0.61 + 10)
  #   x$x.adj <- x$x.adj * x$AL / 101
  # }
  x$x.adj <- x$x.adj * x$AL * x$R / max(x$x.adj)

  # Check consistency
  x$consistent <- TRUE
  twenty.pc = which.max(x$x.adj[x$x.adj <= 0.5*tail(x$x.adj, n=1)]);
  if(x$y.p50[twenty.pc]==0){
    x$consistent <- FALSE
    x$warnings <- c(x$warnings,
          paste("Median of the curve is zero at 50% of the reads, check",
          "parameters and re-run (e.g., decrease -L in nonpareil -T alignment)."))
  }
  if(x$kappa <= 1e-5){
    x$consistent <- FALSE
    x$warnings <- c(x$warnings,
          paste("Redundancy curve too low, check parameters and re-run",
          "(e.g., decrease -L in nonpareil -T alignment)."))
  }
  if(x$y.cov[2]>=1-1e-5){
    x$consistent <- FALSE
    x$warnings <- c(x$warnings,
          paste("Curve too steep, check parameters and re-run",
          "(e.g., increase value of -d in nonpareil)."))
  }
  if(sum(x$y.cov>0 & x$y.cov<0.9) <= 10){
    x$consistent <- FALSE
    x$warnings <- c(x$warnings,
          paste("Insufficient resolution below 90% coverage, check",
          "parameters and re-run (e.g., increase the value of -d in nonpareil)."))
  }

  invisible(x)
}
Nonpareil.fit_model <- function(
      ### Fit the sigmoidal model to the rarefied coverage
      np,
      ### `Nonpareil.Curve` object
      weights.exp
      ### Numeric; See `Nonpareil.curve` for details
      ){
  if(!inherits(np, "Nonpareil.Curve"))
    stop("'np' must inherit from class `Nonpareil.Curve`")

  # Prepare data
  sel <- np$y.cov>0 & np$y.cov<0.9
  data <- list(x=np$x.adj[ sel ], y=np$y.cov[ sel ])
  if(is.na(weights.exp[1])){
    if(np$log.sample==0){ weights.exp <- c(-1.1,-1.2,-0.9,-1.3,-1) }
    else{ weights.exp <- c(0,1,-1,1.3,-1.1,1.5,-1.5,3,-3) }
  }

  # Find the first weight with proper fit
  np$has.model <- FALSE
  weights.i <- 0
  while(!np$has.model & !is.na(weights.exp[weights.i+1])){
    weights.i <- weights.i+1
    suppressWarnings(
      model <- nls(y ~ Nonpareil.f(x, a, b), data=data,
            weights=(np$y.sd[sel]^weights.exp[weights.i]),
            start=list(a=1, b=0.1), lower=c(a=0, b=0), algorithm="port",
             control=nls.control(
                   minFactor=1e-25000, tol=1e-15, maxiter=1024, warnOnly=TRUE))
    )
    tryCatch({ is.conv <- summary(model)$convInfo$isConv },
          error=function(e){ is.conv <- FALSE })
    if(is.conv){
      np$model <- model
      np$has.model <- TRUE
    }
  }

  # Estimate diversity and projections
  if(np$has.model){
    pa <- coef(np$model)[1]
    pb <- coef(np$model)[2]
    if(pa > 1) np$diversity <- (pa-1)/pb
    np$LRstar <- Nonpareil.antif(np$star/100, pa, pb)
    np$modelR <- cor(data$y, predict(np, lr=data$x))
  }else{
    np$warnings <- c(np$warnings,
          "Model didn't converge. Try modifying the values of weights.exp.")
  }

  invisible(np)
}
Nonpareil.coverage_factor <- function(
      ### Factor to transform redundancy into coverage (internal function).
      x
      ### `Nonpareil.Curve` object
      ){
  return(1 - exp(2.23E-2 * x$overlap - 3.5698))
  ### A numeric scalar.
}
Nonpareil.col <- function(
      ### Returns the color of the curve
      x,
      ### `Nonpareil.Curve` or `Nonpareil.Set` object
      alpha=1
      ### Alpha level of the color from 0 to 1
      ){
  if(inherits(x, "Nonpareil.Curve")){
    col <- x$col
  }else if(!inherits(x, "Nonpareil.Set")){
    col <- sapply(x$np.curves, function(np) np$col)
  }else{
    stop("'x' must inherit from class `Nonpareil.Curve` or `Nonpareil.Set`")
  }
  apply(col2rgb(col), 2, function(x) do.call(rgb, as.list(c(x[1:3]/256, alpha))) )
}
Nonpareil.legend <- function(
      ### Generates a legend for Nonpareil plots
      np,
      ### A `Nonpareil.Set` object or a list of `Nonpareil.Curve` objects
      x,
      ### X coordinate, or any character string accepted by legend (e.g.,
      ### 'bottomright').
      y=.3,
      ### Y coordinate.
      ...
      ### Any other parameters supported by legend().
      ){
  if(inherits(np, "Nonpareil.Set")) np <- np$np.curves
  if(!inherits(np, "list"))
    stop("'np' must inherit from `list` or class `Nonpareil.Set`")

  labels <- sapply(np, function(x) x$label)
  cols <- sapply(np, Nonpareil.col)
  legend(x=x, y=y, legend=labels, fill=cols, ...);
}
Nonpareil.add.curve <- function(
      ### Adds a `Nonpareil.Curve` to a `Nonpareil.Set`
      nps,
      ### `Nonpareil.Set` object
      np
      ### `Nonpareil.Curve` object
      ){
  if(!inherits(nps, "Nonpareil.Set"))
    stop("'nps' must inherit from class `Nonpareil.Set`")
  if(!inherits(np, "Nonpareil.Curve"))
    stop("'np' must inherit from class `Nonpareil.Curve`")

  nps$np.curves[[ length(nps$np.curves)+1 ]] <- np
  return(nps)
}
setMethod("+", "Nonpareil.Set", function(e1,e2) Nonpareil.add.curve(e1,e2))

# Model Functions
Nonpareil.f <- function(
      ### Function of the projected model
      x,
      ### Values of sequencing effort (typically in bp)
      a,
      ### Parameter alpha of the Gamma CDF
      b
      ### Parameter beta of the Gamma CDF
      ){
  return(pgamma(log1p(x), a, b))
  ### Predicted values of abundance-weighted average coverage.
}
Nonpareil.antif <- function(
    ### Complement function of `Nonpareil.f`
    y,
    ### Values of abundance-weighted average coverage
    a,
    ### Parameter alpha of the gamma CDF.
    b
    ### Parameter beta of the gamma CDF.
    ){
  return(exp(qgamma(y,a,b))-1)
  ### Estimated sequencing effort.
}

# Main functions
Nonpareil.curve <- structure(function(
      ### Generates a Nonpareil curve from an .npo file
      file,
      ### Path to the .npo file, containing the read redundancy
      plot=TRUE,
      ### Determines if the plot should be produced. If FALSE, it still computes
      ### the coverage and the model
      label=NA,
      ### Name of the dataset. If NA, it is determined by the file name
      col=NA,
      ### Color of the curve. If NA, a random color is assigned (even if
      ### plot=FALSE),
      enforce.consistency=TRUE,
      ### If TRUE, it fails verbosely on insufficient data, otherwise it warns
      ### about the inconsistencies and attempts the estimations
      star=95,
      ### Objective coverage in percentage; i.e., coverage value considered
      ### near-complete
      correction.factor=TRUE,
      ### Should the overlap-dependent (or kmer-length-dependent) correction
      ### factor be applied? If FALSE, redundancy is assumed to equal coverage.
      weights.exp=NA,
      ### Vector of values to be tested (in order) as exponent of the weights
      ### distribution. If the model fails to converge, sometimes manual
      ### modifications in this parameter may help. By default (NA), five
      ### different values are tested in the following order: For linear
      ### sampling, -1.1, -1.2, -0.9, -1.3, -1. For logarithmic sampling (-d
      ### option in Nonpareil), 0, 1, -1, 1.3, -1.1, 1.5, -1.5.
      skip.model=FALSE,
      ### If set, skips the model estimation altogether.
      ...
      ### Any additional parameters passed to `plot.Nonpareil.Curve`
      ){
  # Check parameters and initialize object
  if(is.na(label)){
    label <- basename(file)
    if(substr(label, nchar(label)-3, nchar(label))==".npo"){
      label <- substr(label, 0, nchar(label)-4)
    }
  }
  if(is.na(col)){
    col <- rgb(sample(200,1),sample(200,1),sample(200,1),maxColorValue=255)
  }
  np <- new("Nonpareil.Curve", file=as.character(file), label=label, col=col,
        star=star, has.model=FALSE, call=match.call(), diversity=0)

  # Read metadata (.npo headers)
  np <- Nonpareil.read_metadata(np)

  # Read data (.npo table)
  np <- Nonpareil.read_data(np, correction.factor)
  if(!np$consistent & enforce.consistency){
    for(w in np$warnings) warning(w)
    return(invisible(np))
  }

  # Fit model
  if(!skip.model) np <- Nonpareil.fit_model(np, weights.exp)

  # Plot
  if(plot) plot(np, ...)

  # Warnings
  for(w in np$warnings) warning(w)

  # Return
  invisible(np)
  ### Returns invisibly a `Nonpareil.Curve` object
}, ex=function(){
  # Generate a Nonpareil plot
  file <- system.file("extdata", "LakeLanier.npo", package="Nonpareil")
  np <- Nonpareil.curve(file)

  # Show the estimated values
  print(np)

  # Predict coverage for 20Gbp
  predict(np, 20e9)

  # Obtain the Nd diversity index
  np$diversity
})
Nonpareil.set <- structure(function(
      ### Generates a collection of Nonpareil curves (a `Nonpareil.Set` object)
      ### and (optionally) plots all of them in a single canvas
      files,
      ### Vector with the paths to the .npo files
      col=NA,
      ### Color of the curves (vector). If not passed, values are randomly
      ### assigned. Values are recycled
      labels=NA,
      ### Labels of the curves (vector). If not passed, values are determined by
      ### the filename. Values are recycled
      plot=TRUE,
      ### If TRUE, it generates the Nonpareil curve plots
      plot.opts=list(),
      ### Any parameters accepted by `plot.Nonpareil.Set` as a list
      ...
      ### Any additional parameters accepted by `Nonpareil.curve`
      ){
  files <- as.vector(files)
  y <- new("Nonpareil.Set", call=match.call())
  col <- rep(col, length.out=length(files))
  labels <- rep(labels, length.out=length(files))

  nonpareil.opts <- list(...)
  nonpareil.opts[["plot"]] <- FALSE
  for(i in 1:length(files)){
    nonpareil.opts[["file"]] <- files[i]
    nonpareil.opts[["col"]] <- col[i]
    nonpareil.opts[["label"]] <- labels[i]
    y$np.curves[[ length(y$np.curve)+1 ]] <- do.call("Nonpareil.curve", nonpareil.opts)
  }

  # Plot
  if(plot){
    plot.opts[["x"]] <- y
    y <- do.call("plot", plot.opts)
  }

  # Return
  invisible(y)
  ### Returns invisibly a `Nonpareil.Set` object
}, ex=function(){
  # Generate a Nonpareil plot with multiple curves
  files <- system.file("extdata",
        c("HumanGut.npo","LakeLanier.npo","IowaSoil.npo"), package="Nonpareil")
  col <- c("orange","darkcyan","firebrick4")
  nps <- Nonpareil.set(files, col=col,
        plot.opts=list(plot.observed=FALSE, model.lwd=2))

  # Show the estimated values
  print(nps)

  # Show current coverage (as %)
  summary(nps)[,"C"]*100

  # Extract Nd diversity index
  summary(nps)[,"diversity"]

  # Extract sequencing effort for nearly complete coverage (in Gbp)
  summary(nps)[,"LRstar"]/1e9

  # Predict coverage for a sequencing effort of 10Gbp
  sapply(nps$np.curves, predict, 10e9)
})
Nonpareil.curve.batch <- Nonpareil.set
