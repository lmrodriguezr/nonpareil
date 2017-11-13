\name{Nonpareil.curve}
\alias{Nonpareil.curve}
\title{Nonpareil curve}
\description{Generates a Nonpareil curve from a .npo file.}
\usage{Nonpareil.curve(file, ksize = NULL, overlap = NULL, factor = 1, 
    plotDispersion = NA, returnModelValues = FALSE, returnModelParameters = FALSE, 
    xmax = 1e+13, xmin = 1000, ymax = 1, ymin = 1e-06, xlab = NULL, 
    ylab = NULL, col = NA, r = NA, g = NA, b = NA, new = TRUE, 
    plot = TRUE, libname = NA, modelOnly = FALSE, plotModel = TRUE, 
    plotDiversity = FALSE, curve.lwd = 2, curve.alpha = 0.4, 
    model.lwd = 1, model.alpha = 1, log = "x", data.consistency = TRUE, 
    useValue = "mean", star = 95, read.length = NA, weights.exp = NA, 
    correction.factor = TRUE, ...)}
\arguments{
  \item{file}{Path to the .npo file, containing the read redundancy.}
  \item{ksize}{K-mer size (if -T kmer). By default this value is extracted
from the .npo file headers.}
  \item{overlap}{Value of the '-L' parameter (in Nonpareil, the default is 50). If not
set, it tries to find the value in the .npo file (supported in Nonpareil
>= 2.0), or fails with an error message. Use only with Nonpareil < v2.0.}
  \item{factor}{Multiplier of the sequencing effort. This can be used to express the
sequencing effort in units other than base pairs (bp). For example, to
express sequencing effort as Gbp, use factor=1e-9. This can also affect
the fit of the model, and it's considered EXPERIMENTAL.}
  \item{plotDispersion}{Indicates if (and how) dispersion of the replicates should be plotted.
It requires modelOnly=FALSE to take effect. It can be NA, in which case
no dispersion is plotted, or any of the following strings: 'sd' (one
standard deviation around the mean), 'ci95' (95\% confidence interval),
'ci90' (90\% confidence interval), 'ci50' (50\% confidence interval), 'iq'
(inter-quartile range).}
  \item{returnModelValues}{If TRUE, returns the coordinates of the model as model.x and model.y.}
  \item{returnModelParameters}{If TRUE, returns the model itself as model.}
  \item{xmax}{Maximum sequencing effort to plot.}
  \item{xmin}{Minimum sequencing effort to plot}
  \item{ymax}{Maximum coverage to plot.}
  \item{ymin}{Minimum coverage to plot.}
  \item{xlab}{Label of the X-axis. If NULL, it's set to sequencing effort and the
units (see factor).}
  \item{ylab}{Label of the Y-axis. If NULL, it's set to Estimated average coverage.}
  \item{col}{The color of the curve. If passed, it overrides `r`, `g`, and `b`.}
  \item{r}{Red component of the curve's color. If NA, it's randomly set. If <=1,
it's assumed to be in the range [0,1]; if >1, it's assumed to be in the
range [0,256].}
  \item{g}{Green component of the curve's color. Same as `r`.}
  \item{b}{Blue component of the curve's color. Same as `r`.}
  \item{new}{If FALSE, it attempts to use a previous (active) canvas to plot the
curve.}
  \item{plot}{Determines if the plot should be produced. If FALSE, it still computes
the coverage and the model.}
  \item{libname}{Name of the library. If NA, it's determined by the file name. This is
useful if you plan to call Nonpareil.legend().}
  \item{modelOnly}{If TRUE, the rarefied data is not presented, only the fitted model.}
  \item{plotModel}{If FALSE, the model is not plotted (but it's still computed).}
  \item{plotDiversity}{If TRUE, the diversity estimate is plotted as a small arrow below the
Nonpareil curve.}
  \item{curve.lwd}{Line width of the Nonpareil curve.}
  \item{curve.alpha}{Alpha value (from 0 to 1) of the Nonpareil curve.}
  \item{model.lwd}{Line width of the model.}
  \item{model.alpha}{Alpha value (from 0 to 1) of the model.}
  \item{log}{Axis to plot in logarithmic scale. It can be 'x' (sequencing effort,
default), 'y' (coverage), 'xy' (both logarithmic), or '' (both linear).}
  \item{data.consistency}{If TRUE, it checks the consistency of the data before plotting.}
  \item{useValue}{Controls how the replicates are to be summarized at each point of
sequencing effort. It can be any of: 'mean' (average of the replicates),
'median' (median), 'ub' (upper boundary of the 95\% confidence interval),
'lb' (lower boundary of the 95\% confidence interval), 'q1' (quartile 1),
'q3' (quartile 3). Note that the quartile 2 is also 'median'.}
  \item{star}{Objective coverage (in percentage). By default: 95, which means the
sequencing effort required to reach 95\% average coverage is to be
estimated.}
  \item{read.length}{Length of the reads. Use only with Nonpareil < v2.0.}
  \item{weights.exp}{Vector of values to be tested (in order) as exponent of the weights
distribution. If the model fails to converge, sometimes manual
modifications in this parameter may help. By default (NA), five
different values are tested in the following order: For linear sampling,
-1.1, -1.2, -0.9, -1.3, -1. For logarithmic sampling (-d option in
Nonpareil), 0, 1, -1, 1.3, -1.1, 1.5, -1.5.}
  \item{correction.factor}{Should the overlap-dependent correction factor be applied? If FALSE,
redundancy is assumed to equal coverage.}
  \item{\dots}{Any other parameters accepted by plot().}
}

\value{A list with the following entries:

kappa: "Redundancy" value of the entire dataset.

C: Average coverage of the entire dataset.

LRstar: Estimated sequencing effort required to reach the objective
average coverage (star, 95\% by default).

LR: Actual sequencing effort of the dataset.

modelR: Pearson's R coefficient betweeen the rarefied data and the
projected model.

diversity: Nonpareil-based diversity index. This value's units are the
natural logarithm of the units of sequencing effort (by default,
log-bp), and indicates the inflection point of the fitted model for the
Nonpareil curve. If the fit doesn't converge, or the model is not
estimated, the value is zero (0).

model.x (if retturnModelValues=TRUE): Values of sequencing effort in the
projected model.

model.y (if retturnModelValues=TRUE): Values of average coverage in the
projected model.

model (if returnModelParameters=TRUE): Fitted model.
}

\author{Luis M. Rodriguez-R <lmrodriguezr@gmail.com>}





