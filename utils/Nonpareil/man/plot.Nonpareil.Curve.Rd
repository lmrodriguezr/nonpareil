\name{plot.Nonpareil.Curve}
\alias{plot.Nonpareil.Curve}
\title{plot Nonpareil Curve}
\description{Plot a `Nonpareil.Curve` object}
\usage{\method{plot}{Nonpareil.Curve}(x, col = NA, add = FALSE, new = !add, plot.observed = TRUE, 
    plot.model = TRUE, plot.dispersion = FALSE, plot.diversity = TRUE, 
    xlim = c(1000, 1e+13), ylim = c(1e-06, 1), main = paste("Nonpareil Curve for", 
        x$label), xlab = "Sequencing effort (bp)", ylab = "Estimated Average Coverage", 
    curve.lwd = 2, curve.alpha = 0.4, model.lwd = 1, model.alpha = 1, 
    log = "x", arrow.length = 0.05, arrow.head = arrow.length, 
    ...)}
\arguments{
  \item{x}{`Nonpareil.Curve` object to plot}
  \item{col}{Color of the curve. If passed, it overrides the colors set in the
`Nonpareil.Curve` object}
  \item{add}{If TRUE, it attempts to use a previous (active) canvas to plot the
curve}
  \item{new}{Inverse of `add`}
  \item{plot.observed}{Indicates if the observed (rarefied) coverage is to be plotted}
  \item{plot.model}{Indicates if the fitted model is to be plotted}
  \item{plot.dispersion}{Indicates if (and how) dispersion of the replicates should be plotted.
Supported values are:
FALSE: no dispersion is plotted (default),
'sd': one standard deviation around the mean,
'ci95': 95\% confidence interval,
'ci90': 90\% confidence interval,
'ci50': 50\% confidence interval,
'iq': Inter-quartile range}
  \item{plot.diversity}{If TRUE, the diversity estimate is plotted as a small arrow below the
Nonpareil curve}
  \item{xlim}{Limits of the sequencing effort (X-axis)}
  \item{ylim}{Limits of the coverage (Y-axis)}
  \item{main}{Title of the plot}
  \item{xlab}{Label of the X-axis}
  \item{ylab}{Label of the Y-axis}
  \item{curve.lwd}{Line width of the rarefied coverage}
  \item{curve.alpha}{Alpha value (from 0 to 1) of the rarefied coverage}
  \item{model.lwd}{Line width of the model}
  \item{model.alpha}{Alpha value (from 0 to 1) of the model}
  \item{log}{Axis to plot in logarithmic scale. Supported values are:
'x': sequencing effort (default),
'y': coverage,
'xy': both logarithmic, or
'': both linear}
  \item{arrow.length}{If plot.diversity=TRUE, it determines the length of the arrow to
display the divsersity (as a fraction of the ylim range).}
  \item{arrow.head}{If plot.diversity=TRUE, it determines the length of the arrow head to
display the diversity index (in inches).}
  \item{\dots}{Additional graphical parameters}
}

\value{Retuns invisibly a `Nonpareil.Curve` object (same as `x` input). For
additional details see help for `summary.Nonpareil.Curve`}

\author{Luis M. Rodriguez-R [aut, cre]}





