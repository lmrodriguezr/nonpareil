\name{plot.Nonpareil.Set}
\alias{plot.Nonpareil.Set}
\title{plot Nonpareil Set}
\description{Plot a `Nonpareil.Set` object}
\usage{\method{plot}{Nonpareil.Set}(x, col = NA, labels = NA, main = "Nonpareil Curves", 
    legend.opts = list(), ...)}
\arguments{
  \item{x}{`Nonpareil.Set` object to plot}
  \item{col}{Color of the curves (vector). If passed, it overrides the colors set
in the `Nonpareil.Curve` objects. Values are recycled}
  \item{labels}{Labels of the curves (vector). If passed, it overrides the labels set
in the `Nonpareil.Curve` objects. Values are recycled}
  \item{main}{Title of the plot}
  \item{legend.opts}{Any additional parameters passed to `Nonpareil.legend`. If FALSE, the
legend is not displayed}
  \item{\dots}{Any additional parameters passed to `plot.Nonpareil.Curve`}
}

\value{Returns invisibly a `Nonpareil.Set` object (same as `x` input).}

\author{Luis M. Rodriguez-R [aut, cre]}





