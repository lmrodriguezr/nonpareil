\name{summary.Nonpareil.Curve}
\alias{summary.Nonpareil.Curve}
\title{summary Nonpareil Curve}
\description{Returns a summary of the Nonpareil.Curve results}
\usage{\method{summary}{Nonpareil.Curve}(object, ...)}
\arguments{
  \item{object}{`Nonpareil.Curve` object}
  \item{\dots}{Additional parameters ignored}
}

\value{Returns a matrix with the following values for the dataset:

kappa: "Redundancy" value of the entire dataset.

C: Average coverage of the entire dataset.

LRstar: Estimated sequencing effort required to reach the objective
average coverage (star, 95\% by default).

LR: Actual sequencing effort of the dataset.

modelR: Pearson's R coefficient betweeen the rarefied data and the
projected model.

diversity: Nonpareil sequence-diversity index (Nd). This value's units are
the natural logarithm of the units of sequencing effort (log-bp), and
indicates the inflection point of the fitted model for the Nonpareil
curve. If the fit doesn't converge, or the model is not estimated, the
value is zero (0).
}

\author{Luis M. Rodriguez-R [aut, cre]}





