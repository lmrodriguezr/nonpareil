\name{predict.Nonpareil.Curve}
\alias{predict.Nonpareil.Curve}
\title{predict Nonpareil Curve}
\description{Predict the coverage for a given sequencing effort}
\usage{\method{predict}{Nonpareil.Curve}(object, lr = object$LR, ...)}
\arguments{
  \item{object}{`Nonpareil.Curve` object}
  \item{lr}{Sequencing effort for the prediction (in bp)}
  \item{\dots}{Additional parameters ignored}
}

\value{Returns the expected coverage at the given sequencing effort.}

\author{Luis M. Rodriguez-R [aut, cre]}





