\name{Nonpareil.curve}
\alias{Nonpareil.curve}
\title{Nonpareil curve}
\description{Generates a Nonpareil curve from an .npo file}
\usage{Nonpareil.curve(file, plot = TRUE, label = NA, col = NA, enforce.consistency = TRUE, 
    star = 95, correction.factor = TRUE, weights.exp = NA, skip.model = FALSE, 
    ...)}
\arguments{
  \item{file}{Path to the .npo file, containing the read redundancy}
  \item{plot}{Determines if the plot should be produced. If FALSE, it still computes
the coverage and the model}
  \item{label}{Name of the dataset. If NA, it is determined by the file name}
  \item{col}{Color of the curve. If NA, a random color is assigned (even if
plot=FALSE),}
  \item{enforce.consistency}{If TRUE, it fails verbosely on insufficient data, otherwise it warns
about the inconsistencies and attempts the estimations}
  \item{star}{Objective coverage in percentage; i.e., coverage value considered
near-complete}
  \item{correction.factor}{Should the overlap-dependent (or kmer-length-dependent) correction
factor be applied? If FALSE, redundancy is assumed to equal coverage.}
  \item{weights.exp}{Vector of values to be tested (in order) as exponent of the weights
distribution. If the model fails to converge, sometimes manual
modifications in this parameter may help. By default (NA), five
different values are tested in the following order: For linear
sampling, -1.1, -1.2, -0.9, -1.3, -1. For logarithmic sampling (-d
option in Nonpareil), 0, 1, -1, 1.3, -1.1, 1.5, -1.5.}
  \item{skip.model}{If set, skips the model estimation altogether.}
  \item{\dots}{Any additional parameters passed to `plot.Nonpareil.Curve`}
}

\value{Returns invisibly a `Nonpareil.Curve` object}

\author{Luis M. Rodriguez-R [aut, cre]}




\examples{
# Generate a Nonpareil plot
file <- system.file("extdata", "LakeLanier.npo", package="Nonpareil")
np <- Nonpareil.curve(file)

# Show the estimated values
print(np)

# Predict coverage for 20Gbp
predict(np, 20e9)

# Obtain the Nd diversity index
np$diversity
}
