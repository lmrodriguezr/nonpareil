\name{Nonpareil.set}
\alias{Nonpareil.set}
\title{Nonpareil set}
\description{Generates a collection of Nonpareil curves (a `Nonpareil.Set` object)
and (optionally) plots all of them in a single canvas}
\usage{Nonpareil.set(files, col = NA, labels = NA, plot = TRUE, plot.opts = list(), 
    ...)}
\arguments{
  \item{files}{Vector with the paths to the .npo files}
  \item{col}{Color of the curves (vector). If not passed, values are randomly
assigned. Values are recycled}
  \item{labels}{Labels of the curves (vector). If not passed, values are determined by
the filename. Values are recycled}
  \item{plot}{If TRUE, it generates the Nonpareil curve plots}
  \item{plot.opts}{Any parameters accepted by `plot.Nonpareil.Set` as a list}
  \item{\dots}{Any additional parameters accepted by `Nonpareil.curve`}
}

\value{Returns invisibly a `Nonpareil.Set` object}

\author{Luis M. Rodriguez-R [aut, cre]}




\examples{
# Generate a Nonpareil plot with multiple curves
files <- system.file("extdata",
      c("HumanGut.npo","LakeLanier.npo","IowaSoil.npo"), package="Nonpareil")
col <- c("orange","darkcyan","firebrick4")
nps <- Nonpareil.set(files, col=col,
      plot.opts=list(plot.observed=FALSE, model.lwd=2))

# Show the estimated values
print(nps)

# Show current coverage (as \%)
summary(nps)[,"C"]*100

# Extract Nd diversity index
summary(nps)[,"diversity"]

# Extract sequencing effort for nearly complete coverage (in Gbp)
summary(nps)[,"LRstar"]/1e9

# Predict coverage for a sequencing effort of 10Gbp
sapply(nps$np.curves, predict, 10e9)
}
