\name{Nonpareil.curve.batch}
\alias{Nonpareil.curve.batch}
\title{Nonpareil curve batch}
\description{Generates a collection of Nonpareil curves (a `Nonpareil.Set` object)
and (optionally) plots all of them in a single canvas}
\usage{Nonpareil.curve.batch(files, col = NA, labels = NA, plot = TRUE, 
    plot.opts = list(), ...)}
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

\author{Luis M. Rodriguez-R <lmrodriguezr@gmail.com>}





