\name{Nonpareil.Curve-class}
\Rdversion{1.1}
\docType{class}
\alias{Nonpareil.Curve-class}
\alias{Nonpareil.Curve}
%% \alias{$}
%% \alias{$<-}

\title{Nonpareil.Curve S4 class}
\description{A single Nonpareil curve. This object can be produced by `Nonpareil.curve`
and supports S4 methods `plot`, `summary`, `print`, and `predict`.
For additional details, see help for `summary.Nonpareil.Curve`}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new(Nonpareil.Curve ...)}}
\section{Slots}{
  \describe{
    \item{\code{file}:}{(\code{character}) Input .npo file}
    \item{\code{label}:}{(\code{character}) Name of the dataset}
    \item{\code{col}:}{(\code{character}) Color of the dataset}
    \item{\code{L}:}{(\code{numeric}) Read length}
    \item{\code{AL}:}{(\code{numeric}) Adjusted read length (same as L for alignment)}
    \item{\code{R}:}{(\code{numeric}) Number of reads}
    \item{\code{LR}:}{(\code{numeric}) Effective sequencing effort used}
    \item{\code{overlap}:}{(\code{numeric}) Minimum read overlap}
    \item{\code{ksize}:}{(\code{numeric}) K-mer size (for kmer kernel only)}
    \item{\code{log.sample}:}{(\code{numeric}) Multiplier of the log-sampling (or zero if linear)}
    \item{\code{kernel}:}{(\code{character}) Read-comparison kernel}
    \item{\code{version}:}{(\code{character}) Nonpareil version used}
    \item{\code{x.obs}:}{(\code{numeric}) Rarefied sequencing effort}
    \item{\code{x.adj}:}{(\code{numeric}) Adjusted rarefied sequencing effort}
    \item{\code{y.red}:}{(\code{numeric}) Rarefied redundancy (observed)}
    \item{\code{y.cov}:}{(\code{numeric}) Rarefied coverage (corrected)}
    \item{\code{y.sd}:}{(\code{numeric}) Standard deviation of rarefied coverage}
    \item{\code{y.p25}:}{(\code{numeric}) Percentile 25 (1st quartile) of rarefied coverage}
    \item{\code{y.p50}:}{(\code{numeric}) Percentile 50 (median) of rarefied coverage}
    \item{\code{y.p75}:}{(\code{numeric}) Percentile 75 (3rd quartile) of rarefied coverage}
    \item{\code{kappa}:}{(\code{numeric}) Dataset redundancy}
    \item{\code{C}:}{(\code{numeric}) Dataset coverage}
    \item{\code{consistent}:}{(\code{logical}) Is the data sufficient for accurate estimation?}
    \item{\code{star}:}{(\code{numeric}) Coverage considered 'nearly complete'}
    \item{\code{has.model}:}{(\code{logical}) Was the model successfully estimated?}
    \item{\code{warning}:}{(\code{character}) Warnings generated on consistency or model fitting}
    \item{\code{LRstar}:}{(\code{numeric}) Projected seq. effort for nearly complete coverage}
    \item{\code{modelR}:}{(\code{numeric}) Pearson's R for the estimated model}
    \item{\code{diversity}:}{(\code{numeric}) Dataset Nd index of sequence diversity}
    \item{\code{model}:}{(\code{list}) Fitted sigmoidal model}
    \item{\code{call}:}{(\code{call}) Call producing this object}
  }
}
\section{Methods}{
  \describe{
    \item{$}{\code{signature(x = "Nonpareil.Curve")}: ... }
    \item{$<-}{\code{signature(x = "Nonpareil.Curve")}: ... }
	 }
}

\author{Luis M. Rodriguez-R [aut, cre]}





