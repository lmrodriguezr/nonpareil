\name{Nonpareil.Set-class}
\Rdversion{1.1}
\docType{class}
\alias{Nonpareil.Set-class}
\alias{Nonpareil.Set}
\alias{+,Nonpareil.Set}
%% \alias{$}
%% \alias{$<-}

\title{Nonpareil.Set S4 class}
\description{Collection of `Nonpareil.Curve` objects. This object can be produced by
`Nonpareil.curve.batch` and supports S4 methods `plot`, `summary`, and
`print`}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new(Nonpareil.Set ...)}}
\section{Slots}{
  \describe{
    \item{\code{np.curves}:}{(\code{list}) List of `Nonpareil.Curve` objects}
    \item{\code{call}:}{(\code{call}) Call producing this object}
  }
}
\section{Methods}{
  \describe{
    \item{+}{\code{signature(e1 = "Nonpareil.Set", e2 = "ANY")}: ... }
    \item{$}{\code{signature(x = "Nonpareil.Set")}: ... }
    \item{$<-}{\code{signature(x = "Nonpareil.Set")}: ... }
	 }
}

\author{Luis M. Rodriguez-R [aut, cre]}





