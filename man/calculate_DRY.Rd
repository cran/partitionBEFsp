% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partitioning_functions.R
\name{calculate_DRY}
\alias{calculate_DRY}
\title{Calculate change in relative yield}
\usage{
calculate_DRY(P, M, Q = length(M))
}
\arguments{
\item{P}{biomass of species grown in polyculture}

\item{M}{biomass of species grow in monoculture - note, must include the same species as P, listed in the same order}

\item{Q}{number of species in the community -defaults to length(M), but note that if you are calculating DRY for a large community of Q species of which only N are observed, you should set Q=Q, rather than Q=N.}
}
\value{
a list of changes in relative yields
}
\description{
calculates change in relative yield, DRY, comparing observed relative yield to the expected yield 1/Q
}
\examples{
# Please see package help file (?partitionBEFsp) for examples.
}
