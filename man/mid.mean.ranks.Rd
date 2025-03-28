% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ridits.R
\name{mid.mean.ranks}
\alias{mid.mean.ranks}
\title{calculate mid mean ranks}
\usage{
mid.mean.ranks(
  v1,
  v2,
  maxrank = max(c(v2, v1)),
  minrank = NULL,
  conf.level = 0.05
)
}
\arguments{
\item{v1}{vector of ranks}

\item{v2}{vector of ranks}

\item{maxrank}{max possible rank}

\item{minrank}{min possible rank}

\item{conf.level}{confidence level alpha}
}
\value{
list of marginal mean ranks and mid ranks
}
\description{
calculate mid mean ranks
}
\author{
Mandy Vogel
}
