% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createlms.R
\name{one_iteration}
\alias{one_iteration}
\title{one iteration}
\usage{
one_iteration(
  data.list,
  method,
  prop.fam = 0.75,
  prop.subject = 1,
  age.min = 0,
  age.max = 18,
  age.int = 1/12,
  keep.models = F,
  dist = "BCCGo",
  formula = NULL,
  sigma.df = 3,
  nu.df = 2,
  mu.df = 4,
  tau.df = 2,
  sigma.formula = ~1,
  nu.formula = ~1,
  tau.formula = ~1,
  verbose = F,
  trans.x = F,
  lim.trans = c(0, 1.5),
  method.pb = "ML"
)
}
\arguments{
\item{data.list}{list of dataframes as returned by prepare_data}

\item{method}{use vgam or gamlss}

\item{prop.fam}{proportion of families to be sampled}

\item{prop.subject}{proportion of subject to be sampled}

\item{age.min}{lower bound of age}

\item{age.max}{upper bound of age}

\item{age.int}{stepwidth of the age variable}

\item{keep.models}{indicator whether or not models in each iteration should be kept}

\item{dist}{distribution used for the fitting process, has to be one of BCCGo, BCPEo, BCTo as they are accepted by lms()}

\item{formula}{formula for the location parameter}

\item{sigma.df}{degree of freedem spread parameter}

\item{nu.df}{degree of freedem skewness parameter}

\item{mu.df}{degree of freedem location parameter}

\item{tau.df}{degree of freedem kurtosis parameter}

\item{sigma.formula}{formula for the sigma parameter}

\item{nu.formula}{formula for the nu parameter}

\item{tau.formula}{formula for the tau parameter}

\item{verbose}{whether or not information about sampling will be printed during while iterate}

\item{trans.x}{indicator wether age should be transformed or not}

\item{lim.trans}{limits for the exponent of transformation of age}

\item{method.pb}{GAIC or ML}
}
\value{
list of lists each containing a dataframe of the fitted lms parameter at the given age points and the fitted model
}
\description{
one iteration
}
\details{
function samples families then measurements and fits the model
the function is called inside \code{\link{do_iterations}} and may not called directly
}
\author{
Mandy Vogel
}
