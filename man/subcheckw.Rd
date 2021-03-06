% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Checkw.r
\name{subcheckw}
\alias{subcheckw}
\title{Function for use in checkw}
\usage{
subcheckw(w, indicator, S, L, Scomp)
}
\arguments{
\item{w}{vector of weights}

\item{indicator}{vector describing with which exogenous composite the indicators are connected}

\item{S}{covariance matrix of errors resulling from regession for loadings}

\item{L}{matrix of loadings}

\item{Scomp}{covariance matrix of composites}
}
\value{
out  scalar, sum of squared differences
}
\description{
\code{subcheckw} computes the sum of squared differences of two formulas for the covariancematrix of composites
}
