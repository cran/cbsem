% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CovGSCmodel_ff.r
\name{gscmcovff}
\alias{gscmcovff}
\title{\code{gscmcovff} determines the covariance matrix of a GSC model belonging to scenario ff.}
\usage{
gscmcovff(B, indicatorx, indicatory, wx, wy, Sxixi, R2 = NULL)
}
\arguments{
\item{B}{(q,q) lower triangular matrix describing the interrelations of the latent variables:
b_ij = 1 regression coefficient of eta_j in the regression relation in which eta_i is the depend variable
b_ij = 0 if eta_i does not depend on eta_j in a direct way  (b_ii = 0 !)}

\item{indicatorx}{vector describing with which exogenous composite the X-variables are connected}

\item{indicatory}{vector describing with which endogenous composite the Y-variables are connected}

\item{wx}{vector of weights for building exogenous composites or NULL when loadings are present}

\item{wy}{vector of weights for building endogenous composites or NULL when loadings are present}

\item{Sxixi}{covariance matrix of exogenous composites}

\item{R2}{vector of coefficients of determination for regressions belonging to the structural model}
}
\value{
out list with components
       \tabular{ll}{ S \tab covariance matrix of manifest variables   \cr
         B \tab (q,q) lower triangular matrix with possibly modified coefficients of the structural model \cr
         Scomp \tab covariance matrix of composites   \cr
         wx \tab  vector of weights for building exogenous composites \cr
         wy \tab  vector of weights for building endoogenous composites \cr
        }
}
\description{
\code{gscmcovff} determines the covariance matrix of a GSC model belonging to scenario ff.
}
\examples{
B <- matrix(c(0,0,0, 0,0,0, 0.7,0.4,0),3,3,byrow=TRUE) 
indicatorx <- c(1,1,1,2,2,2)
indicatory <- c(1,1,1)   
Sxixi <- matrix(c(1.0, 0.01, 0.01, 1),2,2)    
wx <- c(0.46, 0.31, 0.32, 0.34, 0.40, 0.37)
wy <- c(0.41, 0.39, 0.37)    
out <- gscmcovff(B,indicatorx,indicatory,wx,wy,Sxixi,R2=NULL)

}
