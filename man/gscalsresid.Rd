% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustergscairls.r
\name{gscalsresid}
\alias{gscalsresid}
\title{For use in clustergscairls, residuals of a gsc-model}
\usage{
gscalsresid(dat, out, indicatorx, indicatory, loadingx, loadingy)
}
\arguments{
\item{dat}{(n,p) data matrix}

\item{out}{list, output from gscals}

\item{indicatorx}{vector describing with which exogenous composite the X-variables are connected}

\item{indicatory}{vector describing with which endogenous composite the Y-variables are connected}

\item{loadingx}{logical TRUE when there are loadings for the X-variables in the model}

\item{loadingy}{logical TRUE when there are loadings for the y-variables in the model}
}
\value{
resid (n,q2) matrix of residuals from structural model, the q2 is the number of endogenous composites .
}
\description{
\code{gscalsresid} computes the residuals of a gsc-model when the parameters and weights are given
}
\examples{
data(simplemodel) 
data(gscalsout)
B <- matrix(c( 0,0,0, 0,0,0, 0.7,0.4,0),3,3,byrow=TRUE)
indicatorx <- c(1,1,1,2,2,2)
indicatory <- c(1,1,1)  
out <- gscalsresid(simplemodel,gscalsout,indicatorx,indicatory,TRUE,TRUE)

}
