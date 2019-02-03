#' Output of covgscmodel for the simplemodel data.
#'
#' A list containing the result of gscmcov for the simplemodel data.
#'
#' @format A list with entries:
#' \describe{
#'   \item{$S}{Covariance matrix of manifest variables} 
#'   \item{$B}{Design matrix of the simple model}
#'   \item{$Scomp}{Covariance matrix of composites}
#'   \item{$wx}{weighting vector for exogenous composites}
#'   \item{$wy}{weighting vector for endogenous composites}
#'   \item{$Sdd}{diagonal covariance matrix of errors for loadings of X-variables}
#'   \item{$See}{diagonal covariance matrix of errors for loadings of Y-variables}
#' } 
"gscmcovout"


 