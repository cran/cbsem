#' Output of gscals for the simplemodel data.
#'
#' A list containing the result of gscals for the simplemodel data.
#'
#' @format A list with entries:
#' \describe{
#'   \item{$Bhat}{estimated esign matrix of the simple model}
#'   \item{$What}{matrix of weights} 
#'   \item{$lambdahat}{mvector of estimated loadings}
#'   \item{$iter}{number of iterations} 
#'   \item{$fehl}{maximal difference of parameter estimates for the last and second last iteration} 
#'   \item{$composit}{data matrix of composites} 
#'   \item{$resid}{data matrix of residuals of the structural model}  
#'   \item{$S}{Covariance matrix of manifest variables}
#'   \item{$ziel}{sum of squared residuals for the final sum}
#'   \item{$fi}{The value of the fit criterion}
#'   \item{$R2}{vector with the coefficients of determination for structural regressions}
#' } 
"gscalsout"
