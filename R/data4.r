#' Simulated data.
#'
#' The data were simulated with two gsc models, both with two
#' exogeneous and one endogeneous composites. The exogeneous and
#' endegeneous composites have three indicators. There are no loadings.
#' The first 50 observations were simulated with one set of path
#' coefficients, the second 50 observations with another set.
#' the last column is the membership of a former clustering (k=2). 
#' 
#' @format A data frame with 10 variables and 50 cases:
#' \describe{
#'   \item{X1,X2,X3}{Indicators of first exogeneous composite} 
#'   \item{X4,X5,X6}{Indicators of second exogeneous composite}
#'   \item{Y1,Y2,Y3}{Indicators of endogeneous composite}
#'   \item{member}{membership of a former clustering}   
#' }  
#'
#' @export
"twoclm"