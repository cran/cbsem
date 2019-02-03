
##
##  covgscmodel.r                                                                  R. Schlittgen 2.12.2017
##
#' Determination of the covariance matrix of a GSC model belonging to scenario 1, scenario 2, scenario 3
#'
#' \code{gscmcov} determines the covariance matrix of a GSC model. This is a wrapper for the functions
#' gscmcovrr, gscmcovfr and gscmcovff
#' @param B   (q,q) lower triangular matrix describing the interrelations of the latent variables:
#'             b_ij = 1 regression coefficient of eta_j in the regression relation in which eta_i is the depend variable
#'             b_ij = 0 if eta_i does not depend on eta_j in a direct way  (b_ii = 0 !)
#' @param  indicatorx vector describing with which exogenous composite the X-variables are connected
#' @param  indicatory vector describing with which endogenous composite the Y-variables are connected
#' @param  lambdax  vector of loadings of indicators for exogenous composites or NULL when there are no loadings for the X-variables in the model
#' @param  lambday  vector of loadings of indicators for endogenous composites or NULL when there are no loadings for the Y-variables in the model
#' @param  wx  vector of weights for building exogenous composites or NULL when loadings are present
#' @param  wy  vector of weights for building endogenous composites or NULL when loadings are present
#' @param  Sxixi covariance matrix of exogenous composites
#' @param  R2 vector of coefficients of determination for regressions belonging to the structural model
#' @return out list with components
#'        \tabular{ll}{ S \tab covariance matrix of manifest variables   \cr
#'          B \tab (q,q) lower triangular matrix with possibly modified coefficients of the structural model \cr
#'          Scomp \tab covariance matrix of composites   \cr
#'          wx \tab  vector of weights for building exogenous composites \cr
#'          wy \tab  vector of weights for building endoogenous composites \cr
#'          Sdd \tab diagonal matrix of variances of errors of X variable loadings or NA  \cr
#'          See \tab  diagonal matrix of variances of errors of Y variable loadings or NA
#'         }
#' @examples
#' Sxixi <- matrix(c(1.0,  0.01,  0.01, 1),2,2)   
#' B <- matrix(c( 0,0,0, 0,0,0, 0.7,0.4,0),3,3,byrow=TRUE) 
#' indicatorx <- c(1,1,1,2,2,2)
#' indicatory <- c(1,1,1)                     
#' lambdax <- c(0.83,0.87,0.87,0.91,0.88,0.82)
#' lambday <- c(0.89,0.90,0.80)  
#' wx <- c(0.46, 0.31, 0.32, 0.34, 0.40, 0.37)
#' wy <- c(0.41, 0.39, 0.37)    
#' out <- gscmcov(B,indicatorx,indicatory,lambdax,lambday,wx=NULL,wy=NULL,Sxixi,R2=NULL)        
#'
#' @export
gscmcov <- function(B,indicatorx,indicatory,lambdax=NULL,lambday=NULL,wx=NULL,wy=NULL,Sxixi,R2=NULL){

Sdd <- NA
See <- NA

if( ( is.numeric(lambdax) ) & ( is.numeric(lambday) )){        

   out <- gscmcovrr(B,indicatorx,indicatory,lambdax=lambdax,lambday=lambday,Sxixi,R2=R2)
   S <- out$S
   Scomp <- out$Scomp
   B <- out$B
   Sdd <- out$Sdd
   See <- out$See

  }else  if (( !is.numeric(lambdax) ) & ( length(wx) > 0 ) & (is.numeric(lambday) ) ){

  out <- gscmcovfr(B,indicatorx,indicatory,lambday,wx,Sxixi,R2=R2)
  S <- out$S
  Scomp <- out$Scomp
  B <- out$B
  Sdd <- out$Sdd
  See <- out$See
  wx <- out$wx

}else{

  out <- gscmcovff(B,indicatorx,indicatory,wx,wy,Sxixi,R2=R2)
  S <- out$S
  Scomp <- out$Scomp
  B <- out$B
  wx <- out$wx
  wy <- out$wy
}

 out <- list(S=S,B=B,Scomp=Scomp,wx=wx,wy=wy,Sdd=Sdd,See=See)
 return(out)
}
