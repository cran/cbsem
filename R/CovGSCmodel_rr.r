##
##  covgscmodel_rr.r                                                                  R. Schlittgen 2.12.2017
##
#' \code{gscmcovrr} determines the covariance matrix of a GSC model belonging to scenario rr.
#' @param B (q,q) lower triangular matrix describing the interrelations of the latent variables:
#'           b_ij = 1 regression coefficient of eta_j in the regression relation in which eta_i is the depend variable
#'           b_ij = 0 if eta_i does not depend on eta_j in a direct way  (b_ii = 0 !)
#' @param  indicatorx vector describing with which exogenous composite the X-variables are connected
#' @param  indicatory vector describing with which endogenous composite the Y-variables are connected
#' @param  lambdax vector of loadings of indicators for exogenous composites
#' @param  lambday  vector of loadings of indicators for endogenous composites
#' @param  Sxixi covariance matrix of exogenous composites
#' @param  R2 vector of coefficients of determination for regressions belonging to the structural model
#' @return out list with components
#'        \tabular{ll}{ S \tab covariance matrix of manifest variables   \cr
#'          B \tab (q,q) lower triangular matrix with possibly modified coefficients of the structural model \cr
#'          Scomp \tab covariance matrix of composites   \cr
#'          Sdd \tab diagonal matrix of variances of errors of X variable loadings    \cr
#'          See \tab  diagonal matrix of variances of errors of Y variable loadings
#'         }
#' @examples
#' Sxixi <- matrix(c(1.0,  0.01,  0.01, 1),2,2)   
#' B <- matrix(c( 0,0,0, 0,0,0, 0.7,0.4,0),3,3,byrow=TRUE) 
#' indicatorx <- c(1,1,1,2,2,2)
#' indicatory <- c(1,1,1)                     
#' lambdax <- c(0.83,0.87,0.87,0.91,0.88,0.82)
#' lambday <- c(0.89,0.90,0.80)   
#' out <- gscmcovrr(B,indicatorx,indicatory,lambdax,lambday,Sxixi,R2=NULL)      
#'
#' @export
gscmcovrr <- function(B,indicatorx,indicatory,lambdax,lambday,Sxixi,R2=NULL){

q <- nrow(B)                                     # number of latent variables
b <- rep(1,q)
equals0 <- function(x) all(x == 0)
q1 <- length(b[apply(B,1,equals0)])    # number of exogenous composites
q2 <- q - q1                                     # number of endogenous composites
p1 <- length(indicatorx)                      # number of indicators of exogenous composites
p2 <- length(indicatory)                       # number of indicators of endogenous composites
p <- p1+p2                                # number of manifest variables

B1 <- B[(q1+1):q,1:q1,drop=FALSE]
B2 <- B[(q1+1):q,(q1+1):q,drop=FALSE]

                                                               # construction of Lx
   Lx <- matrix(indicatorx,p1,q1)
   Lx <- 1*(Lx == matrix(c(1:q1),p1,q1,byrow=T))
   Lx <- Lx*lambdax
                                                                # construction of Ly
   Ly <- matrix(indicatory,p2,q2)
   Ly <- 1*(Ly == matrix(c(1:q2),p2,q2,byrow=T))
   Ly <- Ly*lambday

   Sxx <- Lx%*%Sxixi%*%t(Lx)
   Sdd  <-  diag(1 - diag(Sxx))
   diag(Sxx) <- 1

#  computing covariance matrix of composites an scaling of path coefficients

   Sll <- diag(q)
   Sll[1:q1,1:q1] <- Sxixi

    if (length(R2) > 0){                                      # scaling of path coefficients
       for (m in c((q1+1):(q1+q2))){
         tau <- sqrt(R2[m-q1]/(B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),1:(m-1)]%*%t(B[m,1:(m-1),drop=F])))
         B[m,] <- tau*B[m,]
        }}
       for (m in c((q1+1):(q1+q2))){
        for (j in c(1:(m-1))){
          Sll[j,m] <- Sll[m,j] <- B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),j,drop=F]
       } }
   Setaeta <- Sll[(q1+1):(q1+q2),(q1+1):(q1+q2)]

   Br1 <- B[(q1+1):q,1:q1,drop=FALSE]
   Br2 <- B[(q1+1):q,(q1+1):q,drop=FALSE]
   IB <- solve(diag(q2) - t(Br2))

   Syy <- Ly%*%Setaeta%*%t(Ly)
   See <- diag(1-diag(Syy))
   diag(Syy) <-1

   Sxy <- Lx%*%Sxixi%*%t(Br1)%*%IB%*%t(Ly)
   S <- rbind(cbind(Sxx,Sxy),cbind(t(Sxy),Syy))

out <- list(S=S,B=B,Scomp=Sll,Sdd=Sdd,See=See)
return(out)
}
