##
##  covgscmodel_fr.r                                                                  R. Schlittgen 2.12.2017
##
#' \code{gscmcovfr} determines the covariance matrix of a GSC model belonging to scenario fr. The covariance matrices of the errors
#' are supposed to be diagonal.
#' @param B   (q,q) lower triangular matrix describing the interrelations of the latent variables:
#'             b_ij = 1 regression coefficient of eta_j in the regression relation in which eta_i is the depend variable
#'             b_ij = 0 if eta_i does not depend on eta_j in a direct way  (b_ii = 0 !)
#' @param  indicatorx vector describing with which exogenous composite the X-variables are connected
#' @param  indicatory vector describing with which endogenous composite the Y-variables are connected
#' @param  lambday  vector of loadings of indicators for endogenous composites
#' @param  wx  vector of weights for building exogenous composites
#' @param  Sxixi covariance matrix of exogenous composites
#' @param  R2 vector of coefficients of determination for regressions belonging to the structural model
#' @return out list with components
#'        \tabular{ll}{ S \tab covariance matrix of manifest variables   \cr
#'          B \tab (q,q) lower triangular matrix with possibly modified coefficients of the structural model \cr
#'          Scomp \tab covariance matrix of composites   \cr
#'          wx \tab  vector of weights for building exogenous composites \cr
#'          See \tab  diagonal matrix of variances of errors of Y variable loadings or NA
#'         }
#' @examples
#' Sxixi <- matrix(c(1.0,  0.01,  0.01, 1),2,2)   
#' B <- matrix(c( 0,0,0, 0,0,0, 0.7,0.4,0),3,3,byrow=TRUE) 
#' indicatorx <- c(1,1,1,2,2,2)
#' indicatory <- c(1,1,1)  
#' lambday <- c(0.89,0.90,0.80)    
#' wx <- c(0.46, 0.31, 0.32, 0.34, 0.40, 0.37) 
#' out <- gscmcovfr(B,indicatorx,indicatory,lambday,wx,Sxixi,R2=NULL)      
#'
#' @export
gscmcovfr <- function(B,indicatorx,indicatory,lambday,wx,Sxixi,R2=NULL){

q <- nrow(B)                                         # number of latent variables
b <- rep(1,q)
equals0 <- function(x) all(x == 0)
q1 <- length(b[apply(B,1,equals0)])         # number of exogenous composites
q2 <- q - q1                                          # number of endogenous composites
p1 <- length(indicatorx)                         # number of indicators of exogenous composites
p2 <- length(indicatory)                         # number of indicators of endogenous composites
p <- p1+p2                                           # number of manifest variables


Sdd <- NA
See <- NA
                                                               # construction of Sxx, W1
   Sxx <- matrix(rep(c(p1:0),p1),p1,p1)*0.01 
   Sxx[upper.tri(Sxx)] <- 0
   Sxx <- (Sxx+t(Sxx))/p1
   diag(Sxx) <- 1

   W1 <- matrix(indicatorx,p1,q1)                # construction of W1
   W1 <- 1*(W1 == matrix(c(1:q1),p1,q1,byrow=T))
   P <- cumsum(colSums(W1))
   P <- c(0,P)
   P1 <- P

   W1 <- W1*wx
   W1 <- scale(W1,center=FALSE, scale= sqrt(diag(t(W1)%*%Sxx%*%W1)) )                                    # scaling of weights 
   if (q1 > 1){
   for (j in c(2:q1)){                                     # scaling of off-diagonal of Sxx
       for (i in c(1:(j-1))){
            subSxx <-  Sxx[(P[i]+1):P[i+1],(P[j]+1):P[j+1]]
            f <-  t(W1[(P[i] + 1):P[i+1],i,drop=F])%*%subSxx%*%W1[(P[j] + 1):P[j+1],j,drop=F]
            f <- as.numeric(Sxixi[i,j])/as.numeric(f)
            Sxx[(P[i]+1):P[i+1],(P[j]+1):P[j+1]] <- subSxx*f
            Sxx[(P[j]+1):P[j+1],(P[i]+1):P[i+1]] <- t(subSxx*f)
            }
    }}
   W1 <- scale(W1,center=FALSE, scale= sqrt(diag(t(W1)%*%Sxx%*%W1)) )                                    # scaling of weights 
   
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
   IB <- solve(diag(q2) - t(Br2) )

# construction of Syy, Sxy

   Ly <- matrix(indicatory,p2,q2)                     # construction of Ly
   Ly <- 1*(Ly == matrix(c(1:q2),p2,q2,byrow=T))
   Ly <- Ly*lambday

   Syy <- Ly%*%Setaeta%*%t(Ly)
   See <- diag(1-diag(Syy))
   diag(Syy) <-1

   Sxy <- Sxx%*%W1%*%t(Br1)%*%IB%*%t(Ly)
   S <- rbind(cbind(Sxx,Sxy),cbind(t(Sxy),Syy))

out <- list(S=S,B=B,Scomp=Sll,wx=colSums(W1),See=See)
return(out)
}
