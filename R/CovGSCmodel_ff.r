##
##  covgscmodel_ff.r                                                                  R. Schlittgen 2.12.2017
## 
#' \code{gscmcovff} determines the covariance matrix of a GSC model belonging to scenario ff.
#' @param B   (q,q) lower triangular matrix describing the interrelations of the latent variables:
#'             b_ij = 1 regression coefficient of eta_j in the regression relation in which eta_i is the depend variable
#'             b_ij = 0 if eta_i does not depend on eta_j in a direct way  (b_ii = 0 !)
#' @param  indicatorx vector describing with which exogenous composite the X-variables are connected
#' @param  indicatory vector describing with which endogenous composite the Y-variables are connected
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
#'         }
#' @examples
#' B <- matrix(c(0,0,0, 0,0,0, 0.7,0.4,0),3,3,byrow=TRUE) 
#' indicatorx <- c(1,1,1,2,2,2)
#' indicatory <- c(1,1,1)   
#' Sxixi <- matrix(c(1.0, 0.01, 0.01, 1),2,2)    
#' wx <- c(0.46, 0.31, 0.32, 0.34, 0.40, 0.37)
#' wy <- c(0.41, 0.39, 0.37)    
#' out <- gscmcovff(B,indicatorx,indicatory,wx,wy,Sxixi,R2=NULL)
#'
#' @export
gscmcovff <- function(B,indicatorx,indicatory,wx,wy,Sxixi,R2=NULL){

q <- nrow(B)                                         # number of latent variables
b <- rep(1,q)
equals0 <- function(x) all(x == 0)
q1 <- length(b[apply(B,1,equals0)])         # number of exogenous composites
q2 <- q - q1                                          # number of endogenous composites
p1 <- length(indicatorx)                         # number of indicators of exogenous composites
p2 <- length(indicatory)                         # number of indicators of endogenous composites
p <- p1+p2                                           # number of manifest variables

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
   for (i in c(1:q1)){                                      # scaling of weights
     subSxx <- Sxx[(P[i] + 1):P[i+1],(P[i] + 1):P[i+1]]
     f <-  t(W1[(P[i] + 1):P[i+1],i,drop=F])%*%subSxx%*%W1[(P[i] + 1):P[i+1],i,drop=F]
     W1[(P[i] + 1):P[i+1],i] <- W1[(P[i] + 1):P[i+1],i]/sqrt(as.numeric(f))
     }
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

# construction of Syy, W2

   Syy <- matrix(rep(c(p2:0),p2),p2,p2)*0.01
   Syy[upper.tri(Syy)] <- 0
   Syy <- (Syy+t(Syy))/p2
   diag(Syy) <- 1


   W2 <- matrix(indicatory,p2,q2)                # construction of W2
   W2 <- 1*(W2 == matrix(c(1:q2),p2,q2,byrow=T))
   P <- cumsum(colSums(W2))
   P <- c(0,P)
   P2 <- P
   W2 <- W2*wy
   for (i in c(1:q2)){                                      # scaling of weights
     subSyy <- Syy[(P[i] + 1):P[i+1],(P[i] + 1):P[i+1]]
     f <-  t(W2[(P[i] + 1):P[i+1],i,drop=F])%*%subSyy%*%W2[(P[i] + 1):P[i+1],i,drop=F]
     W2[(P[i] + 1):P[i+1],i] <- W2[(P[i] + 1):P[i+1],i]/sqrt(as.numeric(f))
     }
 
   if (q2 > 1){
   for (j in c(2:q2)){                                     # scaling of off-diagonal of Syy
       for (i in c(1:(j-1))){
            subSyy <-  Syy[(P[i]+1):P[i+1],(P[j]+1):P[j+1]]
            f <-  t(W2[(P[i] + 1):P[i+1],i,drop=F])%*%subSyy%*%W2[(P[j] + 1):P[j+1],j,drop=F]
            f <- as.numeric(Setaeta[i,j])/as.numeric(f)
            Syy[(P[i]+1):P[i+1],(P[j]+1):P[j+1]] <- subSyy*f
            Syy[(P[j]+1):P[j+1],(P[i]+1):P[i+1]] <- t(subSyy*f)
            }
    }}

   Br1 <- B[(q1+1):q,1:q1,drop=FALSE]
   Br2 <- B[(q1+1):q,(q1+1):q,drop=FALSE]
   IB <- solve(diag(q2) - t(Br2))

    Sxy <- Sxx%*%W1%*%t(Br1)%*%IB%*%solve(Setaeta)%*%t(W2)%*%Syy
   S <- rbind(cbind(Sxx,Sxy),cbind(t(Sxy),Syy))

out <- list(S=S,B=B,Scomp=Sll,wx=rowSums(W1),wy=rowSums(W2))
return(out)
}
