##
##  checkw.r                                                             R. Schlittgen 27.11.2017
##
#' Checking composite based SE models if there are weights in accordance with the loadings and the 
#' covariance matrix of the composites
#'
#' \code{checkw} determines if there are sets of weights fulfilling the critical relation for the 
#' covariance matricies of the composites. 
#' 
#' @param B   (q,q) lower triangular matrix describing the interrelations of the latent variables:
#'             b_ij = 1 regression coefficient of eta_j in the regression relation in which eta_i 
#'                      is the depend variable
#'             b_ij = 0 if eta_i does not depend on eta_j in a direct way  (b_ii = 0 !)
#' @param  indicatorx vector describing with which exogenous composite the X-variables are connected
#' @param  indicatory vector describing with which endogenous composite the Y-variables are connected
#' @param  lambdax   vector of loadings  for the X-variables in the model or FALSE
#' @param  lambday  vector of loadings for the Y-variables in the model or FALSE
#' @param  wx vector of weights  for the X-variables in the model or FALSE
#' @param  wy vector of weights  for the Y-variables in the model or FALSE
#' @param  Sxixi covariance matrix of exogenous composites
#' @param  R2 vector of coefficients of determination of structural regression equations
#' @return out list with components
#'        \tabular{ll}{
#'          crit.value \tab vector of length 2 with the values of the optimisation criterion \cr
#'          wx \tab  vector of length p1 of weights for constructing the exogenous composites  \cr
#'          wy \tab  vector of length  p2 of weights for constructing the endogenous composites  \cr
#'         } 
#' @examples    
#' B <- matrix(c(0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,
#'               0,1,1,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0),6,6,byrow=TRUE)
#' indicatorx <- c(1,1,1,1,1)
#' indicatory <- c(1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5) 
#' lambdax <- c(0.73, 0.60, 0.60, 0.77, 0.74) 
#' lambday <- c(0.79, 0.68, 0.60, 0.90, 0.94, 0.80, 0.65, 0.78, 0.78, 0.74, 
#'                      0.77, 0.78, 0.80, 0.84, 0.85, 0.86, 0.23, 0.87) 
#' Sxixi <- matrix(1,1,1)
#' out <- checkw(B,indicatorx,indicatory,lambdax=TRUE,lambday=TRUE,wx=FALSE,wy=FALSE, Sxixi,R2=NULL)
#'
#' @export
checkw <- function(B,indicatorx,indicatory,lambdax=FALSE,lambday=FALSE,wx=FALSE,wy=FALSE,Sxixi,R2=NULL){   

# check, if there is a reflective relation 
if( all(lambdax == FALSE) & all(lambday == FALSE) ){ stop("no reflective relation, nothing to do") }
  
crit.value <- c(NA,NA) 
q <- nrow(B)                                     # number of latent variables
b <- rep(1,q)
equals0 <- function(x) all(x == 0)
q1 <- length(b[apply(B,1,equals0)])     # number of exogeneous composites
q2 <- q - q1                                        # number of endogeneous composites
p1 <- length(indicatorx)                      # number of indicators of exogeneous composites
p2 <- length(indicatory)                      # number of indicators of endogeneous composites
p <- p1 + p2                                       # number of manifest variables     
if( is.numeric(wx)  == FALSE){ wx <- NA }
if( is.numeric(wy)  == FALSE){ wy <- NA } 
 
                                                               # check for Wx if lambdax is given
                                                               # construction of Lx
   Lx <- matrix(indicatorx,p1,q1) 
   Lx <- 1*(Lx == matrix(c(1:q1),p1,q1,byrow=T))   
   Lx <- Lx*lambdax   
   Sxx <- Lx%*%Sxixi%*%t(Lx) 
   Sdd  <-  diag(1 - diag(Sxx))
   diag(Sxx) <- 1
 
   # check of existence of Wx 
   
   W1 <- matrix(indicatorx ,p1,q1)      # construction of W1 
   W1 <- 1*(W1 == matrix(c(1:q1),p1,q1,byrow=T))  
   
   if(is.numeric(wx) == FALSE){
   Wx <- scale(W1,center = F, scale= sqrt(diag(t(W1)%*%Sxx%*%W1))) 
   wstart <- rowSums(Wx)           # construction of wstart 
   out1  <- optim(wstart, subcheckw, gr = NULL, indicatorx,Sdd,Lx,Sxixi, method="L-BFGS-B" ,
                           lower=-1,upper=1)  
   crit.value[1] <- out1$value       
   wx <- out1$par     
   }           
   Wx <- W1*wx
                         
# scaling of path coefficients and computing Setaeta
 
   Sll <- diag(q1+q2)
   Sll[1:q1,1:q1] <- Sxixi
   for (m in c((q1+1):(q1+q2))){  
    if (length(R2) > 0){ 
      tau <- sqrt(R2/(B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),1:(m-1)]%*%t(B[m,1:(m-1),drop=F])) )
      B[m,] <- tau*B[m,] 
      }
      for (j in c(1:m)){ 
             Sll[j,m] <- Sll[m,j] <- B[m,1:(m-1),drop=F]%*%Sll[1:(m-1),j,drop=F] 
         diag(Sll) <- 1  
       }  }     
   Setaeta <- Sll[(q1+1):q,(q1+1):q] 

###################
 
   # check of existence of Wy              
                                                                   # construction of Ly  
   Ly <- matrix(indicatory,p2,q2)
   Ly <- 1*(Ly == matrix(c(1:q2),p2,q2,byrow=T))  
   Ly <- Ly*lambday      
   Syy <- Ly%*%Setaeta%*%t(Ly) 
   See <- diag(1-diag(Syy))
   diag(Syy) <-1   
   
   W2 <- matrix(indicatory ,p2,q2)      # construction of W1 
   W2 <- 1*(W2 == matrix(c(1:q2),p2,q2,byrow=T))  
   if(is.numeric(wy) == FALSE){
     Wy <- scale(W2,center = F, scale= sqrt(diag(t(W2)%*%Syy%*%W2))) 
     wstart <- rowSums(Wy)           # construction of wstart 
     out2  <- optim(wstart, subcheckw, gr = NULL,indicatory,See,Ly,Setaeta, method="L-BFGS-B" , 
                             lower=-1,upper=1)   
     crit.value[2] <- out2$value
     wy <- out2$par  
   } 
   W2 <- W2*wy              
   
################     
  
  out <- list(wx = wx, wy = wy, crit.value = crit.value)
 
return(out)
}
 
#' Function for use in checkw
#'
#' \code{subcheckw} computes the sum of squared differences of two formulas for the covariancematrix of composites
#' 
#' @param w vector of weights
#' @param  indicator  vector describing with which exogenous composite the indicators are connected 
#' @param  S covariance matrix of errors resulling from regession for loadings
#' @param  L matrix of loadings
#' @param  Scomp covariance matrix of composites 
#' @return out  scalar, sum of squared differences
#' 
#' @export
subcheckw <- function(w,indicator,S,L,Scomp){
   p <- length(indicator)
   q <- length(unique(indicator))   
   W  <- matrix(indicator,p,q)      # construction of W 
   W  <- 1*(W  == matrix(c(1:q),p,q,byrow=T)) 
   W <- W*w 
  out <- sum((t(W)%*%(L%*%Scomp%*%t(L) + S)%*%W - Scomp)^2)  
  return(out)  
   }
