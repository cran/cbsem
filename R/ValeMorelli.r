
##
##  ValeMorelli.r                                                                  R. Schlittgen 15.06.2018
##
#'   Functions to generate nonnormal distributed multivariate random vectors with mean=0, var=1 and  given correlations and coefficients of skewness
#'   and excess kurtosis. This is done with the method of Vale & Morelli: The coefficients of the Fleishman transform Y = -c  + bX +cX^2 + dX^3 are computed.
#'   from given skewness gamma[1] = E(Y^3) and kurtosis gamma[2] = E(Y^4) - 3. A indermediate correlation matrix is computed from the desired correlation 
#'   matrix and the Fleishman coefficients. A singular value decomposition of the indermediate correlation matrix is performed and  a matrix of independend 
#'   normal random numbers is generated and transformed into correlated ones. Finally the Fleishman transform is applied to the columns of this data matrix.   
#'
#'   The function are adapted from online support of the SAS system,
#'   URL: support.sas.com/publishing/authors/extras/65378_Appendix_D_Functions_for_Simulating_Data_by_Using_Fleishmans_Transformation.pdf


#' \code{FleishmanIC} produce an initial guess of the Fleishman coefficients from given skewness and kurtosis. It is to use for Newton's algorithm. 
#'  This guess is produced by a polynomial regression. 
#'    
#' @param  skew desired skewness 
#' @param  kurt desired kurtosis 
#' @return par vector with coefficients b,c,d  
#'
#' @examples 
#' out <- FleishmanIC(1,2)
#'
#' @export
FleishmanIC  <- function(skew,kurt){
par <- numeric(3)
par[1] <- 0.95357 - 0.05679*kurt +  0.03520*skew^2 + 0.00133*kurt^2         # c1 = Quad(skew, kurt)   
par[2] <- 0.10007*skew + 0.00844*skew^3                                                   # c2 = Cubic(skew)  
par[3] <- 0.30978 -0.31655 *par[1]                                                               # c3 = Linear(c1) 
return (par)
}          


#' \code{Fleishman} computes the variance, skewness and kurtosis for a given set of of coefficients b,c,d for the Fleishman transform
#'
#' @param  coef vector with the coefficents 
#' @return out vector with coefficients Var,Skew,Kurt  
#'
#' @examples
#' coef <- c( 0.90475830, 0.14721082, 0.02386092) 
#' out <- Fleishman( coef )
#'
#' @export
Fleishman  <- function( coef ){
b <- coef[1] ; c <-coef[2] ; d <- coef[3]
Var <- b^2 + 6*b*d + 2*c^2 + 15*d^2
Skew <- 2*c*(b^2  + 24*b*d +105*d^2 + 2) 
Kurt <- 24*( b*d + c^2*(1+b^2+28*b*d) + d^2*(12+48*b*d +141*c^2 + 225*d^2) ) 
return(c(Var,Skew,Kurt))
}
 

# 
#' \code{FlDeriv}compute the Jacobian of the Fleishman transform  for a given set of coefficients b,c,d   
#'
#' @param  coef vector with the coefficents for the Fleishman transform 
#' @return  J  (3,3) Jacobian matrix of partial derivatives
#'
#' @examples
#' coef <- c( 0.90475830, 0.14721082, 0.02386092)  
#' J <- FlDeriv( coef )
#'
#' @export
FlDeriv <- function( coef ){ 
b <- coef[1]; c <- coef[2]; d <- coef[3] 
b2 <- b^2; c2 <- c^2; d2 <- d^2; bd <- b*d
df1db <- 2*b + 6*d
df1dc <- 4*c
df1dd <- 6*b + 30*d
df2db <- 4*c*(b + 12*d)
df2dc <- 2*(b2 + 24*bd + 105*d2 + 2)
df2dd <- 4*c*(12*b + 105*d)
df3db <- 24*(d + c2*(2*b + 28*d) + 48*d^3)
df3dc <- 48*c*(1 + b2 + 28*bd + 141*d2)
df3dd <- 24*(b + 28*b*c2 + 2*d*(12 + 48*bd + 141*c2 + 225*d2)+ d2*(48*b + 450*d))
J <-  matrix(c(df1db , df1dc , df1dd  ,
                       df2db , df2dc , df2dd ,
                       df3db , df3dc , df3dd),3,3,byrow=T)
return( J )
}


#' \code{NewtonFl} Newton's method to find roots of  the function FlFunc.   
#'
#' @param  target vector with the desired skewness and kurtosis
#' @param  startv vector with initial guess of the coefficents for the Fleishman transform 
#' @param maxIter maximum of iterations 
#' @param converge limit of allowed absolute error
#' @return  out  list with  components                 
#'        \tabular{ll}{
#'          coefficients \tab vector with the approximation to the root    \cr
#'          value \tab vector with differences of root and target  \cr 
#'          iter \tab number of iterations used   \cr  
#'          }
#'
#' @examples
#' skew <- 1; kurt <- 2
#' startv <- c( 0.90475830, 0.14721082, 0.02386092)  
#' out <- NewtonFl(c(skew,kurt),startv) 
#'
#' @export 
NewtonFl <- function(target,startv, maxIter = 100, converge = 1e-12 ){ 
x <- startv
 
FlFunc <- function(coef,target){  return(  Fleishman(coef) - c(1,target) ) } 
f <- FlFunc(x,target)
iter <- 0
while( (iter <= maxIter) & (max(abs(f)) > converge)){
iter <- iter + 1
J <- FlDeriv(x)
delta <- -solve(J, f)                   #  correction vector  
x <- x + delta                          # new approximation  
f <- FlFunc(x,target)
}                                            
if (iter > maxIter){ f <- NA }      # return missing if no convergence 
out <- list(coefficients=x, value=f, iter = iter)
return(out)
}



##########    Vale-Maurelli method of generating multivariate nonnormal data  
#  

#' \code{SolveCorr} Solve the Vale-Maurelli cubic equation to find the intermediate correlation between two normal variables that gives rise to a target
#'  correlation (rho) between the two transformed nonnormal variables.   
#'
#' @param  rho desired correlation of transformed variables
#' @param  coef1 vector with coefficents for the Fleishman transform of the first variable
#' @param  coef2 vector with coefficents for the Fleishman transform of the second variable 
#' @return  root the intermediate correlation
#'
#' @examples
#' rho <- 0.5
#' coef1<-  c( 0.90475830, 0.14721082, 0.02386092)   
#' coef2<-  c( 0.90475830, 0.14721082, 0.02386092)   
#' r <- SolveCorr(rho, coef1, coef2) 
#'
#' @export 
SolveCorr <- function(rho, coef1, coef2){
b1 <- coef1[1]; c1 <- coef1[2]; d1 <- coef1[3]; a1 <- -c1 
b2 <- coef2[1]; c2 <- coef2[2]; d2 <- coef2[3]; a2 <- -c2 
coef <-  c(-rho, (b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2), (2 * c1 * c2), (6*d1*d2))
roots <- polyroot(coef)                                 # solve for zero of cubic polynomial  
roots <-  Re(roots[abs(Im(roots)) < 1e-8])      # extract the real root(s)  
i <-  which.min(abs(roots))                             # return smallest real root 
return ( roots[i] ) 
}
   
#' \code{VMTargetCorr} Given a target correlation matrix, R, and target values of skewness and kurtosis for each marginal distribution, 
#' find the "intermediate" correlation matrix, V     
#'
#' @param  R desired correlation matrix of transformed variables
#' @param  Fcoef either vector with coefficents for the Fleishman transform to be applied to all variables or
#'               (nrow(R),3) matrix with different coefficients 
#' @return  V  the intermediate correlation matrix
#'
#' @examples
#' R <- matrix(c(1, 0.5, 0.3, 0.5 ,1, 0.2 , 0.3, 0.2 , 1),3,3)
#' coef <-  matrix(c( 0.90475830, 0.14721082, 0.02386092,0.78999781,0.57487681,
#'                             -0.05473674,0.79338100, 0.05859729, 0.06363759 ),3,3,byrow=TRUE) 
#' V <- VMTargetCorr(R, coef) 
#'
#' @export 
VMTargetCorr <- function(R, Fcoef){
  if( !is.matrix(Fcoef)){ matrix(Fcoef,nrow(R),3,byrow=TRUE)  }
  V <- diag(nrow(R))  
  for( i in c(2:nrow(R))){ 
     for( j in c(1:(i-1))){ 
     V[i,j] <- SolveCorr(R[i,j], Fcoef[i,], Fcoef[j,]) 
     V[j,i] <- V[i,j] 
     } } 
return (V) 
}

# Simulate data from a multivariate nonnormal distribution such that 


#' \code{rValeMaurelli} Simulate data from a multivariate nonnormal distribution such that
#' 1) Each marginal distribution has a specified skewness and kurtosis
#' 2) The marginal variables have the correlation matrix R     
#'
#' @param  n number of random vectors to be generated
#' @param  R desired correlation matrix of transformed variables
#' @param  Fcoef either vector with coefficents for the Fleishman transform to be applied to all variables or
#'               (nrow(R),3) matrix with different coefficients 
#' @return  X  (n,nrow(R)) data matrix
#'
#' @examples
#' R <- matrix(c(1, 0.5, 0.3, 0.5 ,1, 0.2 , 0.3, 0.2 , 1),3,3)
#' coef <-  matrix(c( 0.90475830, 0.14721082, 0.02386092,0.78999781,0.57487681,
#'                             -0.05473674,0.79338100, 0.05859729, 0.06363759 ),3,3,byrow=TRUE) 
#' V <- rValeMaurelli(50, R, coef) 
#'
#' @export  
rValeMaurelli <- function(n, R, Fcoef){
  
  # adjust correlation matrix 
  p <- nrow(R)
  V <- diag(p)  
  for( i in c(2:p)){ 
     for( j in c(1:(i-1))){ 
     V[i,j] <- SolveCorr(R[i,j], Fcoef, Fcoef) 
     V[j,i] <- V[i,j] 
     } }   
     out <- eigen(V)
     F <- out$vectors%*%diag(sqrt(out$values)) 
      X <- matrix(rnorm(n*p),n,p)     # uncorrelated normals  
      Y <- scale(X%*%t(F) )                 # correlated normals
    for( j in  c(1:p)){ 
       w = Y[,j]                                  # apply Fleishman transformation 
       X[,j] <- -Fcoef[2] + w*(Fcoef[1] + w*(Fcoef[2] + w*Fcoef[3])) 
      }
return(X) 
}                                                                                                      