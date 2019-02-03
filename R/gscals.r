## 
##  gscals.r                                                                  R. Schlittgen 27.11.2017
##
#' Estimating GSC models belonging to scenarios reflective-reflective, formative-reflective and formative-formative  
#'
#' \code{gscals} estimates GSC models  alternating least squares.
#' This leads to estimations of weights for the composites and an overall fit measure.
#'
#' @param  dat (n,p)-matrix, the values of the manifest variables.
#'              The columns must be arranged in that way that the components of refl are (absolutely) increasing.
#' @param B   (q,q) lower triangular matrix describing the interrelations of the latent variables:
#'             b_ij = 1 regression coefficient of eta_j in the regression relation in which eta_i is the depend variable
#'             b_ij = 0 if eta_i does not depend on eta_j in a direct way  (b_ii = 0 !)
#' @param  indicatorx vector describing with which exogenous composite the X-variables are connected
#' @param  indicatory vector describing with which endogenous composite the Y-variables are connected
#' @param  loadingx   logical TRUE when there are loadings for the X-variables in the model
#' @param  loadingy  logical TRUE when there are loadings for the Y-variables in the model
#' @param  maxiter Scalar, maximal number of iterations
#' @param  biascor  Boolean, FALSE if no bias correction is done, TRUE if parametric bootstrap bias correction is  done.
#' @return out list with components
#'        \tabular{ll}{
#'          Bhat \tab (q,q) lower triangular matrix with the estimated coefficients of the structural model \cr
#'          What \tab (n,q) matrix of weights for constructing the composites  \cr
#'          lambdahat \tab vector of length p with the loadings  or 0 \cr
#'          iter \tab number of iterations used   \cr
#'          fehl \tab maximal difference of parameter estimates for the last and second last iteration \cr
#'          composit \tab the data matrix of the composites \cr
#'          resid \tab the data matrix of the residuals of the structural model \cr
#'          S \tab the covariance matrix of the manifest variables \cr
#'          ziel \tab sum of squared residuals for the final sum \cr
#'          fit \tab The value of the fit criterion \cr
#'          R2 \tab vector with the coefficients of determination for all regression equations of the structural model
#'         }
#'
#' @examples
#' data(mobi250)
#' ind <- c(1, 1, 1, 4, 4, 4, 2, 2, 2, 3, 3, 5, 5, 5, 6, 6, 6, 7, 1, 1, 4, 4, 4, 4) 
#' o <- order(ind)
#' indicatorx <- c(1,1,1,1,1)
#' indicatory <- c(1, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 5, 5, 5)   
#' dat <- mobi250[,o]
#' dat <- dat[,-ncol(dat)]
#' B <- matrix(c(0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,
#'               0,1,1,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0),6,6,byrow=TRUE)
#' out <- gscals(dat,B,indicatorx,indicatory,loadingx=TRUE,loadingy=TRUE,maxiter=200,biascor=FALSE)
#'
#' @export
gscals <- function(dat,B,indicatorx,indicatory,loadingx=FALSE,loadingy=FALSE,maxiter=200,biascor=FALSE){

biascorstart <- biascor

p1 <- length(indicatorx)  # number of indicators of exogenous latent v.
p2 <- length(indicatory)  # number of indicators of endogenous latent v.
p <- p1+p2
q1<- length(unique(indicatorx))
q2<- length(unique(indicatory))
q <- q1+q2

n <- nrow(dat)
dat <- scale(dat,center=TRUE,scale=TRUE)
dat <- as.matrix(dat)

b <- rep(1,q)
equals0 <- function(x) all(x == 0)

if(ncol(dat) != p){ stop("length of c(indicatorx,indicatory) does not match number of columns of dat") }
if(!all(B[upper.tri(B,diag=T)]==0)){ stop("B must be lower triangular matrix with zeros on the diagonal") }
if( !(q==nrow(B)) ){ stop("number of composites do not match number of rows of B") }
if( length(b[apply(B,1,equals0)]) != q1 ){ stop("number of exogeneous latent variables do not match number of zero rows of B") }
if( (max(indicatorx) != q1)|(max(indicatory) != q2) ){ stop("inconsistency with indicatorx, indicatory, q1, q2") }

if( (loadingx == TRUE) & (loadingy == TRUE) ){ scenario <- "rr" }
if( (loadingx == FALSE) & (loadingy == TRUE) ){ scenario <- "fr" }
if( (loadingx == FALSE) & (loadingy == FALSE) ){ scenario <- "ff" }

# standardisation of the data

n <- nrow(dat)
dat <- scale(dat,center=TRUE,scale=TRUE)
dat <- as.matrix(dat)

R <- qr.R(qr(dat))
RR <- t(R)%*%R

# computation of starting values

Fehl <- NULL

Bhat <- B*matrix(runif(q^2),q,q)
lambda <- runif(p)
what <- runif(p)
R2 <- rep(NA,q2)

# construction of matrices W1, W2 and Lambdax, Lambday

Lx <- matrix(indicatorx,p1,q1)
Lx <- 1*(Lx == matrix(c(1:q1),p1,q1,byrow=T))

Ly <- matrix(indicatory,p2,q2)
Ly <- 1*(Ly == matrix(c(1:q2),p2,q2,byrow=T))

Lx0 <- Lx                           # for updating of A
Ly0 <- Ly                            # for updating of A
if( scenario == 2) { lambda[1:p1] <- 0 }
if( scenario == 3){ lambda[1:p] <- 0 }
Lambdax = Lx*lambda[1:p1]
Lambday = Ly*lambda[(p1+1):p]

x <- c(indicatorx,q1+indicatory)
W0 <- matrix(x,p,q)
W0 <- 1*(W0 == matrix(unique(x),p,q,byrow=T))
tW0 <- t(W0)
W <- W0*what
tWneu <-  t(W)
W1 <- W[,1:q1, drop=FALSE]
W2 <- W[,(q1+1):(q1+q2), drop=FALSE]

if( scenario == "rr"){                                                       # for reflective-reflective models


  V <- cbind(W2,diag(rep(1,p)))
  U <- W
  A <- t(Bhat[(q1+1):q,, drop=FALSE])
  A <- cbind(A,rbind(t(Lambdax),matrix(0,q2,p1)),rbind(matrix(0,q1,p2),t(Lambday)))
  Aalt <- A
  A0 <- t(B[(q1+1):q,, drop=FALSE])       # for updating of A
  A0 <- cbind(A0,rbind(t(Lx0),matrix(0,q2,p1)),rbind(matrix(0,q1,p2),t(Ly0)))
  Astar <- matrix(0,p,q2+p)
  Astar0 <- rbind(W0[1:p1,1:q1, drop=FALSE]%*%t(B[(q1+1):q,1:q1, drop=FALSE]),W0[(p1+1):p,(q1+1):q, drop=FALSE]%*%t(B[(q1+1):q,(q1+1):q, drop=FALSE]))
  Astar0 <- cbind(Astar0,rbind(W0[1:p1,1:q1, drop=FALSE]%*%t(Lx0),matrix(0,p2,p1)), rbind(matrix(0,p1,p2),W0[(p1+1):p,(q1+1):q, drop=FALSE]%*%t(Ly0)))
 }
else if( scenario == "fr" ){         # for formative-reflective models
  V <- cbind(W2,rbind(matrix(0,p1,p2),diag(rep(1,p2))))
  U <- W
  A <- rbind((t(Bhat[(q1+1):q,1:q1, drop=FALSE])),(t(Bhat[(q1+1):q,(q1+1):q, drop=FALSE])))
  A <- cbind(A,rbind(matrix(0,q1,p2),t(Lambday)))
  Aalt <- A
  A0 <- rbind(t(B[(q1+1):q,1:q1, drop=FALSE]),t(B[(q1+1):q,(q1+1):q, drop=FALSE]))        # für Aufdatierung von A
  A0 <- cbind(A0,rbind(matrix(0,q1,p2),t(Ly0)))
  Astar <- matrix(0,p,q2+p2)
  Astar0 <- rbind(W0[1:p1,1:q1]%*%t(B[(q1+1):q,1:q1, drop=FALSE]), W0[(p1+1):p,(q1+1):q, drop=FALSE]%*%t(B[(q1+1):q,(q1+1):q, drop=FALSE]))
  Astar0 <- cbind(Astar0,rbind(matrix(0,p1,p2),W0[(p1+1):p,(q1+1):q, drop=FALSE]%*%t(Ly0)))
  }
  else if(  scenario == "ff" ){                                              # for formative-formative models
  V <- W2
  U <- W
  A <- t(Bhat[(q1+1):(q1+q2), ,drop=FALSE])
  Aalt <- A
  A0 <- t(B[(q1+1):(q1+q2), ,drop=FALSE])
  tA0 <- t(A0)
  Astar <- matrix(0,p1+p2,q2)
  Astar0 <- rbind(W0[1:p1,1:q1]%*%t(B[(q1+1):q,1:q1, drop=FALSE]), W0[(p1+1):p,(q1+1):q, drop=FALSE]%*%t(B[(q1+1):q,(q1+1):q, drop=FALSE]))
 }
else { stop("only for scenario ff, scenario rr and  scenario fr") }

  fehl <- 1
  W1neu <- 0*W1
  W2neu <- 0*W2

 iter <- 0

# main loop
 while((iter < maxiter)&(fehl > 0.00000000001)){

   iter <- iter+1

   # update of A

   depv <- scale(dat%*%V)     # standardisation of composites instead of weights
   indv <- scale(dat%*%U)
   for(j in c(1:ncol(A))){
     aj <- A0[,j]
     aindj <- c(1:nrow(A))[aj != 0]
     out <- lm.fit(as.matrix(indv[,aindj]),depv[,j])
     A[aindj,j] <-  out$coefficients
   }
   fehlA <- max(abs(Aalt - A))
   Aalt <- A

   for (j in c(1:ncol(Astar0))){                # computation of Astar
     aj <- Astar0[,j]
     aindj <- c(1:nrow(Astar))[aj != 0]
     out <- lm.fit(dat[,aindj,drop=FALSE],depv[,j,drop=FALSE])
     Astar[aindj,j] <-  out$coefficients
   }
   tAstar <- t(Astar)
   tA <- t(A)


   # update of W1,W2

   if ( scenario == "rr"){                                                    # reflective-reflective models
      for(j in c(1:p)){
       aj <- c(1:q)[tW0[,j] == 1]
       out <- lm.fit(tA[,aj,drop=FALSE],tAstar[,j,drop=FALSE])
       tWneu[aj,j] <- out$coefficients
       }
      W1neu[1:p1,1:q1] <- t(tWneu[1:q1,1:p1, drop=FALSE])
      W2neu[(p1+1):p,] <- t(tWneu[(q1+1):q,(p1+1):p, drop=FALSE])
      }
   else if( scenario == "fr" ){                                            # formative-reflective models
      for(j in c(1:p)){
       aj <- c(1:q)[tW0[,j] == 1]
       out <- lm.fit(as.matrix(tA[,aj, drop=FALSE]),tAstar[,j, drop=FALSE])
       tWneu[aj,j] <-  out$coefficients
       }
      W1neu[1:p1,1:q1]  <- t(tWneu[1:q1,1:p1, drop=FALSE])
      W2neu[(p1+1):p,1:q2] <- t(tWneu[(q1+1):q,(p1+1):p, drop=FALSE])
      }
   else{                                                                         # scenario ff
      for(j in c(1:p1)){                                                      # weights for exogenous composites
       aj <- c(1:(q1+q2))[tW0[,j] == 1]
       out <- lm.fit(as.matrix(tA[,aj, drop=FALSE]),tAstar[,j, drop=FALSE])
       tWneu[aj,j] <-  out$coefficients
       }
      W1neu[1:p1,1:q1]  <- t(tWneu[1:q1,1:p1, drop=FALSE])
      W1neu <- scale(W1neu,center=FALSE,scale = apply(dat[,1:p1]%*%W1neu[1:p1,1:q1, drop=FALSE],2,sd))
      WA <- dat%*%as.matrix(cbind(W1neu,W2))%*%A                             # weights for endogenous l.v. 
     for(j in c(1:q2)){
      out  <- lm.fit(dat[,c(indicatorx,q1+indicatory) == (q1+j)],WA[,j])
      W2neu[c(indicatorx,q1+indicatory) ==  (q1+j),j] <- out$coefficients
      W2neu[1:p1,]  <- 0
    }
   }
  
   fehlW <- max(abs(cbind(W1-W1neu,W2-W2neu)))
   fehl <-  max(c(fehlA,fehlW))
   Fehl <- append(Fehl,fehl)

   # actualisation of W1,W2, V,U

   W1 <- W1neu
   W2 <- W2neu

   U <- cbind(W1,W2)
   if (scenario == "rr"){                                                       # reflective-reflective models
      V <- cbind(W2,diag(rep(1,p)))
      }
   else if( scenario == "fr" ){                                                # formative-reflective models
      V <- cbind(W2,rbind(matrix(0,p1,p2),diag(rep(1,p2))))
     }
   else{                                                                           # formative-formative models
      V  <- W2
     }
   }

   ziel <- sum(diag(t(V-U%*%A)%*%RR%*%(V-U%*%A)))
   fit <- 1 - ziel/sum(diag(t(V)%*%RR%*%V))

   # extraction of parameters

   if( scenario == "rr" ){                                                      # scenario  rr
     Bhat[(q1+1):q,1:q1] <- t(A[1:q1,1:q2])
     Bhat[(q1+1):q,(q1+1):q] <- t(A[(q1+1):q,1:q2])
     Lambdax <- t(A[1:q1,(q2+1):(q2+p1),drop=FALSE])
     Lambday <- t(A[(q1+1):q,(q2+p1+1):(q2+p1+p2),drop=FALSE])
     lambda <- c(rowSums(Lambdax),rowSums(Lambday))
   }else if(  scenario == "fr" ){                                             #  scenario fr
     Bhat[(q1+1):q,1:q1] <- t(A[1:q1,1:q2])
     Bhat[(q1+1):q,(q1+1):q] <- t(A[(q1+1):q,1:q2])
     Lambday <- t(A[(q1+1):q,(q2+1):q,drop=FALSE])
     lambda[(p1+1):p] <- rowSums(Lambday)
   }else{                                                                          #  scenario ff
     Bhat[(q1+1):q,1:q1] <- t(A[1:q1,1:q2])
     Bhat[(q1+1):q,(q1+1):q] <- t(A[(q1+1):q,1:q2])
     lambda <- numeric(p)
   }
What <- cbind(W1,W2) 

     Phi <- cor(dat[,1:p1]%*%W1[1:p1,])
     lambdax <- lambda[1:p1]
     lambday <- lambda[(p1+1):p]
    if( equals0(lambda[1:p1]) ){ lambdax <- NULL }
    if( equals0(lambday) ){ lambday <- NULL }
    if(q1 == 1){ wx <- W1[1:p1,] }
    else{  wx <- rowSums(W1[1:p1,]) }
    if(q2 == 1){ wy <- W2[(p1+1):p,] }
    else{  wy <- rowSums(W2[(p1+1):p,]) } 
    S <- gscmcov(Bhat,indicatorx,indicatory,lambdax=lambdax,lambday=lambday,wx=wx,wy=wy,Phi,R2=NULL) 
    S <- S$S 
    E <- eigen(S)                                               # für normalverteilung  
    if( !all(E$values  > 0)){ 
         if(biascor == TRUE){
         warning("Estimates result in a non positive definite covariance matrix of indicators. Bias correction not possibe.") 
         }
         biascor <- FALSE  }
 

 if(biascor == TRUE){
    Bo <- 50
    Bc <- array(NA, dim=c(q,q,Bo))
    Bd <- 0*B 
    lambdab <- 0
    Whatb <- 0     
    F <- E$vectors%*%diag(sqrt(E$values))                                  
    for(bo in c(1:Bo)){  
    X <- matrix(rnorm(n*p),n,p)              # uncorrelated normals  
    dat1 <- scale(X%*%t(F) )                     # correlated normals                          
    outb <- gscals(dat1,B,indicatorx,indicatory,loadingx=loadingx,loadingy=loadingy,maxiter=200,biascor=FALSE)
    Bc[,,bo] <- outb$Bhat
    lambdab <- lambdab + outb$lambdahat
    Whatb <- Whatb + outb$What
    }
    for(bo1 in c(1:q)){
      for(bo2 in c(1:q)){
      Bd[bo1,bo2] <- sum(Bc[bo1,bo2,])/Bo
    }}
    Bhat <- 2*Bhat - Bd
    lambda <- 2*lambda - lambdab/Bo
    What <- 2*What - Whatb/Bo
  }

composit <- as.matrix(scale(dat%*%What))
resid <- matrix(0,n,q2) 
for(j in c(1:q2)){  
   tB <- matrix(Bhat[q1+j,1:(q1+j-1)], (q1+j-1),1) 
   resid[,j] <- composit[,q1+j] -  composit[,1:(q1+j-1),drop=F]%*%tB
   R2[j] <- 1 - var(resid[,j])
   }
 
if(iter == maxiter){ warning("Maximum number of iterations, function may not have converged.", call. = FALSE) }

if( biascorstart != biascor ){ Bhat <- NA*Bhat } 
 
out <- list(Bhat=Bhat,What=What,lambdahat=lambda,iter=iter,fehl=Fehl,composit = composit,resid=resid,S=S,ziel=ziel,fit=fit,R2=R2)
return(out)
}

