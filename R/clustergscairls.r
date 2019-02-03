##
##  clustergscairls.r                                                  R. Schlittgen 29.11.2017
##

#' Clustering gsc-models
#'
#' \code{clustergscairls} clusters data sets in that way that each cluster has a its own
#'  set of coefficients in the gsc-model.
#'
#' @param dat (n,p)-matrix, the values of the manifest variables
#' @param B   (q,q) lower triangular matrix describing the interrelations of the latent variables:
#'             b_ij= 1 regression coefficient of eta_j in the regression relation in which eta_i is
# '            the depend variable
#'             b_ij= 0 if eta_i does not depend on eta_j in a direct way  (b_ii = 0 !)
#' @param  indicatorx vector describing with which exogenous composite the X-variables are connected
#' @param  indicatory vector describing with which endogenous composite the Y-variables are connected
#' @param  loadingx  logical TRUE when there are loadings for the X-variables in the model
#' @param  loadingy  logical TRUE when there are loadings for the Y-variables in the model  
#' @param  k  scalar, the number of clusters to be found
#' @param  minmem number of the cluster's members or FALSE (then ist is set to 2*number of indicators)
#' @param  wieder  scalar, the number of random starts
#'
#' @return out list with components
#'   \tabular{ll}{
#'               member \tab (n,1)-vector, indicator of membership  \cr
#'               Bhat \tab (k,q,q)-array, the path coefficients of the clusters \cr
#'               lambda \tab (p,k)-matrix, the loadings of the clusters \cr
#'               fitall \tab the total fit measure for the structural models only \cr
#'               fit \tab vector of length k, the fit values of the different models \cr
#'               R2 \tab (k,q) matrix, the coefficients of determination for the structural regression equations
#'               }
#' @examples
#' \donttest{ 
#' data(twoclm)
#' dat <- twoclm[,-10]
#' B <- matrix(c( 0,0,0, 0,0,0, 1,1,0),3,3,byrow=TRUE)
#' indicatorx <- c(1,1,1,2,2,2)
#' indicatory <- c(1,1,1)  
#' out <- clustergscairls(dat,B,indicatorx,indicatory,loadingx=FALSE,loadingy=FALSE,2,minmem=6,1) 
#' }
#' @export
clustergscairls <- function(dat,B,indicatorx,indicatory,loadingx=FALSE,loadingy=FALSE,k,minmem=FALSE,wieder){

# standardisation of the data

dat <- scale(dat,center=TRUE,scale=TRUE)
dat <- as.matrix(dat)

n <- nrow(dat)
p1 <- length(indicatorx)                     # number of indicators of exogenous latent v.
p2 <- length(indicatory)                     # number of indicators of endogenous latent v.
p <- p1+p2                                       # number of manifest variables
q1<- length(unique(indicatorx))          # number of  exogenous latent v.
q2<- length(unique(indicatory))          # number of endogenous latent v.
q <- q1+q2                                       # number of latent variables

if(minmem == FALSE){minmem  <- 2*p }                                # minimal number of members fo a cluster
refl <- c(indicatorx,q1+indicatory)
W0 <- matrix(refl,p,q)
W0 <- 1*(W0 == matrix(c(1:q),p,q,byrow=TRUE))
W2 <- W0[,(q1+1):q,drop=FALSE] 
ind <- c(1:n)
indlv <- c(1:q)
indlvreg <- c((q1+1):q)                  # indices of endogenous latent variables
ni <- rep(0,k)

member <- matrix(0,n,wieder)
fit <- numeric(wieder)

#  generate "wieder" solutions with random starts

for( wie in (1:wieder)){

  weights <- matrix(1/(k*100),n,k)
  for( i in c(1:n)){  weights[i,sample(c(1:k),1)] <- 1-(k-1)/(k*100) }
  rest <- matrix(0,n,k)
  coeff <- matrix(NA,length(refl),k)
  param <-  matrix(1,q^2 + 2*length(refl),k)
  differenz <-  rep(1,k)

  B2 <-  array(NA, dim=c(q,q,k))
  lambda <- matrix(NA,p,k)
  ss2 <- matrix(NA,p,k)
  R2 <- matrix(NA,k,q2)
  R2w <- rep(0,q2)
  y <- matrix(0,n,q)
  w <- matrix(0,n,q)
  zaehl <- numeric(k)
  nenn <- numeric(k)

  ite <- 0
  while ((ite<100) & (sum(abs(differenz)) > 0.0001) ){                        # iteration for adjusting the weights
      ite <- ite+1
      oldparam <- param
      for( j in c(1:k)){                                                    # estimation of separate models for the k sets of weights 
        out <- gscals(dat*weights[,j],B,indicatorx,indicatory,loadingx=loadingx,loadingy=loadingy)
                                                                                                        # compute residuals without weighting
       resi <-    gscalsresid(dat,out,indicatorx,indicatory,loadingx,loadingy) 
        rest[,j] <- apply((resi[,1:q2,drop=FALSE])^2,1,sum)                      # only residuals of the inner model
        param[,j] <- c(as.vector(out$Bhat),out$lambdahat,as.vector(rowSums(out$What)))
        differenz[j] <- sum(abs(param[1:(q^2),j]-oldparam[1:(q^2),j]))
      }
      ares <- (rest + 1e-8)^-1
      weights <- ares/apply(ares,1,sum)
  }
  member[,wie] <-  apply(weights,1,which.max)                  # assignment of membership from the models found

  # coefficients of determination / sums of squared residuals ( models without weighting)

  zaehl <- rep(NA,k)
  nenn <- rep(NA,k)
  for( i in c(1:k)){
   indi <- ind[member[,wie] == i]
   dat1 <- dat[indi,]
     out <- gscals(dat1,B,indicatorx,indicatory,loadingx=loadingx,loadingy=loadingy)
     zaehl[i] <- sum(diag(t(out$resid[,1:q2,drop=F])%*%out$resid[,1:q2,drop=F]))  #  fit for stucture model only
     zwi <- dat%*%(W2*rowSums(out$What))
     nenn[i] <- sum(diag(t(zwi)%*%zwi))

   }
  fit[wie] <- 1- sum(zaehl, na.rm = T)/sum(nenn, na.rm = T)
}                                                                                    # end of iteration "wieder"

member <- member[,which.max(fit)]                               # final cluster-membership

# computations for return

zaehl <- rep(NA,k)
nenn <- rep(NA,k)
fit <- rep(NA,k)

for (i in c(1:k)){
 indi <- ind[member==i]
 dat1 <- dat[member == i,]
 ni[i] <- nrow(dat1)
 if( nrow(dat1) > minmem ){
   out <- gscals(dat1,B,indicatorx,indicatory,loadingx=loadingx,loadingy=loadingy)
   R2[i,] <- out$R2
   lambda[,i] <- out$lambdahat
   B2[,,i] <- out$Bhat
   fit[i] <- out$fit
   zaehl[i] <- sum(diag(t(out$resid[,1:q2,drop=F])%*%out$resid[,1:q2,drop=F]))  #  fit for stuctural model only
   zwi <- dat%*%(W2*rowSums(out$What))
   nenn[i] <- sum(diag(t(zwi)%*%zwi))
   # zaehl[i] <- out$ziel[2]
   # nenn[i] <- out$ziel[2]/(1-out$fit)
  }
}

fitall <-  1- sum(zaehl, na.rm = T)/sum(nenn, na.rm = T)

return(list(member=member,Bhat=B2,lambda=lambda,fitall=fitall,fit=fit,R2=R2))
}


#' For use in clustergscairls, residuals of a gsc-model
#'
#' \code{gscalsresid} computes the residuals of a gsc-model when the parameters and weights are given
#'
#' @param  dat   (n,p) data matrix
#' @param  out  list, output from gscals
#' @param  indicatorx vector describing with which exogenous composite the X-variables are connected
#' @param  indicatory vector describing with which endogenous composite the Y-variables are connected
#' @param  loadingx logical TRUE when there are loadings for the X-variables in the model
#' @param  loadingy logical TRUE when there are loadings for the y-variables in the model
#'
#' @return   resid (n,q2) matrix of residuals from structural model, the q2 is the number of endogenous composites .
#'
#' @examples
#' data(simplemodel) 
#' data(gscalsout)
#' B <- matrix(c( 0,0,0, 0,0,0, 0.7,0.4,0),3,3,byrow=TRUE)
#' indicatorx <- c(1,1,1,2,2,2)
#' indicatory <- c(1,1,1)  
#' out <- gscalsresid(simplemodel,gscalsout,indicatorx,indicatory,TRUE,TRUE)
#'
#' @export
gscalsresid <- function(dat,out,indicatorx,indicatory,loadingx,loadingy){
           Bhat <- out$Bhat
           W <- out$What
           lambda <- out$lambdahat
if( (loadingx == TRUE) & (loadingy == TRUE) ){ scenario <- 1 }
if( (loadingx == FALSE) & (loadingy == TRUE) ){ scenario <- 2 }
if( (loadingx == FALSE) & (loadingy == FALSE) ){ scenario <- 3 }

# standardisation of the data

n <- nrow(dat)
dat <- scale(dat,center=TRUE,scale=TRUE)
dat <- as.matrix(dat)

# construction of matricies W1, W2 and Lambdax, Lambday

p1 <- length(indicatorx)                     # number of indicators of exogenous latent v.
p2 <- length(indicatory)                     # number of indicators of endogenous latent v.
p <- p1+p2                                       # number of manifest variables
q1<- length(unique(indicatorx))          # number of  exogenous latent v.
q2<- length(unique(indicatory))          # number of endogenous latent v.
q <- q1+q2                                       # number of latent variables


x <- indicatorx
Lx <- matrix(indicatorx,p1,q1)
Lx <- 1*(Lx == matrix(c(1:q1),p1,q1,byrow=T))
Ly <- matrix(indicatory,p2,q2)
Ly <- 1*(Ly == matrix(c(1:q2),p2,q2,byrow=T))
Lambdax = Lx*lambda[1:p1]
Lambday = Ly*lambda[(p1+1):p]

refl <- c(indicatorx,q1+indicatory)
W1 <- W[,1:q1, drop=FALSE]
W2 <- W[,(q1+1):(q1+q2), drop=FALSE]

if(scenario == 1){                                                          # for scenario1
  V <- cbind(W2,diag(rep(1,p)))
  U <- W
  A <- t(Bhat[(q1+1):(q1+q2),, drop=FALSE])
  A <- cbind(A,rbind(t(Lambdax),matrix(0,q2,p1)),rbind(matrix(0,q1,p2),t(Lambday)))
 }
else if( scenario ==2 ){                                                  # for scenario2
  V <- cbind(W2,rbind(matrix(0,p1,p2),diag(rep(1,p2))))
  U <- cbind(W1,W2)
  A <- rbind((t(Bhat[(q1+1):(q1+q2),1:q1, drop=FALSE])),(t(Bhat[(q1+1):(q1+q2),(q1+1):(q1+q2), drop=FALSE])))
  A <- cbind(A,rbind(matrix(0,q1,p2),t(Lambday)))
}
else{                                                                            # for scenario1
   V <- W2
   U <- W
   A <- t(Bhat[(q1+1):(q1+q2), ,drop=FALSE])
   }
res <- dat%*%(V-U%*%A)

return(res)
}
