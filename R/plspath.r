##
##  plspath.r                                                                  R. Schlittgen 15.06.2018
##  


#' Estimation of pls-path models  
#' 
#' \code{plspath} estimates pls path models using the classical approach formulated in Lohmueller.  
#'  
#' @param  dat (n,p)-matrix, the values of the manifest variables.
#'              The columns must be arranged in that way that the components of refl are (absolutely) increasing 
#' @param B   (q,q) lower triangular matrix describing the interrelations of the latent variables:
#'             b_ij= 1 regression coefficient of eta_j in the regression relation in which eta_i is 
# '            the depend variable
#'             b_ij= 0 if eta_i does not depend on eta_j in a direct way  (b_ii = 0 !)
#' @param  indicatorx (p1,1) vector indicating with which exogenous composite the x-indicators are related.
#' @param  indicatory (p2,1) vector indicating with which endogenous composite the y-indicators are related.
#'              The components of the indicators must be increasing. 
#' @param  modex equals "A" or "B" , the mode for this block of indicators
#' @param  modey equals "A" or "B" , the mode for this block of indicators
#' @param  maxiter Scalar, maximal number of iterations 
#' @param  stdev Boolean Should the standard deviations of the estimates be computed by bootstrap? 
#' @return out list wih components
#'     \tabular{ll}{
#'          Bhat \tab (q,q) lower triangular matrix with the estimated coefficients of the structural model \cr 
#'          eta \tab (n,q)-matrix, the scores of the latent variables \cr
#'          w \tab vector of length p of weights for constructing the latent variables \cr
#'          lambdahat \tab vector of length p with the loadings  \cr
#'          resa \tab (n,?) matrix of residuals from outer model  \cr
#'          resi \tab (n,?) matrix of residuals from inner model  \cr
#'          R2 \tab vector with the coefficients of determination for all regression equations of the structural model \cr
#'          iter \tab number of iterations used   \cr
#'          ret \tab scalar, return code:  \cr
#'              \tab  0 normal convergence  \cr
#'              \tab  1 limit of iterations attained, probably without convergence  \cr
#'          sdev.beta  \tab (q,q) matrix, the standard deviations of path coefficients (when stdev = TRUE) \cr
#'          sdev.lambda \tab vector, the standard deviations of loadings (when stdev = TRUE)  
#'          }
#'
#' @examples
#' data(mobi250)
#' refl <- c(1, 1, 1, 4, 4, 4, 2, 2, 2, 3, 3, 5, 5, 5, 6, 6, 6, 7, 1, 1, 4, 4, 4, 4)  
#' o <- order(refl)
#' dat <- mobi250[,o]
#' dat <- dat[,-ncol(dat)]
#' refl <- refl[o][-length(refl)] 
#' indicatorx <- refl[1:5]
#' indicatory <- refl[-c(1:5)] - 1
#' B <- matrix(c(0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,
#'               0,1,1,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0),6,6,byrow=TRUE)  
#' out <- plspath(dat,B,indicatorx,indicatory,modex="A",modey="A")
#' 
#' @export
plspath <- function(dat,B,indicatorx,indicatory,modex="A",modey="A",maxiter=100,stdev = FALSE){
 
# check of input

if(!(all(B[upper.tri(B, diag = T)] == 0))){ stop("B must be lower triangular matrix") } 
if(length(c(indicatorx,indicatory)) != ncol(dat)){ stop("length of indicatox,indicatory must be equal to columns of dat") } 

p <- ncol(dat)                                    # p = number of manifest variables 
n <- nrow(dat) 
q <- nrow(B)                                     # q = number of latent variables

ind <- c(1:p) 
ind1 <- c(1:q)  
q1 <- length(unique(indicatorx)) 
q2 <- length(unique(indicatory))    
refl <- c(indicatorx,(indicatory + q1))
eta <- matrix(0,n,q) 
lambda <- rep(0,p) 
Bhat <- B
resa <- NULL 
resi <- NULL 
R2 <- rep(NA,q)
tol <-  10^(-7)
ret <- 0 
sdev.lambda <- NULL
sdev.beta <- NULL
 
                                                    
# 0. step: centering and norming of the data  
  
d <- scale(dat,scale = TRUE,center=TRUE)         # centering and scaling of the data    
d <- as.matrix(d) 

# 1. step: choice of weights  

w <- matrix(1,p,1)                                               #  starting vector of outer weights 

# 2. step: outer approximation  

for( g in c(1:q)){
  i <- (g == refl) 
  iind <- ind[g == refl]                                            # indices of manif. variables connected with   
  eta[,g] <- d[,iind,drop=F]%*%w[iind,1]                  # latent variable g 
  eta[,g] <- scale(eta[,g],center=FALSE)                   # standardisation
  co <- cor(cbind(eta[,g],d[,iind]),use = "pairwise.complete.obs")  
  co <- co[c(2:(length(iind)+1)),1] 
  co <- sum(1*((co > 0)-(co <0 )))
if( !is.na(co) ){
  if(co < 0){ eta[,g] <- -1*eta[,g] }                           # mainly positiv correlations
  }
}  

#   main loop; comparison of w and walt   

walt <- 0 
it <- 0 
while((sum(abs(w-walt)) > tol) & (it < maxiter)){    # start of main loop 
 it <- it+1                                                 

# 3. step: saving of w at walt  

 walt <- w 

# 4. step: inner approximation  

 vorz <- var(eta)  
 vorz <- (vorz > 0)-( vorz < 0)  
 eta <- eta %*% (vorz * (B+t(B)))   
 eta <- scale(eta,center=FALSE)

# 5. step: recomputation of outer weights  
 
for( g in c(1:q)){
  iind <- ind[g == refl]                     # indices of manif. variables connected with  latent variable g                            
  if( ((g <= q1) & (modex == "B"))|((g > q1) & (modey == "B"))  ){         # mode B  
    out <- lm.fit(d[,iind],eta[,g])
    w[iind] <- out$coefficients    
  }else{                                             # mode A otherwise
    for( ii in c(1:length(iind))){  
      eg <- matrix(eta[,g],n,1)
      egd <- as.matrix(na.omit(cbind(eg,d[,iind[ii]])))
      out <- lm.fit(egd[,1,drop=F],egd[,-1,drop=F])
      w[iind[ii]] <- out$coefficients        
     } 
   } 
 } 
 
# 6. step: recomputation of the scores of the latent variables  

for( g in c(1:q)){
  iind <- ind[g == refl]                    # indices of manif. variables connected with  latent variable g   
  eta[,g] <- as.matrix(d[,iind])%*%as.matrix(w[iind])
   eta[,g] <- scale(eta[,g],center=FALSE)  
   co <- cor(cbind(eta[,g],d[,iind]))
   co <- co[2:(length(iind)+1),1] 
   co <- sum((co>0)-(co<0)) 
   if (co<0){ eta[,g] <- -1*eta[,g] }
}   
  
}                                                # end of main loop         

# 7. determination of loadings lambda  
 
for( g in c(1:q)){ 
  iind <- ind[g == refl]                    # indices of manif. variables connected with latent variable g  
  if( ((g <= q1) & (modex == "B"))|((g > q1) & (modey == "B"))  ){       # formativ groups 
    out <- lm.fit(d[,iind],eta[,g]) 
    lambda[iind] <- out$coefficients  
    resa <- cbind(resa,out$residuals) 
  }else{                                                                                                 # reflective groups   
    eg <- matrix(eta[,g],n,1)
    out <- lm.fit(eg,d[,iind])      
    lambda[iind] <- as.vector(out$coefficients)  
    resa <- cbind(resa,out$residuals)  
  }
}   


# 8. determination of path coefficients beta and R2  

for( g in c((q1+1):q)){  
   ii <- ind1[ B[g,] != 0]  
   out <- lm.fit(as.matrix(eta[,ii]),as.matrix(eta[,g]))  
   Bhat[g,ii] <- out$coefficients   
   resi <- cbind(resi,out$residuals) 
   R2[g] <- var(out$fitted.values)/var(eta[,g])  
}     

if ((sum(abs(w-walt))>0.0001) & (it >= maxiter)){ ret <- 1 }

if(stdev == TRUE){
    Bo <- 100
    Bc <- array(NA, dim=c(q,q,Bo))
    Bd <- 0*B  
    lambdab <- matrix(NA,p,Bo)  
    for(bo in c(1:Bo)){
      dat1 <- dat[sample(1:n,n,replace=T),]
      out <- plspath(dat1,B,refl)
      Bb[,,bo] <- out$Bhat
      lambdab[,bo] <- out$lambda
      }
    sdev.lambda <- apply(lambdab,1,sd)
    sdev.beta <- matrix(0,q,q)
    for( i in (1:q)){
       for(j in c(1:q)){
         sdev.beta[i,j] <- sd(sdev.beta[i,j,])
       }}
     }   
   
return(list(Bhat=Bhat,eta=eta,w=w,lambdahat=lambda,resa=resa,resi=resi,R2=R2,iter=it,ret=ret,sdev.beta=sdev.beta,sdev.lambda=sdev.lambda) )
} 
