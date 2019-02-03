##
##  boottestgscm.r                                                                  R. Schlittgen 1.12.2017   
##  

#' Testing two segmentations of a GSC model
#' 
#' \code{boottestgscm} computes a confidence interval for the difference of weighted average of averages of coefficients of determination 
#' for two segmentations of a  GSC model 
#' For a one sided alternative hypothesis the error alpha has to be duplicated
#' @param  dat (n,p)-matrix, the values of the manifest variables.    
#'              The columns must be arranged in that way that the components of refl are (absolutely) increasing. 
#' @param B   (q,q) lower triangular matrix describing the interrelations of the latent variables:
#'             b_ij = 1 regression coefficient of eta_j in the regression relation in which eta_i is the depend variable 
#'             b_ij = 0 if eta_i does not depend on eta_j in a direct way  (b_ii = 0 !)
#' @param  indicatorx vector describing with which exogenous composite the X-variables are connected 
#' @param  indicatory vector describing with which endogenous composite the Y-variables are connected 
#' @param  loadingx   logical FALSE when there are no loadings for the X-variables in the model 
#' @param  loadingy  logical FALSE when there are no loadings for the Y-variables in the model 
#' @param member1 vector of length n, indicating the cluster the observation belongs to for the first clustering
#' @param member2 vector of length n, indicating the cluster the observation belongs to for the second clustering
#' @param  alpha scalar, significance level ( = 1 - confidence level )
#' @param  inner Boolean, should a inner bootstrap loop be computed?
#' @return KI vector with the confidence bounds; positive lower limit indicates significant superiority of first
#'         clustering, negative upper limit of  second clustering. 
#'
#' @examples
#' \donttest{ 
#'   data(twoclm)
#'   member1 <- c(rep(1,50),rep(2,50))
#'   member2 <- twoclm[,10]
#'   dat <- twoclm[,-10]
#'   B <- matrix(c( 0,0,0, 0,0,0, 1,1,0),3,3,byrow=TRUE)
#'   indicatorx <- c(1,1,1,2,2,2)
#'   indicatory <- c(1,1,1)   
#'   boottestgscm(dat,B,indicatorx,indicatory,loadingx=FALSE,loadingy=FALSE,
#'                member2,member1,0.1,inner=FALSE)
#' }
#' @export

boottestgscm <- function(dat,B,indicatorx,indicatory,loadingx=FALSE,loadingy=FALSE,member1,member2,alpha,inner=FALSE){ 

n <- nrow(dat)
p <- ncol(dat) 
ind <- c(1:n)
k1 <- length(unique(member1)) 
k2 <- length(unique(member2)) 
 
q <- nrow(B)                                     # number of composites
b <- rep(1,q)
equals0 <- function(x) all(x == 0)
q1 <- length(b[apply(B,1,equals0)])    # number of exogeneous composites  
                                                        # measure for first clustering  
r1 <- averageR2w(dat,B,indicatorx,indicatory,loadingx=loadingx,loadingy=loadingy,member1)
                                                      # measure for second clustering  
r2 <- averageR2w(dat,B,indicatorx,indicatory,loadingx=loadingx,loadingy=loadingy,member2)

Boot <- 1000                                    # outer bootstrap-loop
tb <- numeric(Boot)
ustar <- numeric(Boot)
tstar <- numeric(Boot)
bo <- 0 
while( bo < Boot ){
    bo <- bo+1 
    indb <- sample(ind,n,replace=T)
    datb <- dat[indb,] 
    member1b <- member1[indb] 
    member2b <- member2[indb] 
                                                                       # to ensure that each cluster has enough observations  
    m1 <- m2 <- n
    for( i in c(1:k1)){ m1 <- min(m1,length(member1b[member1b == i] ) ) }
    for( i in c(1:k2)){ m2 <- min(m2,length(member2b[member2b == i] ) ) }
    if (min(c(m1,m2)) < (2*p)){ bo <- bo-1 }
    else{ 
                                                                                           # measure for first and second clustering
       r1b <-  averageR2w(datb,B,indicatorx,indicatory,loadingx=loadingx,loadingy=loadingy,member1b) 
       r2b <-  averageR2w(datb,B,indicatorx,indicatory,loadingx=loadingx,loadingy=loadingy,member2b)
 
          tstar[bo] <- r1b-r2b  
          
          if (inner){                                        # inner bootstrap loop, if inner == TRUE 
            Booti <- 1000        
            tstarstar <- numeric(Booti)  
            boi <- 0
            while(boi < Booti){ 
               boi <- boi+1 
                indbi <- sample(c(1:n),n,replace=TRUE) 
                datbi <- datb[indbi,]  
                member1bi <- member1b[indbi] 
                member2bi <- member2b[indbi] 
                # ensure, that all clusters have enough members 
               m1 <- m2 <- n
               for( i in c(1:k1)){ m1 <- min(m1,length(member1bi[member1bi == i] )) }
               for( i in c(1:k2)){ m2 <- min(m2,length(member2bi[member2bi == i] )) }  
                if (min(c(m1,m2)) < (2*p)){  boi <- boi-1 }
                else{
                                                                           #  measure for the first and second clustering         
                   r1bi <- averageR2w(datbi,B,indicatorx,indicatory,loadingx=loadingx,loadingy=loadingy,member1bi) 
                   r2bi <- averageR2w(datbi,B,indicatorx,indicatory,loadingx=loadingx,loadingy=loadingy,member2bi)
                   if( is.na(r1bi)|is.na(r2bi)  ){ boi <- boi-1  } 
                   else{    tstarstar[boi] = (r1bi-r2bi)     }     
                   } }  
                   ustar[bo] =  sum( 1*(tstarstar <= (2*tstar[bo]-(r1-r2))) )/Booti  
            } 
   }
}    # end of main loop 
tstar <- sort(tstar)  
if(inner){  
  ustar <- sort(ustar) 
  qhatu <- ustar[(Boot+1)*alpha/2] 
  qhato <- ustar[(Boot+1)*(1-alpha/2)]   
  au <-  tstar[min(max(((Boot+1)*qhatu),1)|Boot)]     # index becomes Boot+1 
  ao <-  tstar[max(min(((Boot+1)*qhato)|Boot),1)]  
}else{
  au <- tstar[(Boot+1)*alpha/2] 
  ao <- tstar[(Boot+1)*(1-alpha/2)]  
}
KI <- c(2*(r1-r2)-ao,2*(r1-r2)-au)
return(KI)
}

 
#' For use in boottestgscm. 
#' 
#' \code{averageR2w} computes the weighted average of average of coefficients of determination  for the structural parts 
#' of a segmented GSC model 
#' @param  dat (n,p)-matrix, the values of the manifest variables.    
#'              The columns must be arranged in that way that the components of refl are (absolutely) increasing. 
#' @param B   (q,q) lower triangular matrix describing the interrelations of the latent variables:
#'             b_ij = 1 regression coefficient of eta_j in the regression relation in which eta_i is the depend variable 
#'             b_ij = 0 if eta_i does not depend on eta_j in a direct way  (b_ii = 0 !)
#' @param  indicatorx vector describing with which exogenous composite the X-variables are connected 
#' @param  indicatory vector describing with which endogenous composite the Y-variables are connected 
#' @param  loadingx   logical TRUE when there are loadings for the X-variables in the model 
#' @param  loadingy  logical TRUE when there are loadings for the Y-variables in the model 
#' @param member vector of length n, indicating the cluster the observation belongs to  
#' @return r scalar, 'global' r2 coefficiet of determination
#'
#' @export
averageR2w <- function(dat,B,indicatorx,indicatory,loadingx=FALSE,loadingy=FALSE,member){

 n <- nrow(dat) 
 ind <- c(1:n) 
 mem  <- unique(member)
 k <- length(mem)   
 r <-  numeric(k)
 for(i in  c(1:k)){ 
     ind1 <- ind[member == mem[i]] 
     out <- gscals(dat[ind1,],B,indicatorx,indicatory,loadingx=loadingx,loadingy=loadingy) 
     r[i]  <- mean(out$R2)*length(ind1) 
     }
 r <- sum(r)/n                     # weighted mean of means of coefficients of determination   
return(r)
}

