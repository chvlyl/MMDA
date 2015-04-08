#' Model-based multi-sample approach for metagenomic differential abundance analysis.
#'
#' Differential abundance analysis for metagenomic data
#' @param X the data matrix. The columns must be the markers, the rows must be the samples
#' @param Z the covariates. The length must be the same the the rows of the data matrix X
#' @param tc total number of read counts for each sample
#' @param l marker length
#' @return test.table  The table contains the estimated coefficients and pvalues for b
#' @export
#' @examples
#' mada(Data,Covariates,total.counts,marker.length)

mada <- function(X,Z,tc,l,estimate.phi=TRUE,verbose=FALSE){
    tc <- as.matrix(tc)
    l  <- as.matrix(l)
    X  <- as.matrix(X)
    Z  <- as.matrix(Z)
    if (dim(tc)[1] != nrow(X) | dim(tc)[2] != 1){
      stop('The length of tc is not right!')
    }
    if (dim(l)[1] != ncol(X) | dim(l)[2]!= 1){
      stop('The length of l is not right!')
    }
    if(nrow(X) <= 1){
      stop('Too few samples')
    }
    if(ncol(X) <= 2){
      stop('Twoo few markers')
    }
    if(any(colSums(X,na.rm=TRUE) == 0)){
      stop('Find maker with all zero counts across all samples')
    }
    #if(any(rowSums(X) == 0)){
    #  stop('Find sample with all zero counts across all markers')
    #}
 
    LRT <- likelihood_ratio_test_mada(X,Z,tc,l,estimate.phi,verbose)
    return(list(test.table = LRT$test.table,phi = LRT$phi,
                eta=LRT$eta,b=LRT$b,Etheta=LRT$Etheta,converge=LRT$converge))
}


#' Likelihood ratio test
#'
#' @param X the data matrix. The columns must be the markers, the rows must be the samples
#' @param Z the covariates. The length must be the same the the rows of the data matrix X
#' @param tc total number of read counts for each sample
#' @param l marker length
#' @param estimate.phi. Whether to estimate phi (TRUE) or not (FALSE). 
#' @param verbose. Report the estimation details. 
#' @return test.table, eta, b, phi
#' @export
#' @examples
#' likelihood_ratio_test_mada(Data,Covariates,total.counts,marker.length)

likelihood_ratio_test_mada <- function(X,Z,tc,l,estimate.phi=TRUE,verbose=TRUE){
  b.null.index <- rep(FALSE,ncol(Z)+1)
  opt.H1 <- estimate_parameter_mada(X,Z,tc,l,b.null.index,estimate.phi,verbose)
  test.table <- matrix(NA,ncol=2,nrow=length(b.null.index))
  colnames(test.table) <- c('Est.Coeff','Pvalue')
  if (is.null(colnames(Z))){colnames(Z) <- paste('b',1:ncol(Z),sep='')}
  rownames(test.table) <- c('Intersept',colnames(Z))
  
  for (i in 1:nrow(test.table)){
    b.null.index.H0    <- b.null.index
    b.null.index.H0[i] <- TRUE
    test.table[i,'Est.Coeff'] <- opt.H1$b[i]
    opt.H0 <- estimate_parameter_mada(X,Z,tc,l,b.null.index.H0,estimate.phi,verbose=FALSE)
    ## liklihood ratio
    test.stat <- (-2) * (calculate_likelihood_mada(X,Z,tc,l,opt.H0$eta,opt.H0$b,opt.H0$phi)-
                           calculate_likelihood_mada(X,Z,tc,l,opt.H1$eta,opt.H1$b,opt.H1$phi))
    ## print(test.stat)
    ## p value
    ## cat('Test ', names(opt.H0)[i],' D.F. is',length(opt.H1$par)-length(opt.H0[[i]]$par),'\n')
    pvalue <- 1-pchisq(test.stat,df=1)
    test.table[i,'Pvalue'] <- pvalue
  }
  #print(test.table)
  return(list(test.table = test.table,eta=opt.H1$eta,b=opt.H1$b,
              phi=opt.H1$phi,Etheta=opt.H1$Etheta,converge=opt.H1$converge))
}



#' Model-based multi-sample approach for metagenomic differential abundance analysis.
#'
#' Differential abundance analysis for metagenomic data
#' @param X the data matrix. The columns must be the markers, the rows must be the samples
#' @param Z the covariates. The length must be the same the the rows of the data matrix X
#' @param tc total number of read counts for each sample
#' @param l marker length
#' @return test.table 
#' @export
#' @examples
#' calculate_likelihood_mada(X,Z,tc,l,eta,b,phi)

calculate_likelihood_mada <- function(X,Z,tc,l,eta,b,phi){
  ln <- function(X){log(X,base=exp(1))}
  M <- ncol(X)
  N <- nrow(X)
  mu  <- as.matrix(exp(cbind(1,Z) %*% b))
  ####################
  P1 <- lgamma(rowSums(X, na.rm=TRUE)+eta) - lgamma(eta)
  P2 <- rowSums(X*ln(mu %*% t(phi)),na.rm=TRUE) + eta*ln(eta)
  P3 <- (rowSums(X, na.rm=TRUE)+eta)*ln(rowSums((mu*tc) %*% t(phi*l),na.rm=TRUE) + eta)
  logl <- sum(P1 + P2 - P3, na.rm=TRUE)
  return(logl)
}



#' Estimate parameters
#'
#' @param X the data matrix. The columns must be the markers, the rows must be the samples
#' @param Z the covariates. The length must be the same the the rows of the data matrix X
#' @param tc total number of read counts for each sample
#' @param l marker length
#' @param b.null.index vector of TRUE or FALSE indicating which b to estimate (set as TRUE)
#' @param verbose TRUE or FALSE 
#' @return numeric vector giving number of characters in each element of the 
#'   character vector.  Missing strings have missing length.
#' @seealso \code{\link{nchar}} which this function wraps
#' @export
#' @examples
#' EstimateParameter(X,Z,tc,l,b.null.index,verbose=TRUE)

estimate_parameter_mada <- function(X,Z,tc,l,b.null.index,estimate.phi=TRUE,verbose=TRUE){
  ln <- function(X){log(X,base=exp(1))}
  M <- ncol(X)
  N <- nrow(X)
  eta <- 1
  b <- rep(1,length(b.null.index))
  b[b.null.index] <- 0
  b <- as.matrix(b)
  phi <- as.matrix(rep(1,M))
  ita.EM <- 0
  diff <- 100
  converge <- TRUE
  while(ita.EM < 100 & diff > 10^(-4)){
    ita.EM <- ita.EM + 1
    ### Save the value before each iteration
    eta.temp <- eta
    b.temp   <- b
    phi.temp <- phi
    mu  <- as.matrix(exp(cbind(1,Z) %*% b))
    #############
    
    
    ##################
    ### E-step
    E.theta    <- as.matrix((rowSums(X,na.rm=TRUE) + eta)/(eta/mu + tc*sum(phi*l)))
    E.ln.theta <- as.matrix(digamma(rowSums(X,na.rm=TRUE) + eta) - ln(eta/mu + tc*sum(phi*l)))
    #cat('E.theta',E.theta[1],'E.ln.theta',E.ln.theta[1],'\n')
    ##################
    ### M-step
    if (estimate.phi){ 
      phi <- colSums(X,na.rm=TRUE)/(l*sum(E.theta*tc,na.rm=TRUE))
      phi <- phi/sum(phi)*M
      #phi <- phi/sqrt(sum(phi^2)/M)
    }
    else {
      #phi <- 1:M
      phi <- phi/sum(phi)*M
      #phi <- phi/sqrt(sum(phi^2)/M)
    }
    
    QFunction <- function(para,X,Z,l,tc,b.null.index,E.theta,E.ln.theta){
      b   <- rep(NA,length(b.null.index))
      eta <- para[1]
      b[b.null.index]  <- 0
      b[!b.null.index] <- para[-1]
      b   <- as.matrix(b)
      mu  <- as.matrix(exp(cbind(1,Z) %*% b))
      
      
      Q <- sum((eta-1)*E.ln.theta - eta/mu*E.theta-eta*(ln(mu)-ln(eta))-lgamma(eta),na.rm=TRUE)
      return(-Q)
    }
    opt <- nlminb(c(eta,b[!b.null.index]), QFunction, 
                  lower = c(0.001,rep(-Inf,sum(!b.null.index))), upper = Inf, 
                  X = X, Z = Z, l = l, tc = tc,b.null.index=b.null.index,
                  E.theta=E.theta,E.ln.theta=E.ln.theta)
    
    eta <- opt$par[1]
    b[!b.null.index]   <- opt$par[-1]
    if (verbose){cat(ita.EM,opt$convergence,'eta',eta,'b',b,'\n')}
    
    diff <- max(c(abs(eta-eta.temp),abs(b-b.temp),abs(phi-phi.temp)))
  }
  ## convergence:  An integer code. 0 indicates successful convergence
  if (ita.EM>90 & opt$convergence){converge <- FALSE;print('Not converge');print(opt)}
  Etheta <- as.matrix((rowSums(X,na.rm=TRUE) + eta)/(eta/mu + tc*sum(phi*l)))
  return(list(eta=eta,b=b,phi=phi,Etheta=Etheta,opt=opt,converge=converge))
}  
