# Define some variables shared within the package:
package_env <- new.env(parent = emptyenv())
# color palette for heatmaps:
package_env$YlGnBu5   <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58',"#092d94")
package_env$hmcols    <- colorRampPalette(package_env$YlGnBu5)(256)

# helper functions:

#' convert binary codes to decimal integers
#'
#' @param mat A binary matrix the rows of which are the binary codes to be converted.
#' @param LOG Take the logarithm of the decimal integer (default is \code{TRUE});
#' otherwise, set to \code{FALSE}.
#'
#' @return A vector of numbers after conversion
#' @export
#'
#' @examples
#'
#' bin2dec_vec(rbind(c(1,0,0),c(1,1,0),c(0,0,1),c(1,0,1),c(0,0,0)))
#' bin2dec_vec(rbind(c(1,0,0),c(1,1,0),c(0,0,1),c(1,0,1),c(0,0,0)),LOG=FALSE)
bin2dec_vec <- function(mat,LOG=TRUE){# This is the log version:
  L <- ncol(mat)
  M <- nrow(mat)
  v <- log(as.matrix(2^{(L-1):0},ncol=1))
  diag_v <- matrix(0,nrow=L,ncol=L)
  diag(diag_v) <- v
  tmp_prod <- mat%*%diag_v
  permute_M_vec <- sapply(1:M,function(i){matrixStats::logSumExp(tmp_prod[i,mat[i,]!=0])})
  if (!LOG){
    permute_M_vec <- exp(permute_M_vec)
  }
  #permute_M_vec[rev(order(permute_M_vec))]
  permute_M_vec
}


#' Order a binary matrix by row
#'
#' The order is determined by magnitude of the rows as in a binary system (decreasing)
#'
#' @param mat A binary matrix to be ordered by row
#' @param decreasing \code{TRUE} (default) to order rows in a decreasing order;
#' \code{FALSE} otherwise.
#' @return a list of res - the ordered matrix from large to small, ord - the order
#' @examples 
#' themat <- matrix(c(0,1,0,1,0,0,1,1,1,0,1,1),ncol=3)
#' order_mat_byrow(themat)
#' order_mat_byrow(themat,decreasing=FALSE)
#' @export
order_mat_byrow <- function(mat,decreasing=TRUE){
  L <- ncol(mat)
  M <- nrow(mat)
  v <- log(as.matrix(2^{(L-1):0},ncol=1))
  diag_v <- matrix(0,nrow=L,ncol=L)
  diag(diag_v) <- v
  tmp_prod <- mat%*%diag_v
  permute_M_vec <- sapply(1:M,function(i){matrixStats::logSumExp(tmp_prod[i,mat[i,]!=0])})
  if (decreasing){
    return(list(res = mat[rev(order(permute_M_vec)),,drop=FALSE],
                ord = rev(order(permute_M_vec))))
  } else{
    return(list(res = mat[(order(permute_M_vec)),,drop=FALSE],
                ord = (order(permute_M_vec))))
  }
}

#' flip the matrix to the right form for \code{image()}
#'
#' This function makes the result of \code{image()} to be plotted as shown in a matrix
#'
#' @param m a matrix
#' @return a matrix after image() will look the same as the matrix
#' @export
f <- function(m) {t(m)[,nrow(m):1]}


#' create a string with parts that need to be replaced by a variable's actual values
#' 
#' 
#' @param x quote(something)
#' @param nm_in_x the string to in x to be replaced by val (could be a vector)
#' @param val the value to be inserted at "nm_in_x" (could be a vector)
#' @return an expression with 
#' @export
#' 
#' @examples 
#' mytitle <- quote(paste("(",Y[il],"): ",blabla,collapse=""))
#' xx <- "design"
#' eval_string(mytitle,"blabla",xx)
#' 
#' mytitle <- eval_string(quote(paste(a,"--",b,"--","--",c,"--",d)),
#'    c("a","b","c","d"),c("Thanks","for","using","rewind!"))
#' plot(rnorm(10),rnorm(10),main=mytitle)
eval_string <- function(x,nm_in_x,val){
  if (length(nm_in_x)!=length(val)){stop("[rewind] 'nm_in_x' and 'val' are of different lengths.")}
  string_list <- as.list(val)
  names(string_list) <- nm_in_x
  as.expression(do.call('substitute', list(x, string_list)))
}


#' Sample truncated Beta random variables 
#'
#' @param n sample size
#' @param a first Beta parameter
#' @param b second Beta parameter
#' @param lower lower limit
#' @param upper upper limit
#'
#' @return a vector of truncated Beta random variables of size n.
#' @export
#'
#' @examples
#' hist(trunc_rbeta(100000,1,1,0.2,0.8),freq=FALSE,xlim=c(0,1));abline(h=1.67)
#' hist(trunc_rbeta(100000,2,3,0.2,0.8),freq=FALSE,xlim=c(0,1))
trunc_rbeta <- function(n,a,b,lower=0,upper=1){
  u <- runif(n,pbeta(lower,a,b),pbeta(upper,a,b))
  pmax(pmin(0.999999,qbeta(u,a,b)),0.000001)
}




#' Cuthill McKee (CM) algorithm
#' 
#' Transform sparse matrix into a band matrix
#' 
#' @author Fabian Schmich
#' @param x Input matrix
#' @return Band matrix
#' @export
#' 
#' @seealso \code{\link{mblock}}
#' 
#' @examples
#' 
#' L0 <- 100   # dimension of measurements.
#' M0 <- 30    # true dimension of latent states.
#' frac <- 0.05
#' NITER <- 1
#' maxlen_simu <- rep(NA,NITER)
#' for (iter in 1:NITER){
#'   Q    <- simulate_Q(M0,L0,p = frac)
#'   simu <- list(Q=Q)
#'   block_ind_list <- mblock(simu$Q,PLOT = TRUE)
#'   maxlen_simu[iter] <- max(unlist(lapply(block_ind_list,length)))
#' }
#' #hist(maxlen_simu)
#' table(maxlen_simu)
cuthill_mckee <- function(x) {
  degs <- data.frame(Idx=1:ncol(x), NonZero=apply(x, 1, function(x) length(which(x != 0))))
  R <- degs$Idx[which.min(degs$NonZero)]
  i <- 1
  for (i in 1:ncol(x)) {
    Ai <- setdiff(which(x[R[i],] != 0), R)
    if (length(Ai) > 0) {
      Ai <- Ai[order(degs$NonZero[Ai], decreasing = FALSE)]
      R <- append(R, Ai)
    } else {
      R <- append(R, degs$Idx[-R][which.min(degs$NonZero[-R])])
    }
    i <- i + 1
  }
  rR <- rev(R)
  return(list(rR=rR,rcm=x[rR, rR]))
}

#' obtain the indices of rows of Q that form orthogonal blocks
#'
#' @param Q binary Q matrix (M by L), the rows of which are to be split
#' @param PLOT plot blocks in the band matrix or not
#'
#' @return a list of indices that form blocks
#' @export
#' @seealso \code{\link{cuthill_mckee}}
#' 
mblock <- function(Q, PLOT=FALSE){
  res_RCM <- cuthill_mckee(Q%*%t(Q))
  y <- rep(0,M0)
  for (i in 1:M0){
    y[i] <- sum(res_RCM$rcm[1:i,1:i])-sum(diag(res_RCM$rcm)[1:i])
  }
  if (PLOT){
    plot(y)
    par(mfrow=c(2,2))
    image(t(Q),col=package_env$hmcols)
    image(t(Q)[,res_RCM$rR],col=package_env$hmcols)
    image(f(Q%*%t(Q)),col=package_env$hmcols)
    image(f(res_RCM$rcm),col=package_env$hmcols)
  }
  
  nblock <- sum((diff(y)==0))
  
  if (nblock ==0 ){
    return(list(1:nrow(Q)))
  } else{
    res <- matrix(NA,nrow=2,ncol=nblock)
    s <- 1
    res[1,s] <- 1
    res[2,s] <- which(diff(y)==0)[s]
    if (nblock>1){
      for (s in 2:nblock){
        res[1,s] <- res[2,s-1]+1
        res[2,s] <- which(diff(y)==0)[s]
      }
    }
    res_ind <- list()
    for (s in 1:nblock){
      res_ind[[s]] <- res_RCM$rR[res[1,s]:res[2,s]]
    }
    return(res_ind)
  }
}