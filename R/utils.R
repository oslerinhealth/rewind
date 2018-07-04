# Define some variables shared within the package:
package_env <- new.env(parent = emptyenv())
# color palette for heatmaps:
package_env$YlGnBu5   <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58',"#092d94")
package_env$hmcols    <- colorRampPalette(package_env$YlGnBu5)(256)

# helper functions:

#' convert binary codes to decimal integers
#'
#' @param mat A binary matrix the rows of which are the binary codes to be converted.
#' @param LOG Take the logarithm of the decimal integer (defauly this \code{TRUE});
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
#' @return
#' @export
#'
#' @examples
#' hist(trunc_rbeta(100000,1,1,0.2,0.8),freq=FALSE,xlim=c(0,1));abline(h=1.67)
#' hist(trunc_rbeta(100000,2,3,0.2,0.8),freq=FALSE,xlim=c(0,1))
trunc_rbeta <- function(n,a,b,lower=0,upper=1){
  u <- runif(n,pbeta(lower,a,b),pbeta(upper,a,b))
  pmax(pmin(0.999999,qbeta(u,a,b)),0.000001)
}
