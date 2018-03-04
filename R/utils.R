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
#' The order is determined by magnitude of the rows as in a binary system
#'
#' @param mat A binary matrix to be ordered by row
#' @return a list of res - the ordered matrix from large to small, ord - the order
#' @export
order_mat_byrow <- function(mat){
  L <- ncol(mat)
  M <- nrow(mat)
  v <- log(as.matrix(2^{(L-1):0},ncol=1))
  diag_v <- matrix(0,nrow=L,ncol=L)
  diag(diag_v) <- v
  tmp_prod <- mat%*%diag_v
  permute_M_vec <- sapply(1:M,function(i){logSumExp(tmp_prod[i,mat[i,]!=0])})
  list(res = mat[rev(order(permute_M_vec)),,drop=FALSE],
       ord = rev(order(permute_M_vec)))
}
