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
  M0 <- nrow(Q)
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


#' Get frequency count for each integer
#'
#' tabulate v by 1:maxk
#'
#' @param v the vector of integers to be tabulated
#' @param maxk maximum integer to search within v
#'
#' @return a vector of probabilities that sum to one
#' @export
#'
get_freq <- function(v,maxk = 20){
  sapply(1:maxk,function(k) sum(v==k)/length(v))
}







##
##  For slice sampler for infinite dimension models:
##

#' compute the log of density of inactive states' probability given
#' the previous one (Teh et al., 07')
#'
#' This function is used in sampling a model with infinite number of columns in H
#'
#' @param log_p0 log probability (the argument of this function)
#' @param alpha hyperparameter for Indian Buffet Process (Infinite version)
#' @param t the number of rows in the IBP
#'
#' @return a value corresponding to the log density f(log_p0)
#' @examples
#'
#' library(ars)
#' n_grid  <- 10000
#' p0_grid <- seq(log(0.001),0,len=n_grid)
#' y_den  <- log_f_logden(p0_grid,5,8)
#'
#' y <- ars::ars(n_grid,log_f_logden,log_f_logprima,
#' c(-6,-4,-2,-1),m=4,
#' ub=TRUE,xub=0,alpha=5,t=8)
#' hist(y,breaks="Scott",main='Adaptive Rjection Sampling',freq=FALSE)
#' rug(y)
#' points(density(y),type="l",col="blue")
#' points(p0_grid,exp(y_den-rewind:::logsumexp(y_den)-log(diff(p0_grid)[2])),
#'       type="l",col="red") # <-- true density; log_f_logden is only
#'                          #  correct upto a proportionality constant.
#'
#' legend("topleft",c("Sample Density","True Density"),lty=c(1,1),col=c("blue","red"),
#'       bty="n")
#'
#' @export
log_f_logden <- function(log_p0,alpha,t){
  #if (p0 > p0_minus_one || p0 < 0){return(-Inf)}
  res <- 0
  for (j in 1:t){
    res  <- res+(1/j)*(1-exp(log_p0))^j
  }
  alpha*res+(alpha-1)*log_p0+
    t*log(1-exp(log_p0))+log_p0
}


#' compute the derivative of the log of density of inactive states'
#' probability given the previous one (Teh et al., 07')
#'
#' This function is used in sampling a model with an infinite number of columns in H
#'
#' @inheritParams log_f_logden
#'
#' @return a value corresponding to the derivative of the log density d/du f(log_p0)
#' @seealso \code{\link{log_f_logden}} provides an example
#' @export
log_f_logprima <- function(log_p0,alpha,t){
  #if (p0 > p0_minus_one || p0 < 0){return(-Inf)}
  res <- 0
  for (j in 1:t){
    res  <- res - exp(log_p0)*(1-exp(log_p0))^(j-1)
  }
  alpha-t*exp(log_p0)/(1-exp(log_p0))+
    alpha*res
}



#' Metroplis flipping to improve mixing
#'
#' @param x current value, 0 or 1
#' @param curr_prob current probability of being 1
#'
#' @return a binary value
#' @export
#'
#' @examples
#' metrop_flip(1,0.2)
#' 
#' 
#' p = 0.3
#' x <- rep(NA,1000)
#' x[1] <- 0
#' for (iter in 2:1000){
#'   x[iter] <- metrop_flip(x[iter-1],p)
#' }
#' 
#' plot(x,type="l",main=mean(x))
metrop_flip <- function(x,curr_prob){
  if (x==1){
    if (curr_prob <= 0.5) {
      alpha <- 1
    } else{
      alpha <- exp(log(1-curr_prob) - log(curr_prob))
    }
    res <- rbinom(1,1,prob=1-alpha)
  } else{
    if (curr_prob >= 0.5){
      alpha <- 1
    } else{
      alpha <- exp(log(curr_prob) - log(1-curr_prob))
    }
    res <- rbinom(1,1,prob=alpha)
  } #end flipping.
  res
}







#' adjust cluster ids by specification
#'
#' @param cluster_list a list of clusters, each comprised of subject ids known
#' to fall in the same cluster.
#' @param z,N,t,mylist cluster ids for all subjects, cluster sizes, total # clusters, unique cluster ids (includes 0).
#'
#' @return z,N,t,mylist,c_next
#' @export
#' @examples
#' 
#' # simulate data:
#' L0 <- 100   # dimension of measurements.
#' M0 <- 3     # true dimension of latent states.
#' K0 <- 2^M0     # true number of clusters.
#' options_sim0  <- list(N = 100,    # sample size.
#'                       M = M0,     # true number of machines.
#'                       L = L0,     # number of antibody landmarks.
#'                       K = K0,     # number of true components.
#'                       theta = rep(0.8,L0),  # true positive rates.
#'                       psi   = rep(0.15,L0), # false positive rates.
#'                       #alpha1 = 1, # half of the people have the first machine.
#'                       frac = 0.2, # fraction of positive dimensions (L-2M) in Q.
#'                       #pop_frac = rep(1/K0,K0) # population prevalences.
#'                       #pop_frac = (1:K0)/sum(1:K0) # population prevalences.
#'                       pop_frac = c(rep(2,4),rep(1,4)) # population prevalences.
#'                       #pop_frac = c(rep(0.75/4,4),rep(0.25/4,4))
#' )
#' 
#' simu     <- simulate_data(options_sim0, SETSEED=TRUE)
#' 
#' dat <- simu$datmat
#' t_max <- 40
#' n <- options_sim0$N
#' hc <- hclust(dist(dat),"complete")
#' t  <- floor(t_max/4)
#' z  <- cutree(hc,k = t)
#' mylist <- rep(0,t_max+3); mylist[1:t] <- 1:t  # mylist[1:t] is the list of active cluster IDs.
#' c_next <- t+1
#' N <- rep(0,t_max+3); N[1:t] <- table(z)
#' log_p <- rep(0,n+1)
#' zs <- z
#' S  <- rep(0,n)
#' 
#' # before forcing clusters:
#' list(z=z,N=N,t=t,mylist=mylist,c_next=c_next)
#' 
#' cl_list <- list(c(1:(rle(simu$Z)$lengths[1])), 
#'                 c((1+rle(simu$Z)$lengths[1]):(rle(simu$Z)$lengths[2])))
#' cl_list
#' fix_cluster(cl_list,z,N,t,mylist)
fix_cluster <- function(cluster_list,z,N,t,mylist){
  G <- length(cluster_list)
  for (g in 1:G){
    curr_fixed_ind <- cluster_list[[g]]
    diff_clust_id <- setdiff(names(table(z)),names(table(z[-curr_fixed_ind])))
    c_to_reduce <- (mylist[1:t])[match(z[curr_fixed_ind],mylist[1:t])]
    N[as.numeric(names(table(c_to_reduce)))] <- N[as.numeric(names(table(c_to_reduce)))]- table(c_to_reduce)
    #print(length(diff_clust_id))
    if (length(diff_clust_id)==0){
      c_fixed <- ordered_next(mylist)
      z[curr_fixed_ind] <-  c_fixed
      mylist <- ordered_insert(c_fixed,mylist,t)
      t <- t+1
      c_next <- ordered_next(mylist)
      N[c_fixed] <- length(curr_fixed_ind)
    } else{
      for (ss in 1:length(diff_clust_id)){
        # c_to_reduce <- (mylist[1:t])[match(z[curr_fixed_ind],mylist[1:t])]
        # N[as.numeric(names(table(c_to_reduce)))] <- N[as.numeric(names(table(c_to_reduce)))]- table(c_to_reduce)
        mylist <- ordered_remove(as.numeric(diff_clust_id)[ss],mylist,t)
        t <- t-1
        if(ss==1){
          c_fixed <- ordered_next(mylist)
          z[curr_fixed_ind] <- c_fixed
          mylist <- ordered_insert(c_fixed,mylist,t)
          t <- t+1
          N[c_fixed] <- length(curr_fixed_ind)
        }
        c_next <- ordered_next(mylist)
      }
    }
  }
  list(z=z,N=N,t=t,mylist=mylist,c_next=c_next)
}
