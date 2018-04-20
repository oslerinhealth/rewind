#' Compute the posterior similarity matrix (probability that i and j are in same cluster).
#'
#' This function is to evaluate the recovery of clusters
#' @param n total number of subjejcts
#' @param z a matrix of posterior samples, rows for subejcts, columns for MCMC samples kept across
#' multiple iterations
#' @param n_keep the number of columns in z
#'
#' @return a matrix of co-clustering frequencies
#' @export
coclust_mat <- function(n,z,n_keep){
  res     <- matrix(0,n,n)
  count <- 0
  for (index in 1:n_keep){
    for (i in 1:n){
      for (j in 1:n){
        if (z[i,index]==z[j,index]){
          res[i,j] <- res[i,j]+1
        }
      }
    }
    count <- count+1
  }
  res/count
}

#' Frobenius norm: ||AA'-BB'||_F
#'
#'
#' @param A,B two matrixes of the same number of rows
#'
#' @return a positive number
#' @export
#'
distAB <- function(A,B){
  norm(A%*%t(A) - B%*%t(B),"F")
}

#' Convert the cluster indicators to the co-activation patterns.
#'
#' @param z a vector of integers indicating cluster assignments for all subjects
#'
#' @return an N by N binary; an 1 represents that two subjects are in the same
#' cluster; 0 otherwise.
#' @export
#'
#' @examples
#' z2comat(c(1,1,2,2,3,4,5,6,5,7))
z2comat <- function(z){ # for Dahl (2006): simple least square approach to estimate clusters.
  n <- length(z)
  res <- matrix(NA,nrow=n,ncol=n)
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      res[i,j] <- res[j,i] <- (z[i]==z[j])+0
    }
  }
  diag(res) <- 1
  res
}


#' Compute the product Bernoulli probability table
#'
#' @param p A vector of M probabilities
#' @param LOG Default is TRUE; logarithm of probabilities.
#'
#' @return a vector of probabilities that sum to one
#' @export
#'
#' @examples
#'
#' compute_table(c(0.1,0.2,0.3,0.4))
compute_product_Bern_table <- function(p,LOG=TRUE){
  M <- length(p)
  pat <- as.matrix(expand.grid(rep(list(0:1),M)),ncol=M)
  pat <- pat[order(rowSums(pat)),]
  res <- rowSums(pat%*%diag(log(p))+(1-pat)%*%diag(log(1-p)))
  if (!LOG){res <- exp(res)}
  names(res) <- apply(pat,1,paste,collapse="")
  res
}



#' Plot population "etiology" fractions
#'
#' @param p_samp A matrix of posterior samples of Bernoulli probabilities (
#' rows for latent state dimensions, columns for MCMC iterations.)
#' @param state_nm A string vector for the names of the latent state dimensions
#'
#' @return a figure with population etiology fraction information (
#' NB: need to add more info.)
#' @export
#' @import stats
plot_population_fractions <- function(p_samp,state_nm){
  M <- nrow(p_samp)
  log_eti_samp  <- apply(p_samp,2,compute_product_Bern_table)
  eti_samp_mean <- rowMeans(exp(log_eti_samp))

  pat <-  as.matrix(expand.grid(rep(list(0:1),M)),ncol=M)
  pat <-  pat[order(rowSums(pat)),]

  multiplier <- 2
  mycex      <- 1.5
  plot(rep(1:nrow(pat), each = ncol(pat)),
       rep(-rev(1:ncol(pat)), nrow(pat)),
       axes = FALSE, ann = FALSE,
       pch = ifelse(t(pat), 19, 1), cex = mycex,asp=2, xpd = NA,lwd=0.2,
       col = rep(c("grey","black")[1+(eti_samp_mean>quantile(eti_samp_mean,0.75))],
                 each=ncol(pat)))
  mtext(expression(paste("Population Etiology Fraction by Infection Patterns (",
                         pi[eta],")",
                         collapse="")),side = 3,line=-3)
  axis(2,at=-rev(1:ncol(pat)),
       labels=state_nm,las=2,xpd = NA,cex.axis=mycex)

  segments(1:nrow(pat),exp(apply(log_eti_samp,1,quantile,0.025))*ncol(pat)*multiplier,
           1:nrow(pat),exp(apply(log_eti_samp,1,quantile,0.975))*ncol(pat)*multiplier)
  segments(1:nrow(pat),exp(apply(log_eti_samp,1,quantile,0.25))*ncol(pat)*multiplier,
           1:nrow(pat),exp(apply(log_eti_samp,1,quantile,0.75))*ncol(pat)*multiplier,lwd=4,col="dodgerblue2")
  points(1:nrow(pat),rowMeans(exp(log_eti_samp))*ncol(pat)*multiplier,pch=18,xpd = NA)
  points(1:nrow(pat),exp(apply(log_eti_samp,1,quantile,0.025))*ncol(pat)*multiplier,lwd=6,pch="-")
  points(1:nrow(pat),exp(apply(log_eti_samp,1,quantile,0.975))*ncol(pat)*multiplier,lwd=6,pch="-")
  axis(2,at=ncol(pat)*c(0,0.25,0.5,1)*multiplier,
       labels=c(0,0.25,0.5,1),las=2,xpd = NA)

  mtext(text = expression(paste("", pi[eta],sep="")),line = 3,side=2,
        cex=2,las=2,padj=-8)
  #mtext(text = expression(paste(eta[1],",..., ", eta[M],collapse="")),line = 1,side=2,
  #      cex=1,adj = 0.5,las=0)
  legend("topright",legend = c("95% CI","50% CI"),col=c("black","dodgerblue2"),lwd=c(1,4),
         bty="n")

  # # This is the best estimated Q:
  # #
  # NROW_Q_PLOT <- ncol(dat_rlcm_case) #sum(rowSums(Q_merged)!=0)
  # Q_PLOT <- diag(NROW_Q_PLOT)
  # #<- f(order_mat_byrow(Q_merged[rowSums(Q_merged)!=0,,drop=FALSE])$res)
  # image(seq(1,nrow(pat),length=ncol(dat_rlcm_case)),-rev(1:NROW_Q_PLOT)-NROW_Q_PLOT-1,
  #       Q_PLOT,
  #       main="Best Q (merged & ordered)",col=hmcols,
  #       xlab="Dimension (1:L)",
  #       ylab="Latent State (1:M)",yaxt="n",cex.lab=1.2,add=TRUE)
  # for (k in 1:(NROW_Q_PLOT+1)){
  #   abline(h=-k+0.5-NROW_Q_PLOT-1,
  #          lty=2,col="grey",lwd=2)
  # }
  #
  # mtext(text = expression(paste(Q[1.],",..., ", Q[M.],collapse="")),line = 1,side=2,
  #       cex=1,adj = 0.15,las=0)

}



#' Plot individual probabilities
#'
#' @param p_samp A vector of posterior probabilities for each individual
#' @param state_nm A string vector for the names of the latent state dimensions
#' @param thedata a vector of binary measurements for a subject; of identical
#' length to `state_nm`
#'
#' @return a figure with population etiology fraction information (
#' NB: need to add more info.)
#' @export
#' @import stats
plot_individual_pred <- function(p_samp,state_nm,thedata){
  # p_samp = apply(H_pat_res,1,table)[[i]]/sum(apply(H_pat_res,1,table)[[i]])
  # state_nm = analysis_list
  # thedata = dat_rlcm_case[i,,drop=FALSE]

  M   <- length(state_nm)
  pat <-  as.matrix(expand.grid(rep(list(0:1),M)),ncol=M)
  pat <-  pat[order(rowSums(pat)),]
  pat_with_prob <- cbind(pat,rep(0,nrow(pat)))
  pat_with_prob[match(as.numeric(names(p_samp)),round(exp(bin2dec_vec(pat)))),M+1] <- p_samp

  multiplier <- 2
  mycex      <- 1.5
  plot(rep(1:nrow(pat), each = M),
       rep(-rev(1:M), nrow(pat)),
       axes = FALSE, ann = FALSE,
       pch = ifelse(t(pat), 19, 1), cex = mycex,asp=2, xpd = NA,lwd=0.2,
       col = rep(c("grey","black")[1+(1:nrow(pat))%in%match(as.numeric(names(p_samp)),round(exp(bin2dec_vec(pat))))],
                 each=M))
  # mtext(expression(paste("Individual Etiology Fractions by Infection Patterns (P{",
  #                        eta[i],"=",eta,"|Data})",
  #                        collapse="")),side = 3,line=-3)
  points(rep(-1,M),-rev(1:M),pch = ifelse(t(matrix(thedata,nrow=1)), 19, 1), cex = mycex,lwd=0.2,
         col="blue")
  axis(2,at=-rev(1:M),
       labels=state_nm,las=2,xpd = NA,cex.axis=mycex)

  points(1:nrow(pat),pat_with_prob[,M+1]*M*multiplier,pch=18,xpd = NA,type="h")
  axis(2,at=M*c(0,0.25,0.5,1)*multiplier,
       labels=c(0,0.25,0.5,1),las=2,xpd = NA)

  mtext(expression(paste("P{",
                         eta[i],"=",eta,"|Data}",
                         collapse="")),side = 3,line=0,cex=2)
  legend("topright",legend = c("95% CI","50% CI"),col=c("black","dodgerblue2"),lwd=c(1,4),
         bty="n")
}
