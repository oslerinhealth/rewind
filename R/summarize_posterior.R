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
#' compute_product_Bern_table(c(0.1,0.2,0.3,0.4))
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
#' @import stats graphics
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
#' @param whoisthis a character string, e.g., subject name.
#' @param ... other parameters for adjusting the individual plots.
#' length to `state_nm`
#'
#' @return a figure with population etiology fraction information (
#' NB: need to add more info.)
#' @export
#' @import stats graphics
plot_individual_pred <- function(p_samp,state_nm,thedata,whoisthis,...){
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
       rep(-(1:M), nrow(pat)),
       axes = FALSE, ann = FALSE,
       pch = ifelse(t(pat), 19, 1), cex = mycex, xpd = NA,lwd=0.2,
       col = rep(c("grey","black")[1+(1:nrow(pat))%in%match(as.numeric(names(p_samp)),round(exp(bin2dec_vec(pat))))],
                 each=M),...)
  # mtext(expression(paste("Individual Etiology Fractions by Infection Patterns (P{",
  #                        eta[i],"=",eta,"|Data})",
  #                        collapse="")),side = 3,line=-3)
  points(rep(-1,M),-(1:M),pch = ifelse(t(matrix(thedata,nrow=1)), 19, 1), cex = mycex,lwd=0.2,
         col="blue")
  axis(2,at=-rev(1:M),
       labels=state_nm,las=2,xpd = NA,cex.axis=mycex)

  points(1:nrow(pat),pat_with_prob[,M+1]*M*multiplier,pch=18,xpd = NA,type="h")
  axis(2,at=M*c(0,0.25,0.5,1)*multiplier,
       labels=c(0,0.25,0.5,1),las=2,xpd = NA)

  thetitle <- bquote(paste(.(whoisthis),": P{",
                                eta[i],"=",eta,"|Data}",
                                collapse=""))
  mtext(thetitle,side = 3,line=0,cex=2)
}


#' Merge latent state profile matrix (H_star) and Q matrix
#'
#' This function
#' \itemize{
#' \item merges identical rows in H: pseudo-cluster ---> scientific clusters;
#' \item merges partner latent states (identical non-zero columns in H);
#' merge the corresponding rows in Q by \code{pmax}.
#' }
#'
#' @param H_star_redun A binary matrix of rows \code{t_max+3} by \code{m_max};
#' It may be redundant because it has rows corresponding to empty pseudo-clusters,
#' and if the entire matrix is non-zero, the zero columns mean inactive latent states
#' @param mylist Cluster labels; A vector of length \code{t_max+3} comprised of
#' integers that need not be consecutive.
#' @param t The number of non-empty pseudo-clusters
#' @param Q The \code{m_max} by \code{L} Q matrix
#' @param VERBOSE Default to \code{FALSE}: no print of the merged matrices; Otherwise
#' set to \code{TRUE}.
#' @param z_pseudo Default is \code{NULL}.
#' A vector of pseudo-cluster indicators to be merged into
#' scientific cluster indicators \code{z_sci}.
#' @param skip_Q Default is \code{FALSE} - merge Q rows with partner states;
#' Set to \code{TRUE} if \code{Q} is given.
#'
#' @return A list comprised of two elements named:
#' \code{H_star_merge} and \code{Q_merge} and \code{z_sci} if \code{z_pseudo}
#' is provided.
#' @export
#'
#' @examples
#'
#'# the third latent state is inactive among non-empty pseudo-clusters;
#'# the 2nd and 4th latent states are partners.
#'# merge 1 and 2 pseudo-clusters.
#'  H_star_redun <- matrix(c(0,1,0,1,
#'                           0,1,0,1,
#'                           1,0,0,0,
#'                           0,0,0,0,
#'                           0,0,0,0,
#'                           0,0,0,0,
#'                           0,0,0,0),nrow=7,ncol=4,byrow=TRUE)
#' mylist <- c(1,2,3,5,0,0,0)
#' t <- 4
#' Q <- simulate_Q(4,100,p=0.1)
#' z_pseudo <- c(1,1,3,5,2,2,2,1,3,3,5,1,2,2,1)
#' merge_H_Q(H_star_redun,mylist,t,Q,TRUE,z_pseudo)
#'
#' # if all zero columns in H_star exists, they will be ignored.
#'
merge_H_Q <- function(H_star_redun,mylist,t,Q,VERBOSE=FALSE,z_pseudo=NULL,skip_Q=FALSE){
  H_star  <- H_star_redun[mylist[1:t],,drop=FALSE] #<-- focus on Eta vectors
  # associated with a pseudo-cluster;
  ind_zero_col <- which(colSums(H_star)==0)
  if (nrow(H_star)>1 && sum(H_star)>0 && length(ind_zero_col)>0){H_star <- H_star[,-ind_zero_col]} #<<<<- removed zero columns here.
  # merge rows (pseudo clusters to scientific clusters defined by \bEta_j):
  pat_H_star    <- apply(H_star,1,paste,collapse="")
  curr_merge    <- merge_map(pat_H_star,unique(pat_H_star))
  #<-- can get the mapping from pseudo clusters to scientific clusters.
  H_star_merge  <- curr_merge$uniq_pat
  if (!is.null(z_pseudo)){
    t_tilde <- nrow(H_star_merge)
    z_sci <- rep(NA,length(z_pseudo))
    sci_lab <- 0
    for (i in 1:t_tilde) {
      sci_lab <- sci_lab+1
      used_cluster_label <- mylist[1:t]
      ind_labels_to_merge <- which(curr_merge$map==i)
      z_sci[z_pseudo%in%used_cluster_label[ind_labels_to_merge]] <- sci_lab
    }
  }
  string_merge1 <- NULL
  if (VERBOSE && nrow(H_star_merge)<nrow(H_star)){
    string_merge1 <- paste0(">> absorbed ",nrow(H_star)-nrow(H_star_merge)," pseudo clusters ---> ", nrow(H_star_merge)," scientific clusters.\n")
    cat(string_merge1)
  }

  # merge columns (combine factors that are present or absent at the same time;
  # partner latent states):
  pat_H_star_merge <- apply(t(H_star_merge),1,paste,collapse="")
  curr_merge_col <- merge_map(pat_H_star_merge,unique(pat_H_star_merge))
  # <-- can get the mapping from partner machines to final merged machines.
  H_star_merge  <- t(curr_merge_col$uniq_pat)
  string_merge2 <- NULL
  if (VERBOSE && ncol(H_star_merge)<ncol(H_star)){
    string_merge2 <- paste0(">> absorbed ",ncol(H_star)-ncol(H_star_merge)+1,"` partner` latent states ---> ", ncol(H_star_merge)," latent states. \n")
    cat(string_merge2)
  }
  if (!skip_Q){
    if (nrow(H_star)>1 && length(ind_zero_col)>0){
      Q_merge <- merge_Q(Q[-ind_zero_col,],curr_merge_col$map)
    } else {
      Q_merge <- merge_Q(Q,curr_merge_col$map)
    }
    if (!is.null(z_pseudo)){
      return(list(H_star_merge=H_star_merge, Q_merge=Q_merge,z_sci=z_sci))
    } else{
      return(list(H_star_merge=H_star_merge, Q_merge=Q_merge))
    }
  } else{
    return(list(H_star_merge=H_star_merge, z_sci=z_sci))
  }
}

#' Merge redundant H_star by column
#'
#' This function removes zero columns if thre are more than two rows in H_star
#' and if \code{H_star_redun} is not trivially all zero
#'
#' @param H_star_redun A binary matrix of rows \code{t_max+3} by \code{m_max};
#' It may be redundant because it has rows corresponding to empty pseudo-clusters,
#' and if the entire matrix is non-zero, the zero columns mean inactive latent states
#' @param VERBOSE Default to \code{FALSE}: no print of the merged matrices; Otherwise
#' set to \code{TRUE}.
#'
#' @return a list of one element named \code{H_star_merge}
#' @export
#'
#' @examples
#'
#'# the third latent state is inactive among non-empty pseudo-clusters;
#'# the 2nd and 4th latent states are partners.
#'# merge 1 and 2 pseudo-clusters.
#'  H_star_redun <- matrix(c(0,1,0,1,
#'                           0,1,0,1,
#'                           1,0,0,0,
#'                           0,0,0,0,
#'                           0,0,0,0,
#'                           0,0,0,0,
#'                           0,0,0,0),nrow=7,ncol=4,byrow=TRUE)
#' Q <- simulate_Q(4,100,p=0.1)
#' merge_H_col(H_star_redun,TRUE)
#'
merge_H_col <- function(H_star_redun,VERBOSE=FALSE){
  #H_star_redun <- out$H_star_samp[,,56]
  H_star  <- H_star_redun #<-- focus on Eta vectors
  # associated with a pseudo-cluster;
  ind_zero_col <- which(colSums(H_star)==0)
  if (nrow(H_star)>1 && sum(H_star)>0 && length(ind_zero_col)>0){H_star <- H_star[,-ind_zero_col]} #<<<<- removed zero columns here.

  # merge columns (combine factors that are present or absent at the same time;
  # partner latent states):
  pat_H_star_merge <- apply(t(H_star),1,paste,collapse="")
  curr_merge_col <- merge_map(pat_H_star_merge,unique(pat_H_star_merge))
  # <-- can get the mapping from partner machines to final merged machines.
  H_star_merge  <- t(curr_merge_col$uniq_pat)
  string_merge2 <- NULL
  if (VERBOSE && ncol(H_star_merge)<ncol(H_star)){
    string_merge2 <- paste0(">> absorbed ",ncol(H_star)-ncol(H_star_merge)+1,"` partner` latent states ---> ", ncol(H_star_merge)," latent states. \n")
    cat(string_merge2)
  }
  # if (nrow(H_star)>1 && length(ind_zero_col)>0){
  #   Q_merge <- merge_Q(Q[-ind_zero_col,],curr_merge_col$map)
  # } else {
  #   Q_merge <- merge_Q(Q,curr_merge_col$map)
  # }

  #return(list(H_star_merge=H_star_merge,Q_merge=Q_merge))
  return(list(H_star_merge=H_star_merge))
}

#' Merge mapping
#'
#' This function combines pseudo-clusters by the current H_star samples.
#'
#' @param pseudo_pat a vector of character strings for
#' binary patterns that might not be unique
#' @param uniq_pat a vector of character strings of distinct binary codes
#'
#' @return a list \itemize{
#' \item \code{map} a vector of integer of identical length to \code{pseudo_pat}; takes values
#' from 1 to \code{length(uniq_pat)}
#' \item \code{uniq_pat} a matrix of unique binary patterns (# rows =  \code{length(uniq_pat)},
#' # columns = number of 1/0s for each element in \code{uniq_pat})
#' }
#' @export
#'
#' @examples
#' #' # simulate data:
#' L0 <- 100
#' options_sim0  <- list(N = 200,  # sample size.
#'                      M = 3,   # true number of machines.
#'                      L = L0,   # number of antibody landmarks.
#'                      K = 8,    # number of true components.,
#'                      theta = rep(0.8,L0), # true positive rates
#'                      psi   = rep(0.01,L0), # false positive rates
#'                      alpha1 = 1 # half of the people have the first machine.
#')
#'
#' #image(simulate_data(options_sim0,SETSEED = TRUE)$datmat)
#'
#' simu     <- simulate_data(options_sim0, SETSEED=TRUE)
#' tmp <-  simu$H_star[c(1,1,2,2,2,3,4,5,6,7,8),]
#' uid <- unique(apply(tmp,1,paste,collapse=""))
#' merge_map(apply(tmp,1,paste,collapse=""),uid)
merge_map <- function(pseudo_pat,uniq_pat){
  if (length(pseudo_pat)<length(uniq_pat)){stop("==[rewind] more pseudo patterns than unique patterns.==")}
  map <- sapply(1:length(pseudo_pat),function(i){which(uniq_pat==pseudo_pat[i])})
  uniq_pat <- do.call("rbind",lapply(sapply(uniq_pat,strsplit,""),as.numeric))
  list(map=map,uniq_pat=uniq_pat)
}


#' merge Q matrix by rows
#'
#' Some rows of Q correspond to partner factors that are present or absent together;
#' It is of scientific interest to combine them by taking the maximum for each column
#' of Q among these rows.
#'
#' @param Q A Q matrix (row for ACTIVE factors that might be partners, columns for dimension of multivariate binary data)
#' @param map_id a vector taking possibly duplicated values in {1,...,M^+}, where M^+ is the number
#' of active factors. \code{map_id=c(1,1,2,2,2,3)} means factor 1 and 2 are partner factors, factor 3 to 5 are another group
#' of partner factors.
#'
#' @return A Q matrix with merged rows (by taking maximum within each group of partner factors)
#' @export
#'
#' @examples
#'
#' Q <- simulate_Q(6,100,0.1)
#' map_id <- c(1,1,2,2,3,2)
#' Q_merge <- merge_Q(Q,map_id)
#' par(mfrow=c(1,2))
#' image(Q,main="before merging")
#' image(Q_merge, main="after merging")
merge_Q <- function(Q,map_id){
  res <- matrix(NA,nrow=length(unique(map_id)),ncol=ncol(Q))
  for (i in unique(map_id)){
    row_id_to_merge <- which(map_id==i)
    res[i,] <- apply(Q[row_id_to_merge,,drop=FALSE],2,max)
  }
  res
}

