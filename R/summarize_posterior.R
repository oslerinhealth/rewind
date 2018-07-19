###### Compute the posterior co-clustering probability matrix (probability that i and j 
###### are clustered together).
######
###### This function is to evaluate the recovered clusters
###### 
###### @param z a matrix of posterior samples, with subjects and MCMC samples in 
###### the rows and columns, respectively.
######
###### @return a matrix of empirical co-clustering frequencies based on the
###### posterior samples
###### 
###### @examples 
###### z2comat(matrix(c(1,1,2,2,3,4,5,6,5,7),ncol=1))
###### 
###### @export
#####z2comat <- function(z){
#####  if(!is.matrix(z)){z <- matrix(z,ncol=1)}
#####  n_inference <- ncol(z)
#####  n <- nrow(z)
#####  res     <- matrix(0,n,n)
#####  count <- 0
#####  for (index in 1:n_inference){
#####    for (i in 1:n){
#####      for (j in 1:n){
#####        if (z[i,index]==z[j,index]){
#####          res[i,j] <- res[i,j]+1
#####        }
#####      }
#####    }
#####    count <- count+1
#####  }
#####  res/count
#####}

#' Frobenius norm: ||AA'-BB'||_F
#'
#' @param A,B two matrixes of the same number of rows
#'
#' @return a positive number
#' @export
distAB <- function(A,B){
  norm(A%*%t(A) - B%*%t(B),"F")
}

#' Compute the product Bernoulli probability table
#'
#' This function maybe useful for marginalization or just to visualize the
#' probability table.
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

#' Plot population "etiology" fractions (for known Q only)
#'
#' @param p_samp A matrix of posterior samples of Bernoulli probabilities (
#' rows for latent state dimensions, columns for MCMC iterations.)
#' @param state_nm A string vector for the names of the latent state dimensions
#' @param mycex The cex parameter for the size of the latent state names
#' plotted on the left margin; Default is 1.
#' @param ... other parameters for adjusting the individual plots.
#'
#' @return a figure with population etiology fraction information (
#' NB: need to add more info.)
#' @export
#' @import stats graphics
plot_population_fractions <- function(p_samp,state_nm,mycex=1,...){
  M <- nrow(p_samp)
  log_eti_samp  <- apply(p_samp,2,compute_product_Bern_table)
  eti_samp_mean <- rowMeans(exp(log_eti_samp))
  
  pat <-  as.matrix(expand.grid(rep(list(0:1),M)),ncol=M)
  pat <-  pat[order(rowSums(pat)),]
  
  multiplier <- 2
  plot(rep(1:nrow(pat), each = ncol(pat)),
       rep(-(1:ncol(pat)), nrow(pat)),
       axes = FALSE, ann = FALSE,
       pch = ifelse(t(pat), 19, 1), cex = mycex,xpd = NA,lwd=0.2,
       col = rep(c("grey","black")[1+(eti_samp_mean>quantile(eti_samp_mean,0.75))],
                 each=ncol(pat)),...)
  mtext(expression(paste("Population Etiology Fraction by Infection Patterns (",
                         pi[eta],")",
                         collapse="")),side = 3,line=-3)
  axis(2,at=-(1:ncol(pat)),
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
  #       main="Best Q (merged & ordered)",col=package_env$hmcols,
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
#' @param whoisthis a character string, e.g., subject name.
#' @param thedata a vector of binary measurements for a subject; of identical
#' @param mycex the cex parameter for the names of the latent states. Default is 1.
#' @param ... other parameters for adjusting the individual plots.
#' length to `state_nm`
#'
#' @return a figure with population etiology fraction information (
#' NB: need to add more info.)
#' @export
#' @import stats graphics
plot_individual_pred <- function(p_samp,state_nm,whoisthis,thedata=NULL,mycex=1,...){
  # p_samp = apply(H_pat_res,1,table)[[i]]/sum(apply(H_pat_res,1,table)[[i]])
  # state_nm = analysis_list
  # thedata = dat_rlcm_case[i,,drop=FALSE]
  
  M   <- length(state_nm)
  pat <-  as.matrix(expand.grid(rep(list(0:1),M)),ncol=M)
  pat <-  pat[order(rowSums(pat)),]
  pat_with_prob <- cbind(pat,rep(0,nrow(pat)))
  pat_with_prob[match(as.numeric(names(p_samp)),round(exp(bin2dec_vec(pat)))),M+1] <- p_samp
  
  multiplier <- 2
  plot(rep(1:nrow(pat), each = M),
       rep(-(1:M), nrow(pat)),
       axes = FALSE, ann = FALSE,
       pch = ifelse(t(pat), 19, 1), cex = mycex, xpd = NA,lwd=0.2,
       col = rep(c("grey","black")[1+(1:nrow(pat))%in%match(as.numeric(names(p_samp)),round(exp(bin2dec_vec(pat))))],
                 each=M),...)
  # mtext(expression(paste("Individual Etiology Fractions by Infection Patterns (P{",
  #                        eta[i],"=",eta,"|Data})",
  #                        collapse="")),side = 3,line=-3)
  if (!is.null(thedata)){
    if (length(thedata)!=M){
      stop("[rewind] observed data and latent states have unequal dimensions.Set thedata to NULL and retry.")
    }
    points(rep(0,ncol(thedata)),-(1:ncol(thedata)),pch = ifelse(t(matrix(thedata,nrow=1)), 15, 0), cex = mycex,lwd=0.5,
           col="blue",xpd=NA)
    arrows(0,-M-2,0,-M-0.5,lwd=3,angle=20,length=0.1,col="blue")
    text(0,-M-2.5,"observed",pos=4,col="blue",cex=1.5)
  }
  axis(2,at=-(1:M),
       labels=state_nm,las=2,xpd = NA,cex.axis=mycex)
  
  points(1:nrow(pat),pat_with_prob[,M+1]*M*multiplier,pch=18,xpd = NA,type="h",lwd=2)
  axis(2,at=M*c(0,0.25,0.5,1)*multiplier,
       labels=c(0,0.25,0.5,1),las=2,xpd = NA)
  
  thetitle <- bquote(paste(.(whoisthis),": P{",
                           eta[i],"=",eta,"|Data}",
                           collapse=""))
  mtext(thetitle,side = 3,line=0,cex=2)
}

#' Plot posterior sampling chain and the histogram on the right margin 
#' 
#' NB: currently works for integer y; need to make it general.
#'
#' @param x numeric vector; usually iteration ids
#' @param y the values at each x value
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#'
#' @return a plot with sampling chain at the bottom left and histogram
#' on the right margin
#' @export
#'
#' @examples
#' x <- 1:1000
#' y <- sample(1:10,1000,prob = c(1,2,3,4,5,5,4,3,2,1), replace=TRUE)
#' scatterhist(x,y)
#' 
scatterhist = function(x, y, xlab="", ylab=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  plot(x,y,type="l",cex.lab=1.5)
  par(mar=c(0,3,1,1))
  #barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  plot.new()
  par(mar=c(3,0,1,1))
  barplot(table(y) , axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0,
        at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=2, outer=TRUE, adj=0,
        at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}


#' Merge latent state profile matrix (H_star) and Q matrix
#'
#' This function
#' \itemize{
#' \item merges identical rows in H: pseudo-cluster ---> scientific clusters;
#' \item merges partner latent states (identical non-zero columns in H);
#' merge the corresponding rows in Q by \code{pmax}. Even though identifiability
#' condition may prevent non-identical columns in H, because we specify
#' a maximum latent state dimension \code{m_max}, MCMC may produce identical 
#' columns in some iterations.
#' }
#'
#' @param H_star_redun A binary matrix of rows \code{t_max+3} by \code{m_max};
#' It may be redundant because it has rows corresponding to empty pseudo-clusters,
#' and if the matrix has at least one non-zero element, the zero columns indicate 
#' inactive latent states
#' @param mylist Cluster labels; A vector of length \code{t_max+3} comprised of
#' integers that need not be consecutive.
#' @param t The number of non-empty pseudo-clusters
#' @param Q The \code{m_max} by \code{L} Q matrix
#' @param VERBOSE Default to \code{FALSE}: no print of the merged matrices; Otherwise
#' set to \code{TRUE}.
#' @param z_pseudo Default is \code{NULL}.
#' A vector of pseudo-cluster indicators to be merged into scientific cluster 
#' indicators \code{z_sci}.
#' @param skip_Q Default is \code{FALSE} - merge the rows of Q that
#' correspond to partner states; Set to \code{TRUE} if \code{Q} is given 
#' in which case only rows of \code{H_star_redun} will be merged (keep
#' the zero columns in H_star_redun, and not merging partner states).
#'
#' @return A list comprised of two elements named:
#' \code{H_star_merge} and \code{Q_merge} and \code{z_sci} if \code{z_pseudo}
#' is provided.
#' @export
#'
#' @examples
#'# the third latent state is inactive among non-empty pseudo-clusters;
#'# the 2nd and 4th latent states are partners.
#'# merge 1st and 2nd pseudo-clusters.
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
#' # Q is unknown and sampled
#' merge_H_Q(H_star_redun,mylist,t,Q,TRUE,z_pseudo)
#' # Q is known:
#' merge_H_Q(H_star_redun,mylist,t,Q,TRUE,z_pseudo,TRUE)
merge_H_Q <- function(H_star_redun,mylist,t,Q,VERBOSE=FALSE,z_pseudo=NULL,skip_Q=FALSE){
  H_star  <- H_star_redun[mylist[1:t],,drop=FALSE] #<-- focus on Eta vectors
  # associated with a pseudo-cluster;
  # if all zero columns in H_star exists, they will be ignored:
  ind_zero_col <- which(colSums(H_star)==0)
  if (!skip_Q){
    if (nrow(H_star)>1 && sum(H_star)>0 && length(ind_zero_col)>0){
      H_star <- H_star[,-ind_zero_col,drop=FALSE];
      string_removed_zero_column <-  
        paste0(">> ignored ",length(ind_zero_col),
               " `non-active` latent states ---> ", ncol(H_star),
               " latent states. \n")
      if (VERBOSE){
        cat(string_removed_zero_column)
      }
    } #<- removed zero columns here.
  }
  # merge rows (pseudo clusters to scientific clusters defined by \bEta_j):
  pat_H_star    <- apply(H_star,1,paste,collapse="")
  curr_merge    <- merge_map(pat_H_star,unique(pat_H_star))
  #<-- get the mapping from pseudo clusters to scientific clusters.
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
    string_merge1 <- paste0(">> absorbed ",
                            nrow(H_star)-nrow(H_star_merge),
                            " pseudo clusters ---> ", 
                            nrow(H_star_merge)," scientific clusters.\n")
    cat(string_merge1)
  }
  
  if (skip_Q){
    H_star_row_merge <- H_star_merge
  }
  
  # merge columns (combine factors that are present or absent at the same time;
  # partner latent states):
  pat_H_star_merge <- apply(t(H_star_merge),1,paste,collapse="")
  curr_merge_col <- merge_map(pat_H_star_merge,unique(pat_H_star_merge))
  # <-- can get the mapping from partner machines to final merged machines.
  H_star_merge  <- t(curr_merge_col$uniq_pat)
  string_merge2 <- NULL
  if (VERBOSE && ncol(H_star_merge)<ncol(H_star)){
    string_merge2 <- paste0(">> absorbed ",
                            ncol(H_star)-ncol(H_star_merge),
                            " `partner` latent states ---> ", 
                            ncol(H_star_merge)," latent states. \n")
    cat(string_merge2)
  }
  if (!skip_Q){
    if (nrow(H_star)>1 && length(ind_zero_col)>0){
      Q_merge <- merge_Q(Q[-ind_zero_col,,drop=FALSE],curr_merge_col$map)
    } else {
      Q_merge <- merge_Q(Q,curr_merge_col$map)
    }
    if (!is.null(z_pseudo)){
      return(list(H_star_merge=H_star_merge, Q_merge=Q_merge,z_sci=z_sci))
    } else{
      return(list(H_star_merge=H_star_merge, Q_merge=Q_merge))
    }
  } else{
    if (!is.null(z_pseudo)){
      return(list(H_star_merge=H_star_row_merge, Q_merge=Q,z_sci=z_sci))
    } else{
      return(list(H_star_merge=H_star_row_merge, Q_merge=Q))
    }
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
    string_merge2 <- paste0(">> absorbed ",ncol(H_star)-ncol(H_star_merge)+1," `partner` latent states ---> ", ncol(H_star_merge)," latent states. \n")
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
#'                       M = 3,   # true number of machines.
#'                       L = L0,   # number of antibody landmarks.
#'                       K = 8,    # number of true components.,
#'                      theta = rep(0.8,L0), # true positive rates
#'                      psi   = rep(0.01,L0), # false positive rates
#'                      alpha1 = 1, # half of the people have the first machine.
#'                      frac = 0.2, # fraction of positive dimensions (L-2M) in Q.
#'                      #pop_frac = rep(1/K0,K0) # population prevalences.
#'                      #pop_frac = (1:K0)/sum(1:K0) # population prevalences.
#'                      pop_frac = c(rep(2,4),rep(1,4)) # population prevalences.
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
#' @param Q A Q matrix (row for ACTIVE factors that might be partners, 
#' columns for dimension of multivariate binary data)
#' @param map_id a vector taking possibly duplicated values in {1,...,M^+}, where M^+ is the number
#' of active factors. \code{map_id=c(1,1,2,2,2,3)} means factor 1 and 2 are partner factors, factor 3 to 5 are another group
#' of partner factors.
#'
#' @return A Q matrix with merged rows (by taking maximum within each group of partner factors)
#' NB: does this work for other restricted LCM models? (DINO - ok; what
#' about two latent states that never co-exist? DINA - okay by duality. For
#' general RLCM, the log-linear representation will only have q_jk*q_jk',
#' after merging, the term will be q_jk'', where k'' is merged version of k and k'.
#' The population model will be find because if there is no person with k or k' only,
#' the parameter for q_jk or q_jk' alone could not be estimated - no differentials
#' in response probability levels.)
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
  #print(map_id)
  #print(Q)
  res <- matrix(NA,nrow=length(unique(map_id)),ncol=ncol(Q))
  for (i in unique(map_id)){
    row_id_to_merge <- which(map_id==i)
    res[i,] <- apply(Q[row_id_to_merge,,drop=FALSE],2,max)
  }
  res
}

#' Merge H and Q
#'
#' NB: need to add another option to skip merge Q when Q is known, in
#' which case, no col merging, no Q merging, just simple row-wise merging.
#' Currently function merge_H_Q can deal with this. So just need to use that 
#' functionality by specifying skip_Q = TRUE.
#'
#' @param outres output from \code{\link{sampler}}
#' @param Q Q can be specified to be known. Default is \code{NULL}.
#'
#' @return a new output with extra elements:
#' \code{H_star_merge_samp}, \code{Q_merge_samp},
#' \code{z_sci_samp}, \code{col_merged_H_star_samp}.
#' @export
postprocess_H_Q <- function(outres,Q=NULL){
  n_kept <- dim(outres$H_star_samp)[3]
  curr_m_max <- ncol(outres$H_star_samp)
  curr_L <- dim(outres$theta_samp)[1]
  Q_merge_samp <- array(0,c(curr_m_max,curr_L,n_kept))
  H_star_merge_samp <- array(0,c(dim(outres$H_star_samp)[1],curr_m_max,n_kept))
  col_merged_H_star_samp <- H_star_merge_samp
  z_sci_samp <- matrix(0,nrow=dim(outres$z_samp)[1],ncol=n_kept)
  for (iter in 1:n_kept){
    if (is.null(Q)){
      merged_res <- merge_H_Q(outres$H_star_samp[,,iter],
                              outres$mylist[,iter],
                              outres$t_samp[iter],
                              outres$Q_samp[,,iter],
                              VERBOSE=FALSE,
                              outres$z_samp[,iter],FALSE)
    } else{
      merged_res <- merge_H_Q(outres$H_star_samp[,,iter],
                              outres$mylist[,iter],
                              outres$t_samp[iter],
                              Q,
                              VERBOSE=FALSE,
                              outres$z_samp[,iter],TRUE)
    }
    
    merged_H <- merged_res$H_star_merge
    merged_Q <- merged_res$Q_merge
    H_star_merge_samp[1:nrow(merged_H),1:ncol(merged_H),iter] <- merged_H
    Q_merge_samp[1:nrow(merged_Q),1:ncol(merged_Q),iter] <- merged_Q
    z_sci_samp[,iter]    <- merged_res$z_sci
    if (is.null(Q)){
      col_merged_H_star_samp[,1:ncol(merged_H),iter] <- merge_H_col(
        outres$H_star_samp[,,iter])$H_star_merge # <-- may still have extra zeros.
    } else{
      col_merged_H_star_samp[,1:ncol(merged_H),iter] <- 
        outres$H_star_samp[,,iter]
    }
  }
  outres$H_star_merge_samp <- H_star_merge_samp # <-- could have zero columns.
  outres$Q_merge_samp      <- Q_merge_samp # <-- could have zero rows.
  outres$z_sci_samp        <- z_sci_samp # <-- add scientific indicators.
  outres$col_merged_H_star_samp <- col_merged_H_star_samp # <-- for visualizing marginal probabilities of latent states.
  outres
}


#' Summarize individual latent state patterns given Q
#'
#' This function takes in posterior samples of individual latent states,
#' orders latent states according to the rows of Q and output 
#' \itemize{
#' \item 1) the
#' marginal probabilities of each individual's latent states being active
#' and
#' \item 2) the posterior samples of individual latent states (after ordering
#' by row of Q and converting the binary codes to decimal numbers)
#' }
#'
#' @param H_star_samp A latent state profile array with elements indexed by
#' (pseudo-cluster, latent state, MCMC iteration for inference)
#'
#' @param z_samp A pseudo-cluster indicator matrix with elements indexed by (individual, MCMC
#' iteration for inference)
#' @param Q Q matrix; known but could have rows to be ordered
#'
#' @return
#' \itemize{
#' \item \code{marg_prob} the
#' marginal probabilities of each individual's latent states being active
#' \item \code{H_pat_ordered_samp} the posterior samples of individual
#'  latent states (after ordering
#' by row of Q and converting the binary codes to decimal numbers)
#' }
#' @export
summarize_latent_state_given_Q <- function(H_star_samp, z_samp, Q){
  if (sum(rowSums(Q)==0)>0){stop("[rewind] Some rows of Q are all zeros. Please remove these rows in Q and corresponding
                                 columns in H_star_samp and retry.")}
  n     <- dim(z_samp)[1]
  M     <- nrow(Q)
  n_inf <- dim(H_star_samp)[3]
  H_res <- matrix(0,nrow=n,ncol=M)
  H_pat_res <- matrix(0,nrow=n,ncol=n_inf)
  for (iter in 1:n_inf){
    tmp_not_ordered <- H_star_samp[z_samp[,iter],,iter]
    tmp_ordered <- tmp_not_ordered[,order_mat_byrow(Q)$ord,drop=FALSE]
    H_pat_res[,iter] <- round(bin2dec_vec(tmp_ordered,LOG=FALSE))
    H_res <- (tmp_ordered+H_res*(iter-1))/iter
  }
  list(marg_prob = H_res,H_pat_ordered_samp = H_pat_res)
}



#' Plot data, design matrix or co-clustering probabilities
#' 
#'
#' @param x data, design matrix or coclustering probabilities
#' @param simu simulated objects; see \code{\link{simulate_data}}
#' @param type default to "data"; could be "design" or "cocluster"
#' @param title_val the title of the plot; default is just \code{type}
#'
#' @return a plot
#' @export
#'
#' @examples
#' # simulate data:
#' L0 <- 100   # dimension of measurements.
#' M0 <- 3     # true dimension of latent states.
#' K0 <- 8     # true number of clusters.
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
#' plot_rewind(simu$datmat,simu,"data","Simulated Data")
plot_rewind <- function(x,simu,type="data",title_val=type){
  L <- ncol(x)
  N <- nrow(x)
  M <- ncol(simu$H_star)
  K <- nrow(simu$H_star)
  if (type=="data" | type =="design"){
    image(1:L,1:N,
          f(x),
          col=package_env$hmcols,
          xlab="Dimension (1:L)",yaxt="n",
          ylab="Subject (1:N)",cex.lab=1.2)
    axis(side=2,at=seq(0,N,by=10)+1,
         labels=c(rev(seq(0,N,by=10)[-1]),"1"),las=2,cex.axis=1.2)
    if (type=="data"){mytitle <- quote(paste("(",Y[il],"): ",tt,collapse=""))}
    if (type=="design"){mytitle <- quote(paste("(",Gamma[il],"): ",tt,collapse=""))}
    mtext(eval_string(mytitle,"tt",title_val),3,1)
    for (k in 1:K){
      abline(h=N-cumsum(rle(simu$Z)$lengths)[k]+0.5,
             lty=2,col="grey",lwd=2)
    }
  } else if (type =="cocluster"){
    image(1:N,1:N,
          f(x),col=package_env$hmcols,
          main=title_val,
          xlab="Subject (1:N)",
          ylab="Subject (1:N)",cex.lab=1.2,yaxt="n")
    axis(side=2,at=seq(0,N,by=10)+1,
         labels=c(rev(seq(0,N,by=10)[-1]),"1"),las=2,cex=1.2)
    for (k in 1:K){
      abline(h=cumsum(rle(rev(simu$Z))$lengths)[K-k+1]+0.5,lty=2)
      abline(v=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
    }
  } else{
    stop("[rewind] Requested 'type' not available. Please reselect.")
  }
}


