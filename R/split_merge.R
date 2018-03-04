# Define helper functions:

#' insert a cluster ID prior to "index"; mylist keeps a list of cluster IDs
#'
#' This function is called when adding a new cluster with more than one subjects
#'
#' @param index an integer of cluster ID to be inserted
#' @param mylist a vector of integer cluster IDs
#' @param t the effective number of clusters
#'
#' @examples
#' ordered_insert(4,c(1:3,5:10),7)
#' ordered_insert(4,c(1:10),7)
#'
#' @return an vector of cluster IDs with an inserted ID \code{index}.
#' @export
ordered_insert <- function(index,mylist,t){
  j <- t
  while (j>0 && mylist[j]>index){
    mylist[j+1] <- mylist[j]
    j <- j-1
  }
  mylist[j+1] <- index # all elements greater than index have been moved one position to the right,
  # leaving space for index.
  mylist
}


#' remove an index from a vector of cluster ID
#'
#' This function is called when removing a cluster.
#'
#' @inheritParams ordered_insert
#'
#' @examples
#' ordered_remove(3,c(1:3,5:10),7)
#' ordered_remove(3,c(1:10),7)
#'
#'
#'# some quick examples:
#'
#' list0 <- c(1:7,0,0,0)
#' k <- 1
#' while (list0[k]==k){ k=k+ 1}; cm = k
#' while (list0[k]==k+1){ k =k+ 1}; ci = k+1
#' while (list0[k]==k+2) {k =k+ 1};  cj = k+2
#'
#' a1 <- ordered_remove(3,list0,7)
#' a2 <- ordered_remove(5,a1,6)
#' a3 <- ordered_insert(cm,a2,5)
#'
#' # end of examples.
#'
#' @export
ordered_remove <- function(index,mylist,t){
  for (j in 1:t){ # it stopped at t, which might be smaller than length of mylist.
    if (mylist[j] >= index){
      mylist[j] <- mylist[j+1]
    }
  }
  mylist # only t-1 matters now.
}


#' Find the next available cluster ID to use:
#'
#' @param mylist a vector of integer cluster IDs
#' @examples
#' ordered_next(c(1:3,5:10))
#' ordered_next(c(1:3,5:6,9:10))
#' @return an index to use next for an observation during complete Gibbs scan
#' @export
ordered_next <- function(mylist){
  j <- 1
  while (mylist[j]==j) {
    j <- j+1
  }
  j # identify the first gap.
}

#' Restricted Gibbs scan to update assignment indicators in S.
#'
#' Function to perform Gibbs scan only for the observations associated
#' with either i or j (randomly chosen), but not including the cluster indicator
#' for i or j. This function does not change the cluster observation i or j currently is in,
#' because the restricted Gibbs scan only updates other observations in S by randomly
#' assigning them to the cluster with i or j.
#'
#' @param Y multivariate binary data (row for subjects, column for dimensions)
#' @param zsa cluster assignment indicators prior to updating (transition from in MH)
#' @param zsb cluster assignment indicators after updating (transition to in MH; same below.)
#'            both zsa and zsb will be set to zs when active is true, because we want the zs to be
#'            modified/updated.
#' @param cia,cja cluster indicator for subject i,j prior to updating
#' @param cib,cjb cluster indicator for subject i,j after updating
#' @param ni,nj the number of subects clustered with i (j) - including i (j).
#' @param i,j two distinct subejcts' indices
#' @param S the set of indices clustered with either i or j (including i and j)
#' @param ns size of S
#' @param b gamma in the dirichlet distribution in MFM.
#' @param active TRUE for performing a split update (both updating cluster indicators and computing the transition probabilities), FALSE otherwise (only for calculating
#' proposal density; used in q(Z|Z_merge) = q(Z|Z_launch)).
#' @param Q a Q matrix of dimension M by L
#' @param p a vector of machine (factor) prevalences of length M
#' @param theta a vector of true positive rates of length L
#' @param psi a vector of false positive rates of length L
#'
#'
#' @return \itemize{
#' \item \code{log_p} the log probability of transitioning from zsa to zsb (product of all subjects in S);
#' \item \code{ni,nj} the number of subjects assigned to cluster for i,j
#' \item \code{zsb} updated assignment indicators for subjects in S
#' }
#'
#' @export
restricted_gibbs <- function(Y,zsa,zsb,cia,cib,cja,cjb,ni,nj,i,j,S,ns,b,active,Q,p,theta,psi){
  # ??need to update the observations to be in sync with zs.
  # ?? i must be different from j?
  log_p <- 0 # variable for log transition probabilities.
  for (ks in 1:ns){ # beging iteration over observations in S.
    k <- S[ks]
    if (k!=i && k!=j){
      if (zsa[k]==cia){# <--mod # if the current subject to be updated beldongs to cluster of i:
        ni <- ni-1 # <--mod
      } else{
        nj <- nj-1 # <--mod
      }
      # compute probability to assign observation k to the cluster to which obs i belongs:
      Li <- log_marginal(rbind(Y[(zsa==cia)[-k],,drop=FALSE],Y[k,]),Q,p,theta,psi) - log_marginal(Y[(zsa==cia)[-k],,drop=FALSE],Q,p,theta,psi) # from restaurant process.
      Lj <- log_marginal(rbind(Y[(zsa==cja)[-k],,drop=FALSE],Y[k,]),Q,p,theta,psi) - log_marginal(Y[(zsa==cja)[-k],,drop=FALSE],Q,p,theta,psi)
      Pi <- exp(log(ni+b)+Li- matrixStats::logSumExp(c(log(ni+b)+Li,log(nj+b)+Lj))) # the (ni+b) also comes from restaurant process.

      # if we need to update the assignment indicators:
      if (active){
        zsb[k] <- ifelse(stats::runif(1) < Pi, cib, cjb)
      }
      if (zsb[k]==cib){
        ni <- ni+1
        log_p <-log_p+log(Pi)
      } else{
        nj <- nj+1
        log_p <- log_p + log(1-Pi)
      }
    }
  } #end iteration over observations in S.
  return(list(log_p=log_p,ni=ni,nj=nj,zsb=zsb)) # zsb will be updated and passed to the next intermediate Gibbs scan. zsb indicates which observations are used for calculating log marginal likelihood.
}

#' Split-merge type update of clusterings
#'
#' This function uses restricted_gibbs to arrive at a launch state and then propose
#' a split-merge move. Based on the computed acceptance probability, this function accepts
#' or rejects the proposal arriving at an old or new state for the cluster indicators.
#'
#' @param Y multivariate binary data (row for subjects, column for dimensions)
#' @param z original cluster indicators
#' @param zs temporary variable used for split-merge assignments
#' @param S temporary set of indicies for restricted Gibbs
#' @param mylist the cluster ID list
#' @param N a vector of the numbers of subjects in the list of clusters
#' @param t total number of non-empty clusters
#' @param b gamma parameter in the mixture of finite mixture formulation
#' @param log_v coefficients in MFM (check Miller and Harrison 2016 JASA)
#' @param n_split the number of intermediate Gibbs scan to arrive at the launch split state (usually 5)
#' @param Q a Q matrix of dimension M by L
#' @param p a vector of machine (factor) prevalences of length M
#' @param theta a vector of true positive rates of length L
#' @param psi a vector of false positive rates of length L
#' @return Returns the values at the end of the current iteration \itemize{
#' \item \code{t} the number of clusters;
#' \item \code{z} the cluster assignment indicators for all subjects (taking values from \code{mylist});
#' \item \code{N} the number of subjects for each cluster; please refer to mylist to get the matching cluster IDs;
#' \item \code{mylist} the non-empty and empty cluster IDs.
#' }
#'
#' @export
split_merge <- function(Y,z,zs,S,mylist,N,t,b,log_v,n_split,Q,p,theta,psi){
  # for non-conjugate sampler, there needs to be n_merge
  n <- nrow(Y)
  # randomly choose a pair of indices:
  rand_pair <- sample(n,2,replace=FALSE)
  i <- rand_pair[1]
  j <- rand_pair[2]
  ci0 <- z[i]
  cj0 <- z[j] # original states.

  # set S[1],...,S[ns] to the indices of the observations in clusters ci0 and cj0:
  ns <- 0
  for (k in 1:n){
    if (z[k]==ci0 ||  z[k]==cj0){
      ns <- ns+1
      S[ns] <- k
    }
  }

  # find available cluster IDs for split and merge parameters:
  k <- 1
  while (mylist[k]==k){ k=k+ 1};     cm = k
  while (mylist[k]==k+1){ k =k+ 1};  ci = k+1
  while (mylist[k]==k+2) {k =k+ 1};  cj = k+2

  # merge state:
  # for (ks in 1:ns){
  #   #??
  # }

  # randomly choose the launch split state:
  zs[i] <- ci; ni <- 1
  zs[j] <- cj; nj <- 1

  for (ks in 1:ns){# start with a uniformly chosen split:
    k <- S[ks]
    if (k !=i && k!= j){
      if (stats::runif(1)<0.5){
        zs[k] <- ci
        ni<- ni+1 # need adjoin function to collect observations in the same cluster.
      } else{
        zs[k] <- cj
        nj <- nj+1
      }
    }
  }

  for (rep in 1:n_split){#make several restricte Gibb scan moves:
    res_intermediate_GS <- restricted_gibbs(Y,zs,zs,ci,ci,cj,cj,ni,nj,i,j,S,ns,b,TRUE,Q,p,theta,psi) # make sure arguments are updated.
    log_p <- res_intermediate_GS$log_p
    ni <- res_intermediate_GS$ni
    nj <- res_intermediate_GS$nj
    zs <- res_intermediate_GS$zsb
  }

  # make MH proposal:
  if (ci0==cj0){#propose a split:
    # make one final sweep and compute the probability density:
    res_final_GS <- restricted_gibbs(Y,zs,zs,ci,ci,cj,cj,ni,nj,i,j,S,ns,b,TRUE,Q,p,theta,psi) # make sure arguments are updated.
    log_prop_ab <- res_intermediate_GS$log_p
    ni <- res_intermediate_GS$ni
    nj <- res_intermediate_GS$nj
    zs <- res_intermediate_GS$zsb

    # probability of transitioning from the merge state to original state:
    log_prop_ba <- 0 #log(1). Because given current split state, with probability 1 we get the merged state.

    # compute MH acceptance probability
    log_prior_b <- log_v[t+1]+lgamma(ni+b)+lgamma(nj+b)-2*lgamma(b)
    log_prior_a <- log_v[t] + lgamma(ns+b)-lgamma(b)
    log_lik_ratio <- log_marginal(Y[zs==ci,,drop=FALSE],Q,p,theta,psi)+log_marginal(Y[zs==cj,,drop=FALSE],Q,p,theta,psi)-
      log_marginal(Y[z==ci0,,drop=FALSE],Q,p,theta,psi)
    p_accept <- min(1,exp(log_prop_ba-log_prop_ab+log_prior_b-log_prior_a+log_lik_ratio))

    # accept or reject:
    if (stats::runif(1)<p_accept){#accept split:
      for (ks in 1:ns){
        z[S[ks]] <- zs[S[ks]]
      }
      mylist <- ordered_remove(ci0,mylist,t)
      mylist <- ordered_insert(ci,mylist,t-1)
      mylist <- ordered_insert(cj,mylist,t)
      N[ci0] <- 0
      N[ci]  <- ni
      N[cj]  <- nj
      t <- t+1
    }
  } else{ # propose a merge:
    # probability of transitioning to merge state:
    log_prop_ab <- 0  # log(1)

    # compute probability density of going from split launch state to original state:
    res_imaginary_GS <- restricted_gibbs(Y,zs,z,ci,ci0,cj,cj0,ni,nj,i,j,S,ns,b,FALSE,Q,p,theta,psi)
    log_prop_ba <- res_imaginary_GS$log_p
    ni <- res_imaginary_GS$ni
    nj <- res_imaginary_GS$nj

    # compute acceptance probability:
    log_prior_b <- log_v[t-1]+lgamma(ns+b)-lgamma(b)
    log_prior_a <- log_v[t]+lgamma(ni+b)+lgamma(nj+b)-2*lgamma(b)
    log_lik_ratio <- log_marginal(Y[S,,drop=FALSE],Q,p,theta,psi) -
      log_marginal(Y[z==ci0,,drop=FALSE],Q,p,theta,psi)-log_marginal(Y[z==cj0,,drop=FALSE],Q,p,theta,psi) # computed for original (not launch state) and proposed states.
    p_accept <- min(1,log_prop_ba-log_prop_ab+log_prior_b-log_prior_a+log_lik_ratio)

    #accept or reject:
    if (stats::runif(1)<p_accept){
      for (ks in 1:ns){
        z[S[ks]] <- cm
      }
      mylist <- ordered_remove(ci0,mylist,t)
      mylist <- ordered_remove(cj0,mylist,t-1)
      mylist <- ordered_insert(cm,mylist,t-2)
      N[cm] <- ns
      N[ci0] <- 0
      N[cj0] <- 0
      t <- t-1
    }
  }

  return(list(t=t,z=z,N=N,mylist=mylist)) # be extremely clear about what variables got updated.
}













