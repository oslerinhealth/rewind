# initializations:

#' Randomly generate a Q matrix within the identifying conditions
#'
#' This function is used to initialize the MCMC sampling chain within identifiability
#' cosntraints
#'
#' @param M the number of factors (machines)
#' @param L the dimension of the binary data
#' @param p the Bernoulli probability of an one in the Q matrix
#' @examples
#'
#' simulate_Q(3,100)
#'
#' @return a binary matrix of dimension M by L
#' @export
simulate_Q <- function(M,L,p=0.1){
  Q_col_ordered <-cbind(diag(1,M),diag(1,M),matrix(rbinom(M*(L-2*M),1,p),nrow=M))
  Q_col_ordered[which(rowSums(Q_col_ordered[,-(1:(2*M))]) == 0), sample((2*M+1):(L),1)] <- 1
  Q_col_ordered[,sample(1:L)]
}

#' Randomly generate a Q matrix within the identifying conditions driven by data
#'
#' This function is used to initialize the MCMC sampling chain within identifiability
#' cosntraints; but this is smart - because it will not assign a one to a dimension with
#' few ones in the data
#'
#' @param M the number of factors (machines)
#' @param dat a matrix of binary data (rows for observations, columns for dimensions)
#' @param p the Bernoulli probability of an one in the Q matrix
#' @examples
#'
#' simu     <- simulate_data(options_sim0, SETSEED=TRUE)
#' simu_dat <- simu$datmat
#' simulate_Q_dat(5,simu_dat)
#'
#' @return a binary matrix of dimension M by L
#' @export
simulate_Q_dat <- function(M,dat,p=0.1){
  #dat <- simu_dat
  ind_all <- which(colSums(dat)>nrow(dat)/5) # this initial value is very important!. The extra challenge is to known which column is truly all zeros.
  res <- matrix(0,nrow=M,ncol=ncol(dat))
  for (m in 1:M){
    res[m,sample(ind_all,ceiling(p*length(ind_all)))] <- 1
  }
  res
}


#' main function for MCMC sampling designed for boolean matrix factorization
#'
#'
#' This function perform iterations over MCMC sampling chain.
#' NB: 1) add flexibility to specify other parameters as fixed. 2) sample component-specific
#' parameters. 3) sample other model parameters. 4) add timing and printing functionality.
#' 5) add posterior summary functions.
#'
#' @param dat multivariate binary data (row for observations and column for dimensions4)
#' @param model_options Specifying assumed model options:
#' \itemize{
#' \item \code{n} The number of subjects.
#' \item \code{t_max} The maximum (guessed) number of clusters in the data during
#' posterior inference
#' \item \code{m_max} The maximum (guessed) number of factors during
#' the posterior inference
#' \item \code{a_theta, a_psi} hyperparameters for true and false positive rates;
#' a_theta and a_psi are both a vector of length two.
#' \item \code{log_v} The charaster string representing the prior
#' distribution for the number of true clusters, e.g.,
#' \code{"function(k) {log(0.1) + (k-1)*log(0.9)}"}. It is usually pre-computed log of the coefficients in Mixture of Finite Mixtures
#' (Miller and Harrison, 2017, JASA). Use this code:
#' \code{coefficients(eval(parse(text=model_options0$log_pk)),
#' model_options0$gamma,
#' model_options0$n,
#' model_options0$t_max+1)}
#' }
#' The following are used if one needs to prefix a few unnown parameters to
#' their respective true or other values
#' \itemize{
#' \item \code{Q} Q matrix
#' \item \code{theta} a vector of true positive rates
#' \item \code{psi} a vector of false positive rates
#' \item \code{p} a vector of machine prevalences
#' }
#' options for specifying data, sample size, max cluster number,
#' coefficients in MFM (Miller and Harrison 2017 JASA), gamma parameter in the MFM
#' Dirchlet prior, number of intermediate Gibbs scan to arrive at the launch state,
#' and other hyperparamter specification if needed, n_total for total number of
#' mcmc iterations and n_keep for the number of samples kept for posterior inference.
#' Note that the options involve other parameters for sampling hyperparameters such as
#' alpha in the Indian Buffet Process.
#' @param mcmc_options Options for MCMC sampling:
#' \itemize{
#' \item \code{n_total} total number of MCMC iterations
#' \item \code{n_keep} number of iterations kept
#' \item \code{n_split} the number of restricted Gibbs scan to arrive at a launch state;
#' see \link{restricted_gibbs}
#' \item \code{block_update_Q} update columns of Q or not. Then no constraints
#' about identifiability is imposed upon Q at any iterations.
#' }
#'
#' @return posterior samples for quantities of interest. It is a list comprised of the following elements
#' \itemize{
#' \item \code{t_samp}
#' \item \code{z_samp}
#' \item \code{N_samp}
#' \item \code{keepers} indices of MCMC samples kept for inference
#' \item \code{H_star_samp}
#' \item \code{alpha_samp}
#' } The following are recorded if they are not fixed in a priori:
#' \itemize{
#' \item \code{Q_samp}
#' \item \code{theta_samp}
#' \item \code{psi_samp}
#' \item \code{p_samp}
#' }
#' @export
sampler <- function(dat,model_options,mcmc_options){
  # dat = simu_dat
  # model_options = model_options0
  # mcmc_options = mcmc_options0
  #
  n <- nrow(dat) # number of observations
  L <- ncol(dat) # number of dimensions (protein landmarks).

  n_total <- mcmc_options$n_total # total number of mcmc iterations.
  n_keep  <- mcmc_options$n_keep  # toral number of samples kept for posterior inference.
  keepers <- seq(ceiling(n_total/n_keep),n_total,len=n_keep)
  keep_index <- 0                 # the index to keep during MCMC inference.
  n_split <- mcmc_options$n_split # number of intermediate GS scan to arrive at launch state.
  block_update_Q <- mcmc_options$block_update_Q # TRUE for updating columns of Q. FALSE otherwise.

  # set options (need to fill in as specific paramters are needed):
  t_max <- model_options$t_max # maximum allowable number of clusters.
  m_max <- model_options$m_max # maximum number of machines (truncated to m_max).
  b     <- model_options$b # gamma parameter in MFM - hyperparameter for Dirichlet distn. needs to be sampled?????
  log_v <- model_options$log_v #coefficients for MFM.

  if (is.null(model_options$Q)){
    Q_samp <- array(0,c(m_max,ncol(dat),n_keep))
    #Q      <- simulate_Q(m_max,L) # random initialization, could be impproved by a smarter starting Q matrix.
    Q      <- simulate_Q_dat(m_max,dat,0.5) # random initialization, could be impproved by a smarter starting Q matrix.
  }else{
    Q     <- model_options$Q # suppose we are given Q.
  }

  H_star_enumerate <- as.matrix(expand.grid(rep(list(0:1), m_max)),ncol=m_max) # all binary patterns for machine usage profiles. 2^m_max of them.

  # initialize the sampling chain:
  t <- 1        # number of clusters.
  z <- rep(1,n) # z[i] is the cluster ID for observation i.
  mylist <- rep(0,t_max+3); mylist[1] <- 1  # mylist[1:t] is the list of active cluster IDs.
  # mylist is maintained in increasing order for 1:t, and
  # is 0 after that.
  c_next <- 2 # an available cluster ID to be used next.
  N <- rep(0,t_max+3); N[1] <- n # N[c] is the size of cluster c. Note that N is different from n.

  # H_star <- matrix(0,nrow=t_max+3,ncol=m_max) # parameter for clusters. In our application, the machine usage profiles.
  # NOTE: need to consider SMARTER initial values for this. How to sample H_star? Not needed in conjugate algorithm.

  log_p <- rep(0,n+1) # probability for sampling z; n+1 because at most there could be n clusters;
  # and sometimes one needs to assign an observation to a new cluster - by current bookkeeping,
  # we are not efficient in memory and create a new cluster ID - hence +1.
  zs <- rep(1,n)  # temporary cluster indicators for split-merge assignments.
  S  <- rep(0,n)  # temporary variable for indicies used during split-merge step.

  log_Nb <- log(1:n)+b # the multipliers needed when assigning an observation to an existing cluster. Restaurant process stuff.

  # variables for samples kept:
  t_samp <- rep(0,n_total)
  N_samp <- matrix(0,nrow=t_max+3,ncol=n_total)
  z_samp <- matrix(0,nrow=n,ncol=n_keep) # posterior samples of cluster indicators.
  H_star_samp <- array(0,c(t_max+3,m_max,n_total))# component-specific machine profiles.

  if (is.null(model_options$theta)){
    theta_samp <- matrix(0,nrow=ncol(dat),ncol=n_keep)
    a_theta <- model_options$a_theta
    theta <- sapply(1:L,function(i){rbeta(1,a_theta[1],a_theta[2])}) # initialization.
  } else{
    theta   <- model_options$theta # suppose we are given Q; use specified theta.
  }
  if (is.null(model_options$psi)){
    psi_samp   <- matrix(0,nrow=ncol(dat),ncol=n_keep)
    a_psi      <- model_options$a_psi
    psi <- sapply(1:L,function(i){rbeta(1,a_psi[1],a_psi[2])}) # initialization.
  } else{
    psi     <- model_options$psi # suppose we are given Q; use specified psi.
  }

  if (is.null(model_options$p0)){
    p_samp <- matrix(0,nrow=m_max,ncol=n_total)
    p      <- rep(0.5,m_max) # initialization.
  } else{
    p      <- model_options$p0 # fix prevalence.
  }

  if (is.null(model_options$alpha)){
    alpha_samp <- rep(0,n_total)  # hyperparameter for machine prevalence.
    alpha      <- m_max
  }

  cat("==[rewind] Start MCMC: ==\n")
  for (iter in 1:n_total){
    if (iter%%10==0){cat("==[rewind] Iteration: ", iter, "==\n");
      cat("==[rewind]factor indicators (machines usages) for t=",t," clusters:==\n")
      print(H_star[,colSums(H_star)!=0,drop=FALSE])
      image(Q[colSums(H_star)!=0,drop=FALSE,])
    }

    # update cluster indicators z for all subjects - one complete Gibbs scan to refine clusters:
    for (i in 1:n){ # iterate over subjects:
      # remove obs i from its current cluster:
      c    <- z[i]
      N[c] <- N[c]-1
      if (N[c]>0){ # if there are observations left after removal of i, use c_next:
        c_prop <- c_next
      } else{# if i was by itself, use its cluster ID again.
        c_prop <- c
        mylist <- ordered_remove(c,mylist,t)
        t <- t-1 # total number of clusters is now reduced by one.
      }
      # compute probability for Gibbs updating - the probability of assigning a subject
      # to a cluster.
      for (j in 1:t){
        cc <- mylist[j]
        log_p[j] <- log_Nb[N[cc]]+log_marginal(rbind(dat[(z==cc)[-i],,drop=FALSE],dat[i,]),Q,p,theta,psi)-
          log_marginal(dat[(z==cc)[-i],,drop=FALSE],Q,p, theta,psi) # existing cluster.
      }
      log_p[t+1] <- log_v[t+1]-log_v[t] + log(b) + log_marginal(dat[i,,drop=FALSE],Q,p,theta,psi) # new cluster.

      j <- sample(t+1,1,prob = exp(log_p[1:(t+1)]))

      # add obs i to its new cluster:
      if (j<=t){
        c <- mylist[j]
      } else{
        c      <- c_prop
        mylist <- ordered_insert(c,mylist,t)
        t      <- t+1
        c_next <- ordered_next(mylist)
        if(t>t_max) {stop("==[rewind]Sampled t has exceeded t_max. Increast t_max and retry.==")}
      }
      z[i] <- c
      N[c] <- N[c] + 1

    }# END iteration over subjects.

    # if (use_splitmerge){}
    res_split_merge <- split_merge(dat,z,zs,S,mylist,N,t,b,log_v,n_split,Q,p,theta,psi)
    t <- res_split_merge$t
    z <- res_split_merge$z
    N <- res_split_merge$N
    mylist <- res_split_merge$mylist
    c_next <- ordered_next(mylist)
    if(t>t_max) {stop("==[rewind]Sampled t has exceeded t_max. Increast t_max and retry.==")}
    cat("==[rewind]current number of subjects per cluster: t=",t," clusters:", N[N!=0],"==\n")

    # # update xi_star matrix of dimension K times L:
    # #if (is.null(model_options$Q)){
    # xi_samp <- array(0,c(t_max+3,L,n_total)) # add this to the results.
    # for (j in 1:t){
    #   xi_samp[mylist[j],,iter]  <- update_xi(dat[z==mylist[j],,drop=FALSE],p,theta,psi)
    # }
    # xi_star <- xi_samp[mylist[1:t],,iter]
    # #}

    # update machine usage profile given clusters: (check here for sampling H)
    for (j in 1:t){
      # for (m in 1:m_max){
      #   H_star_samp[mylist[j],m,iter] <- 0
      #   L0            <- log_full(dat[z==mylist[j],,drop=FALSE],matrix(H_star_samp[mylist[j],,iter],nrow=1),Q,p,theta,psi)
      #   H_star_samp[mylist[j],m,iter] <- 1
      #   L1            <- log_full(dat[z==mylist[j],,drop=FALSE],matrix(H_star_samp[mylist[j],,iter],nrow=1),Q,p,theta,psi)
      #   curr_eta_p    <- exp(L1-matrixStats::logSumExp(c(L0,L1)))
      #   H_star_samp[mylist[j],m,iter] <-rbinom(1,1,curr_eta_p)
      #
      #   #H_star_samp[mylist[j],m,iter] <- metrop_flip(H_star_samp[mylist[j],m,iter],curr_eta_p)
      # }
      curr_pattern_log_p   <- log_full(dat[z==mylist[j],,drop=FALSE],H_star_enumerate,Q,p,theta,psi)
      curr_pattern_p       <- exp(curr_pattern_log_p-matrixStats::logSumExp(curr_pattern_log_p))
      curr_ind             <- sample(1:(2^m_max),1,prob=curr_pattern_p)
      H_star_samp[mylist[j],,iter] <- H_star_enumerate[curr_ind,]
    }

    H_star <- matrix(H_star_samp[mylist[1:t],,iter],ncol=m_max)

    # update Q:
    if (is.null(model_options$Q)){
      #if(iter%%5==0){
      if (block_update_Q){
      Q <- update_Q_col_block(dat,Q,H_star,z,t,mylist,p,theta,psi)
      } else{
      #}else{
      Q <- update_Q(dat,Q,H_star,z,t,mylist,p,theta,psi)
      }
      #}
      #Q <- update_Q_no_H(dat,Q,z,t,mylist,p,theta,psi)
      # if (sum(colSums(H_star)==0)>0){
      #   Q[colSums(H_star)==0,] <- simulate_Q_dat(sum(colSums(H_star)==0),dat)# # reinitialize if there are unused machines vs sample from prior?. Now we can see
      # }
    }

    # update true/false positive rates - theta/psi:
    if (is.null(model_options$theta) && is.null(model_options$psi)){
      res_update_pr <- update_positive_rate(dat,H_star_samp[z,,iter],Q,a_theta,a_psi)
      theta <- res_update_pr$theta
      psi   <- res_update_pr$psi
    }

    # update hyperparameter for {p_m}:
    if (is.null(model_options$alpha)){
      alpha <- update_alpha(H_star,t,m_max)
      alpha_samp[iter] <- alpha
    }

    # update prevalence vector {p_m}:
    if (is.null(model_options$p0)){
      p    <- update_prevalence(H_star,alpha,m_max)
      p_samp[,iter]    <- p
    }

    #
    # Record MCMC results:
    #
    t_samp[iter] <- t
    for (j in 1:t){
      N_samp[mylist[j],iter] <- N[mylist[j]]
    }

    if (iter==keepers[keep_index+1]){
      keep_index <- keep_index +1
      for (i in 1:n){
        z_samp[i,keep_index] <- z[i]
      }
      if (is.null(model_options$theta)){
        theta_samp[,keep_index] <- theta # of length L.
      }
      if (is.null(model_options$psi)){
        psi_samp[,keep_index]   <- psi
      }
      if (is.null(model_options$Q)){
        Q_samp[,,keep_index] <- Q
      }
    }
  }# END mcmc iterations.

  res <- list(t_samp=t_samp,N_samp=N_samp,z_samp=z_samp,keepers=keepers,
              H_star_samp = H_star_samp,
              alpha_samp=alpha_samp)
  if (is.null(model_options$theta)){
    res$theta_samp <- theta_samp # of length L.
  }
  if (is.null(model_options$psi)){
    res$psi_samp   <- psi_samp
  }
  if (is.null(model_options$Q)){
    res$Q_samp <- Q_samp
    #res$xi_samp <- xi_samp
  }
  if (is.null(model_options$p0)){
    res$p_samp <- p_samp
  }
  res
}
