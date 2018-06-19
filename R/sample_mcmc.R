# initializations:

#' Randomly generate a Q matrix within the identifying conditions
#'
#' This function initializes Q (if unknown) for MCMC sampling within identifiability
#' constraints
#'
#' @param M latent state dimension
#' @param L dimension of the binary responses
#' @param p Bernoulli probability of 1 in the Q matrix (except two diagonal matrices)
#' @examples
#'
#' simulate_Q(3,100)
#'
#' @return a binary matrix of dimension M by L
#' @export
simulate_Q <- function(M,L,p=0.1){
  Q_col_ordered <-cbind(diag(1,M),diag(1,M),matrix(stats::rbinom(M*(L-2*M),1,p),nrow=M))
  Q_col_ordered[which(rowSums(Q_col_ordered[,-(1:(2*M))]) == 0), sample((2*M+1):(L),1)] <- 1
  # if a row sums to two, add an extra "1".
  Q_col_ordered[,sample(1:L)]
}

#'Randomly generate a Q matrix within the identifying condition (data-driven)
#'
#'This function initializes Q (if unknown) during MCMC sampling chain within identifiability
#'cosntraints. It is a warm start - because it will not assign a one to a dimension with
#'few ones in the data
#'
#' @param M latent state dimension
#' @param dat binary data matrix (rows for observations, columns for dimensions)
#' @param p Bernoulli probability of 1 in the Q matrix (except two diagonal matrices)
#' @param frac A threshold - this function only initializes the dimensions with at least
#' \code{frac}*100\% observed frequencies. Default: 20\%.
#'
#' @examples
#' # simulate data:
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
#' simu     <- simulate_data(options_sim0, SETSEED=TRUE)
#' simu_dat <- simu$datmat
#' simulate_Q_dat(5,simu_dat)
#'
#' @return a binary matrix of dimension M by L
#' @export
simulate_Q_dat <- function(M,dat,p=0.1,frac=1/4){
  ind_all <- which(colSums(dat)>nrow(dat)*frac) # this initial value helps MCMC sampling!.
  res <- matrix(0,nrow=M,ncol=ncol(dat))
  for (m in 1:M){
    res[m,sample(ind_all,ceiling(p*length(ind_all)))] <- 1
  }
  res
}

# compute the likelihood (in the cpp file):


# update parameters:

# eta_star <- as.matrix(expand.grid(rep(list(0:1), options_sim0$M)),ncol=options_sim0$M)
# sum(abs(log_full0(Y,eta_star,Q, p, theta, psi)-
# log_full(Y,eta_star,Q, p, theta, psi))) # perhaps write this into a unit test.

#' update true and false positive rates WITHOUT constraints
#'
#' This function samples the true positive rates (theta) and false positive rates (psi)
#'
#' @param Y binary data matrix
#' @param H a matrix whose rows are individual specific latent states (all subjects)
#' @param Q Q matrix
#' @param a_theta,a_psi hyperparameters for priors on theta and psi, respectively.
#'
#' @return a list of true positive rates theta and false positive rates psi, each of length L.
#' @export
update_positive_rate <- function(Y,H,Q,a_theta,a_psi){
  xi       <- (H%*%Q>0.5)+0 # <-- NB: for DINO model. Essentially the design matrix
  # in the manuscript. Should be modified for general restricted LCM.
  psi_a1   <- colSums((1-xi)*Y)+a_psi[1,]
  psi_a2   <- colSums((1-xi)*(1-Y))+a_psi[2,]
  theta_a1 <- colSums(xi*Y)+a_theta[1,]
  theta_a2 <- colSums(xi*(1-Y))+a_theta[2,]
  psi_samp   <- sapply(1:ncol(Q),function(i) {stats::rbeta(1,psi_a1[i],psi_a2[i])})

  # theta_all_a1 <- sum(colSums(xi*Y))+a_theta[1,1]
  # theta_all_a2 <- sum(colSums(xi*(1-Y)))+a_theta[2,1]
  # theta_samp   <- rep(stats::rbeta(1,theta_all_a1,theta_all_a2), ncol(Q))

  #theta_samp <- sapply(1:ncol(Q),function(i) {stats::rbeta(1,theta_a1[i],theta_a2[i])})
  theta_samp <- sapply(1:ncol(Q),function(i) {#stats::rbeta(1,theta_a1[i],theta_a2[i])
    u <- runif(1,pbeta(psi_samp[i],theta_a1[i],theta_a2[i]),1)
    #u <- runif(1,pbeta(0.5,theta_a1[i],theta_a2[i]),1) # set the lower limit to be 0.5.
    max(min(0.999999,qbeta(u,theta_a1[i],theta_a2[i])),0.000001)
    #cat(u," - u\n")
    #cat(res,"\n")
    #res
  })
  list(theta=theta_samp,psi=psi_samp)
}

#' sample alpha - hyperparameter for latent state prevalences p_m
#'
#' Function to sample the parameters from a grid (this is used only for model
#' with pre-specified latent state dimension M)
#'
#'@param H_star the matrix of latent state profiles across clusters
#'@param t number of clusters
#'@param M latent state dimension
#'@param a,b hyperparameter for Beta distribution over reparameterized alpha.
#'@param show_density Default to FALSE - hide the full conditional density of alpha
#'given other unknown parameters; TRUE otherwise.
#'
#'@return an updated alpha value (positive)
#'@export
update_alpha <- function(H_star,t,M,a=1,b=1,show_density=FALSE) {
  th    <- seq(0.001,.999,length=5000) # grid over (0,1) after reparametrization.
  sm    <- apply(H_star,2,sum,na.rm=T)
  part1 <- 0
  for (m in 1:M){
    part1 <- part1+lgamma(th/(1-th)/M+sm[m])
  }
  tmp    <- stats::dbeta(th,a,b,log=T)+M*log(th/(1-th))+part1-M*lgamma(th/(1-th)/M+t+1)
  tmp    <- exp(tmp-max(tmp))
  rtheta <- sample(th,1,prob=tmp)
  if (show_density){
    graphics::plot(th,tmp,type="l")
    graphics::abline(v=rtheta,col="red",lty=2)
  }
  (rtheta)/(1-rtheta)
}

#log_full(simu_dat,as.matrix(expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))),
#         Q,p,theta,psi)

#' update latent state prevalences
#'
#' This function updates the prevelance paramters
#'
#' @param H_star latent state profiles for all clusters
#' @param alpha hyperparameter
#' @param M dimension of the latent state vector (number of columns for H_star)
#'
#' @return vector of length M (each between 0 and 1).
#' @export
update_prevalence <- function(H_star,alpha,M){
  n1_star <- apply(H_star,2,sum,na.rm=T)
  n0_star <- apply(1-H_star,2,sum,na.rm=T)
  sapply(1:ncol(H_star),function(i) {stats::rbeta(1,n1_star[i]+alpha/M, n0_star[i]+1)})
  # based on finite dimension IBP construction.
}

#' Update the current element of Q_ml or not
#'
#' Function to test whether an update of the current element of Q_ml is needed.
#' This function must be used for the constrained Gibbs sampler.
#'
#' @param Q a matrix with rows corresponding to latent states and columns
#'          for number of binary measurements
#' @param k,l row and column indices for checking wether to update Q_kl
#' @return A logical value. TRUE for updating in constrained Gibbs sampler, FALSE
#' for skipping the Gibbs update.
#' @references Chen, Y., Culpepper, S. A., Chen, Y., and Douglas, J. (2017). Bayesian estimation of the DINA Q matrix.
#' @export
do_update_Q_one <- function(Q,k,l){ # this saves lots of speed:
  test <- rep(NA,3)
  ind_unit_Q <- colSums(Q)==1
  test[1] <- equal_unit_vec(Q[,l],k)
  test[2] <- (sum(Q[k,])==3) && (Q[k,l]==1)
  Q_cand  <- Q[,ind_unit_Q,drop=FALSE]
  test[3] <- (Q[k,l]==0) && sum(Q[,l])==1 &&
    (sum(sapply(1:ncol(Q_cand), function(c) equal_unit_vec(Q_cand[,c],which(Q[,l]==1))))==2)
  !any(test)
}


#' Compute the full conditional probability of Q_ml given the rest of parameters
#'
#' Function to compute this log full conditional density, which will be used
#' in \link{update_Q}
#'
#' @param Yl a column of the multivariate binary data
#' @param eta_star a matrix of latent state profiles, rows for clusters,
#' columns for M latent states
#' @param Ql l-th column of the Q matrix
#' @param thetal,psil The true, false positive rate (both between 0 and 1)
#'
#' @return log conditional probability of Q_ml given other unknown parameters
#' @export
log_pr_Qml_cond <- function(Yl,eta_star,Ql,thetal,psil){
  n1 <- sum(Yl)
  n0 <- sum(1-Yl)
  xil <- (eta_star%*%matrix(Ql,ncol=1) > 0.5)+0 # NB: works for DINA model;
  # need revision for general restricted LCMs.
  PR_mat <- (1-xil)*(n1*log(psil)+n0*log(1-psil))+
    xil*(n1*log(thetal)+n0*log(1-thetal)) # conditional distribution given the
  # design matrix.
}

#' update the Q matrix element-wise
#'
#' This function updates Q matrix by Gibbs sampler (with option to do constarined
#' updates within the identifiability constraint)
#' NB: - do we need M here? do we modify M after collapsing partner machines.
#'     - need to speed up.
#' @param Y binary data matrix
#' @param Q_old the Q matrix from the last scan in Gibbs sampler (of dimension M by L)
#' @param H matrix of latent state profiles, rows for N subjects,
#'          columns for M latent states
#' @param z a vector of individual cluster indicators
#' @param t the number of clusters at a particular iteration
#' @param mylist the ordered list that keeps the cluster ids in increasing order
#' @param p Latent state prevalences; a vector of length identical to the columns of H.
#' @param theta,psi True and false positive rates. Both are vectors of length L
#' @param constrained Default to FALSE; Set to TRUE if doing constrained Gibbs sampling
#' with identifiability constraints.
#'
#' @return A Q matrix with all its elements updated.
#' @export
update_Q <- function(Y,Q_old,H,z,t,mylist,p,theta,psi,constrained=FALSE){
  M  <- nrow(Q_old)
  N  <- nrow(Y)
  L  <- ncol(Y)
  if (constrained){ # constrained Gibbs sampling:
    for (k in sample(1:M,replace=FALSE)){
      for (l in sample(1:L,replace=FALSE)){ # begin iteration over elements.
        if (do_update_Q_one(Q_old,k,l)){ # begin an update if needed.
          L0 <- L1   <- 0
          Q_old[k,l] <- 0
          #interact   <- - 0.1 # <--- negative value to encourage sparse patterns per dimension.
          #vec_interact <- c(Q_old[k,])
          #vec_interact <- c(Q_old[k,],Q_old[-k,l])
          #ising_tmp  <- sum(outer(vec_interact,vec_interact,"*")[upper.tri(outer(vec_interact,vec_interact,"*"))])
          for (j in 1:t){ #Yl,eta_star,Ql,thetal,psil
            L0 <- L0 + log_pr_Qml_cond(Y[z==mylist[j],l,drop=FALSE],
                                       H[j,,drop=FALSE],Q_old[,l],theta[l],psi[l])
          }
          L0 <- L0 # + dbinom(sum(vec_interact),L-2*M,0.1,log=TRUE) # <-- encourage sparse columns.
          #L0 <- L0  - log(1+exp(-interact*ising_tmp)) # <-- encourage sparse columns.
          Q_old[k,l] <- 1
          #vec_interact <- c(Q_old[k,])
          #vec_interact <- c(Q_old[k,],Q_old[-k,l])
          #ising_tmp <- sum(outer(vec_interact,vec_interact,"*")[upper.tri(outer(vec_interact,vec_interact,"*"))])
          for (j in 1:t){
            L1 <- L1 + log_pr_Qml_cond(Y[z==mylist[j],l,drop=FALSE],
                                       H[j,,drop=FALSE],Q_old[,l],theta[l],psi[l])
          }
          L1 <- L1  #+ dbinom(sum(vec_interact),L-2*M,0.1,log=TRUE)# <-- encourage sparse columns.
          #L1 <- L1 - log(1+exp(-interact*ising_tmp)) # <-- encourage sparse columns.
          curr_prob <- exp(L1- matrixStats::logSumExp(c(L0,L1))) 
          #print(curr_prob)
          #Q_old[k,l] <- metrop_flip(Q_old[k,l],curr_prob) # <-- if doing metroplized flipping.
          Q_old[k,l] <- stats::rbinom(1,1,prob = curr_prob)
        }# end an update if needed.
      }
    }  #end iteration over elements.
  } else { # Gibbs update without constraints of Q within a set.
    for (k in sample(1:M,replace=FALSE)){
      for (l in sample(1:L,replace=FALSE)){ # begin iteration over elements.
        L0 <- L1 <- 0
        Q_old[k,l] <- 0
        for (j in 1:t){ #Yl,eta_star,Ql,thetal,psil
          L0 <- L0 + log_pr_Qml_cond(Y[z==mylist[j],l,drop=FALSE],
                                     H[j,,drop=FALSE],Q_old[,l],theta[l],psi[l])
        }
        Q_old[k,l] <- 1
        for (j in 1:t){
          L1 <- L1 + log_pr_Qml_cond(Y[z==mylist[j],l,drop=FALSE],
                                     H[j,,drop=FALSE],Q_old[,l],theta[l],psi[l])
        }
        curr_prob <- exp(L1- matrixStats::logSumExp(c(L0,L1)))
        #print(curr_prob)

        #Q_old[k,l] <- metrop_flip(Q_old[k,l],curr_prob)
        Q_old[k,l] <- stats::rbinom(1,1,prob = curr_prob)
      }
    }  #end iteration over elements.
  }
  Q_old
}

#' Compute the full conditional probability of column Q_l given other unknown parameters
#'
#' Function to compute the log full conditional density for all patterns of the
#' l-th column of Q. This function is used in \link{update_Q_col_block}.
#'
#' @param Yl a column of the binary data matrix
#' @param eta_star a matrix of machine usage indicators, rows for clusters, columns for M machines
#' @param Ql_enumerate l-th column of the Q matrix
#' @param thetal,psil True and false positive rates. Both are vectors of length L
#'
#' @return a vector of log conditional probability of column Q_l taking each of
#' \code{2^nrow(Ql_enumerate)} given other unknown parameters,
#' the dimension will be identical to \code{2^nrow(Ql_enumerate)}.
#' @export
log_pr_Qml_cond_enum <- function(Yl,eta_star,Ql_enumerate,thetal,psil){
  n1 <- sum(Yl)
  n0 <- sum(1-Yl)
  xil <- (eta_star%*%matrix(Ql_enumerate,nrow=ncol(eta_star)) > 0.5)+0 # NB: only works for DINO model; need to revise
  # for other general restricted LCMs.
  PR_mat <- (1-xil)*(n1*log(psil)+n0*log(1-psil))+
    xil*(n1*log(thetal)+n0*log(1-thetal))
}

#' Block update the Q matrix by column
#'
#' This function updates Q matrix by Gibbs sampler (without identifiability
#' constraints)
#'
#' @param Y binary data matrix
#' @param Q_old the Q matrix from the last scan in Gibbs sampler (of dimension M by L)
#' @param H matrix of latent state profiles, rows for N subjects, columns for M latent states
#' @param z a vector of individual cluster indicators
#' @param t the number of clusters at an iteration
#' @param mylist the ordered list that keeps the cluster ids in increasing order
#' @param p latent state prevalences; a vector of length identical to the columns of H.
#' @param theta,psi True and false positive rates, each of length L
#'
#' @return A Q matrix with all elements updated.
#' @export
update_Q_col_block <- function(Y,Q_old,H,z,t,mylist,p,theta,psi){
  M  <- nrow(Q_old)
  N  <- nrow(Y)
  L  <- ncol(Y)
  Ql_enum  <- t(as.matrix(expand.grid(rep(list(0:1), M)),ncol=M))
  for (l in sample(1:L,replace=FALSE)){ # begin iteration over columns:
    L_enum  <- rep(0,2^M)
    for (j in 1:t){
      L_enum <- L_enum + log_pr_Qml_cond_enum(Y[z==mylist[j],l,drop=FALSE],
                                              H[j,,drop=FALSE],Ql_enum,theta[l],psi[l])
      # for all patterns of each column, we compute the likelihood and sum over clusters.
    }
    curr_prob <- exp(L_enum-matrixStats::logSumExp(L_enum))
    Q_old[,l] <- Ql_enum[,sample(1:2^M,1,prob = curr_prob)]
  }
  Q_old
}

#' MCMC sampling for pre-specified latent state dimension M
#'
#' This function performs MCMC sampling with user-specified options.
#' NB: 1) add flexibility to specify other parameters as fixed. 2) sample component-specific
#' parameters. 3) sample other model parameters. 4) add timing and printing functionality.
#' 5) add posterior summary functions.
#' 6) edit verbose contents.
#'
#' @param dat binary data matrix (row for observations and column for dimensions)
#' @param model_options Specifying assumed model options:
#' \itemize{
#' \item \code{n} The number of subjects.
#' \item \code{t_max} The maximum (guessed) number of clusters in the data during
#' the posterior inference
#' \item \code{m_max} For a model with pre-specified number of factors, \code{m_max};
#' In an infinite dimension model, \code{m_max} is
#' the maximum (guessed) latent state dimension during the posterior inference
#' (see slice_sampler to come); one can increase this number if
#' this pacakge recommends so in the printed message;
#' \item \code{a_theta, a_psi} hyperparameters for true and false positive rates;
#' a_theta and a_psi are both a vector of length two.
#' \item \code{a_alpha, b_alpha} Just for infinite latent state dimension model  -
#' Gamma hyperparameter for the hyperprior on \code{alpha}.
#' (see slice_sampler to come)
#' \item \code{log_v} The charaster string representing the prior
#' distribution for the number of true clusters, e.g.,
#' \code{"function(k) {log(0.1) + (k-1)*log(0.9)}"}. We pre-computed
#' log of the coefficients in Mixture of Finite Mixtures
#' (Miller and Harrison, 2017, JASA). Use this code:
#' \code{coefficients(eval(parse(text=model_options0$log_pk)),
#' model_options0$gamma,
#' model_options0$n,
#' model_options0$t_max+1)}
#' \item \code{the_partition} Put a vector of z here if not to update
#' the partition/clustering. The z must be 1) from 1 to t, which is
#' the number of groups, 2) consecutive from 1 to t.
#' }
#' The following are used if one needs to pre-specify a few unknown parameters to
#' their respective true or other values
#' \itemize{
#' \item \code{Q} Q matrix
#' \item \code{theta} a vector of true positive rates
#' \item \code{psi} a vector of false positive rates
#' \item \code{p} a vector of latent state prevalences
#' \item \code{alpha} For pre-specified latent state dimension, the hyperparameter
#' for \code{Beta(alpha/m_max,1)} (can set to \code{m_max});
#' For infinite dimension model, the hyperparameter for IBP (can set to 1).
#' }
#' Options for specifying data, sample size, max cluster number,
#' coefficients in MFM (Miller and Harrison 2017 JASA), Gamma parameter in the MFM
#' Dirchlet prior, number of intermediate Gibbs scan to arrive at the launch state,
#' and other hyperparamter specification if needed, \code{n_total} for total number of
#' MCMC iterations and \code{n_keep} for the number of samples kept for posterior inference.
#' Note that the options involve other parameters for sampling hyperparameters such as
#' alpha in the Indian Buffet Process.
#' @param mcmc_options Options for MCMC sampling:
#' \itemize{
#' \item \code{n_total} total number of MCMC iterations
#' \item \code{n_keep} number of iterations kept
#' \item \code{n_burnin} the number of burnin iterations to be discarded
#' for drawing the posterior inference.
#' \item \code{n_split} the number of restricted Gibbs scan to arrive at a launch state;
#' see \link{restricted_gibbs}
#' \item \code{print_mod} print intermediate model fitting information
#' \item \code{constrained} update the Q matrix with identifiability constraints (if \code{TRUE});
#' otherwise, set to \code{FALSE}.
#' \item \code{block_update_H} update rows of H (if \code{TRUE}) or not
#' (if \code{NULL} or \code{FALSE} - must be so for slice_sampler to come).
#' \item \code{block_update_Q} update columns of Q (if \code{TRUE}) or not
#' (if \code{NULL} or \code{FALSE} - must be so for slice_sampler to come).
#' Then no identifiability constraint is imposed upon Q at any iterations.
#' \item \code{ALL_IN_ONE_START} \code{TRUE} for putting all subjects in one cluster,
#' \code{FALSE} by starting from a hierechical clustering (complete linkage) and cut
#' to produce floor(t_max/4) clusters. Consider this as a warm start.
#' \item \code{MORE_SPLIT} when getting to launch state of the partition,
#' \code{TRUE} for biasing towards split; FALSE for uniformly choose a pair (i,j)
#' and then deciding to merge (if they belong to distinct clusters)
#' or split (if they belong to the identical cluster)
#' \item \code{hmcols} color palette.
#' }
#'
#' @example /inst/example/simulation_fixed_M.R
#'
#' @return posterior samples for quantities of interest.
#' It is a list comprised of the following elements:
#' \itemize{
#' \item \code{t_samp}
#' \item \code{z_samp}
#' \item \code{N_samp}
#' \item \code{keepers} indices of MCMC samples kept for inference; will be used
#' to see which ones to use for inference after \code{n_burnin}.
#' \item \code{H_star_samp} \code{t_max+3} by \code{m_max} binary matrix
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

  # # <<<<------ testing data:
  # dat <- simu_dat
  # model_options <- model_options0
  # mcmc_options  <- mcmc_options0
  # # <<<<------ end of testing data.

  # # <<<<------ testing data:
  # dat <- simu_dat
  # model_options <- model_options_refit_given_partition
  # mcmc_options  <- mcmc_options_refit_given_partition
  # # <<<<------ end of testing data.

  n <- nrow(dat) # number of observations
  L <- ncol(dat) # number of dimensions (e.g.,protein landmarks).

  ALL_IN_ONE_START <- mcmc_options$ALL_IN_ONE_START # start the chain with all individuals in one cluster.
  n_total <- mcmc_options$n_total # total number of mcmc iterations.
  n_keep  <- mcmc_options$n_keep  # toral number of samples kept for posterior inference.
  keepers <- seq(ceiling(n_total/n_keep),n_total,len=n_keep)
  keep_index <- 0                 # the index to keep during MCMC inference.
  n_split    <- mcmc_options$n_split # number of intermediate GS scan to arrive at launch state.

  do_update_partition <- is.null(model_options$the_partition)
  block_update_Q <- !is.null(mcmc_options$block_update_Q) && mcmc_options$block_update_Q # TRUE for updating columns of Q. FALSE otherwise.
  block_update_H <- !is.null(mcmc_options$block_update_H) && mcmc_options$block_update_H # TRUE for updating rows of H. FALSE otherwise.
  constrained    <- !is.null(mcmc_options$constrained) && mcmc_options$constrained
  hmcols <- mcmc_options$hmcols


  # set options (need to fill in as specific paramters are needed):
  t_max <- model_options$t_max # maximum allowable number of clusters.
  m_max <- model_options$m_max # maximum number of machines (truncated to m_max).
  b     <- model_options$b     # Gamma parameter in MFM - hyperparameter for Dirichlet distn. needs to be sampled?????
  log_v <- model_options$log_v # coefficients for MFM.

  is_identity_Q <- NULL
  if (is.null(model_options$Q)){
    Q_samp <- array(0,c(m_max,L,n_keep))
    Q_merge_samp <- Q_samp
    Q      <- simulate_Q_dat(m_max,dat,0.1,min(max(colMeans(dat)),0.3)) # random initialization; warm start.
  }else{
    Q      <- model_options$Q # if Q is given.
    is_identity_Q <- (nrow(Q)==ncol(Q)) && sum(abs(Q-diag(nrow(Q))))<1e-3
  }

  if (is.null(model_options$Q) | block_update_H){ # <-- if Q may be non-diagonal; or if we block update H.
    H_star_enumerate <- as.matrix(expand.grid(rep(list(0:1), m_max)),ncol=m_max) # all binary patterns for latent state profiles. 2^m_max of them.
  }


  # initialize the sampling chain: (improve this by warm start!)
  if (ALL_IN_ONE_START){
    t <- 1        # number of clusters.
    z <- rep(1,n) # z[i] is the cluster ID for observation i.
    mylist <- rep(0,t_max+3); mylist[1] <- 1  # mylist[1:t] is the list of active cluster IDs.
    # mylist is maintained in increasing order for 1:t, and
    # is 0 after that.
    c_next <- 2                    # an available cluster ID to be used next.
    N <- rep(0,t_max+3); N[1] <- n # N[c] is the size of cluster c. Note that N is different from n.

    log_p <- rep(0,n+1) # probability for sampling z; n+1 because at most there could be n clusters;
    # and sometimes one needs to assign an observation to a new cluster - by current bookkeeping,
    # we are not efficient in memory and create a new cluster ID - hence +1.
    zs <- rep(1,n)  # temporary cluster indicators for split-merge assignments.
    S  <- rep(0,n)  # temporary variable for indicies used during split-merge step.
  } else{
    hc <- hclust(dist(dat),"complete")
    t  <- floor(t_max/4)
    z  <- cutree(hc,k = t)
    mylist <- rep(0,t_max+3); mylist[1:t] <- 1:t  # mylist[1:t] is the list of active cluster IDs.
    c_next <- t+1
    N <- rep(0,t_max+3); N[1:t] <- table(z)
    log_p <- rep(0,n+1)
    zs <- z
    S  <- rep(0,n)
  }



  log_Nb <- log(1:n)+b # the multipliers needed when assigning an observation to an existing cluster. Restaurant process stuff.

  # variables for samples kept:
  t_samp <- rep(0,n_keep)
  N_samp <- matrix(0,nrow=t_max+3,ncol=n_keep)
  z_samp <- matrix(0,nrow=n,ncol=n_keep) # posterior samples of cluster indicators.
  H_star_samp <- array(0,c(t_max+3,m_max,n_keep))# component-specific machine profiles.
  mylist_samp <- matrix(0,nrow=length(mylist),ncol=n_keep) # need it to retrieve the correct latent state vector.

  if (is.null(model_options$theta)){
    theta_samp <- matrix(0,nrow=L,ncol=n_keep)
    a_theta    <- model_options$a_theta # <--------------------------- can be changed into a 2 by L matrix.
    theta   <- sapply(1:L,function(i){stats::rbeta(1,a_theta[1,i],a_theta[2,i])}) # initialization.
  } else{
    theta   <- model_options$theta # use specified true positive rates.
  }
  if (is.null(model_options$psi)){
    psi_samp   <- matrix(0,nrow=L,ncol=n_keep)
    a_psi      <- model_options$a_psi  # <--------------------------- can be changed into a 2 by L matrix.
    psi <- sapply(1:L,function(i){stats::rbeta(1,a_psi[1,i],a_psi[2,i])}) # initialization.
  } else{
    psi     <- model_options$psi # use specified false positive rates.
  }

  if (is.null(model_options$p0)){
    p_samp <- matrix(0,nrow=m_max,ncol=n_keep)
    p      <- rep(0.5,m_max) # initialization.
  } else{
    p      <- model_options$p0 # fix prevalence.
  }

  if (is.null(model_options$alpha) && is.null(model_options$p0)){
    alpha_samp <- rep(0,n_keep)  # hyperparameter for machine prevalence.
    alpha      <- m_max           # initialization.
  }else{
    alpha    <- model_options$alpha # fix alpha.
  }

  if (!do_update_partition){
    z <- model_options$the_partition #<-- must be from 1 to the total number of clusters.
    t <- length(unique(z))
    # recommend using the best scientific cluster indicators
    # after an initial fit without knowing the partition.
    mylist[1:t] <- 1:t
  }

  if (!is.null(mcmc_options$predictive_samp) && mcmc_options$predictive_samp){
    Y_pred_samp <- array(NA,c(n,L,n_keep))
  }
  cat("[rewind] Begin MCMC for model with #latent states (M =", m_max, ")\n")
  if(!is.null(model_options$Q) && is_identity_Q){cat("> identity Q.\n")}
  for (iter in 1:n_total){
    VERBOSE <- iter%%mcmc_options$print_mod==0
    if (!do_update_partition){
      N <- sapply(1:t,function(uu){sum(z==uu)})
    }else{
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

        if (is.null(is_identity_Q)){
          for (j in 1:t){
            cc       <- mylist[j]
            log_p[j] <- log_Nb[N[cc]]+log_marginal(rbind(dat[(z==cc)[-i],,drop=FALSE],dat[i,]),H_star_enumerate,Q,p,theta,psi)-
              log_marginal(dat[(z==cc)[-i],,drop=FALSE],H_star_enumerate,Q,p, theta,psi) # existing cluster.
          }
          log_p[t+1] <- log_v[t+1]-log_v[t] + log(b) + log_marginal(dat[i,,drop=FALSE],H_star_enumerate,Q,p,theta,psi) # new cluster.
        } else{
          if (is_identity_Q){
            for (j in 1:t){
              cc       <- mylist[j]
              log_p[j] <- log_Nb[N[cc]]+log_marginal_Q_identity(rbind(dat[(z==cc)[-i],,drop=FALSE],dat[i,]),p,theta,psi)-
                log_marginal_Q_identity(dat[(z==cc)[-i],,drop=FALSE],p, theta,psi) # existing cluster.
            }
            log_p[t+1] <- log_v[t+1]-log_v[t] + log(b) + log_marginal_Q_identity(dat[i,,drop=FALSE],p,theta,psi) # new cluster.
          }
        }
        j <- sample(t+1,1,prob = exp(log_p[1:(t+1)]))

        # add obs i to its new cluster:
        if (j<=t){
          c <- mylist[j]
        } else{
          c      <- c_prop
          mylist <- ordered_insert(c,mylist,t)
          t      <- t+1
          c_next <- ordered_next(mylist)
          if(t>t_max) {stop("[rewind] Sampled t has exceeded t_max. Increast t_max and retry.")}
        }
        z[i] <- c
        N[c] <- N[c] + 1

      }# END iteration over subjects.
    }# END IFELSE for doing partition update or not.

    if (do_update_partition){
      # if (use_splitmerge){}
      MORE_SPLIT <- mcmc_options$MORE_SPLIT
      res_split_merge <- split_merge(dat,z,zs,S,mylist,N,t,b,log_v,n_split,Q,p,theta,psi,MORE_SPLIT)
      t <- res_split_merge$t
      z <- res_split_merge$z
      N <- res_split_merge$N
      mylist <- res_split_merge$mylist
      c_next <- ordered_next(mylist)
    }
    if(t>t_max) {stop("[rewind] Sampled t has exceeded t_max. Increast t_max and retry.")}

    ## Simulate xi for posterior predictive checking:
    # #if (is.null(model_options$Q)){
    # xi_samp <- array(0,c(t_max+3,L,n_total)) # add this to the results.
    # for (j in 1:t){
    #   xi_samp[mylist[j],,iter]  <- update_xi(dat[z==mylist[j],,drop=FALSE],p,theta,psi)
    # }
    # xi_star <- xi_samp[mylist[1:t],,iter]
    # #}

    # update latent state profile given clusters: (check here for sampling H)
    H_star_redun <- matrix(0,nrow=t_max+3,ncol=m_max) # prior to merging redundant rows or columns.
    for (j in 1:t){
      if (block_update_H){
        curr_pattern_log_p   <- log_full(dat[z==mylist[j],,drop=FALSE],H_star_enumerate,Q,p,theta,psi)
        curr_pattern_p       <- exp(curr_pattern_log_p-matrixStats::logSumExp(curr_pattern_log_p))
        curr_ind             <- sample(1:(2^m_max),1,prob=curr_pattern_p)
        H_star_redun[mylist[j],] <- H_star_enumerate[curr_ind,]
      } else{
        for (m in 1:m_max){
          H_star_redun[mylist[j],m] <- 0
          L0            <- log_full(dat[z==mylist[j],,drop=FALSE],
                                    H_star_redun[mylist[j],,drop=FALSE],Q,p,theta,psi)
          H_star_redun[mylist[j],m] <- 1
          L1            <- log_full(dat[z==mylist[j],,drop=FALSE],
                                    H_star_redun[mylist[j],,drop=FALSE],Q,p,theta,psi)
          curr_eta_p    <- exp(L1-matrixStats::logSumExp(c(L0,L1)))
          H_star_redun[mylist[j],m] <-stats::rbinom(1,1,curr_eta_p)
          #H_star_samp[mylist[j],m,iter] <- metrop_flip(H_star_samp[mylist[j],m,iter],curr_eta_p)
        }
      }
    }

    H_star <- H_star_redun[mylist[1:t],,drop=FALSE] # don't use Z to subset H_star!!! Do so for H_star_samp.

    # update Q:
    if (is.null(model_options$Q)){
      if (block_update_Q && constrained){
        stop("[rewind] 'block_update_Q' and 'constrained' in 'model_options' must NOT be TRUE at the same time. Set one to FALSE and retry.")}
      if (block_update_Q){
        Q <- update_Q_col_block(dat,Q,H_star,z,t,mylist,p,theta,psi)
      } else{
        Q <- update_Q(dat,Q,H_star,z,t,mylist,p,theta,psi,constrained)
      }
    }

    if (VERBOSE){
      cat("\n[rewind] iteration ", iter, ":\n");
      cat("------------ \n");
      if (do_update_partition && !is.null(MORE_SPLIT) && MORE_SPLIT){cat("> Biasing towards splitting.\n")}
      cat("> Latent state profiles for t=",t," pseudo-clusters of sizes: ",N[N!=0],"\n")
      cat("> Merging identical rows (pseudo-clusters) and columns (partner latent states):\n")
      cat(">> H^* Before:\n")
      print_mat <- H_star; rownames(print_mat) <- paste(c("pseudo-cluster",rep("",t-1)),1:t,sep=" ");
      colnames(print_mat) <- paste(c("factor(machine)",rep("",ncol(print_mat)-1)),1:ncol(print_mat),sep=" ")
      print(print_mat) # <-- removed all zero columns.
      if (is.null(model_options$Q)){
        merged_res <- merge_H_Q(H_star,1:t,t,Q,TRUE)
      }
      if (is.null(model_options$alpha) && is.null(model_options$p0)){
        cat("> Finite IBP hyperparameter: alpha = ", alpha,"\n")
      }
    }

    # update true/false positive rates - theta/psi:
    if (is.null(model_options$theta) && is.null(model_options$psi)){
      res_update_pr <- update_positive_rate(dat,H_star_redun[z,,drop=FALSE],Q,a_theta,a_psi)
      theta <- res_update_pr$theta
      psi   <- res_update_pr$psi
    }

    # update hyperparameter for {p_m}:
    if (is.null(model_options$alpha) && is.null(model_options$p0)){
      alpha <- update_alpha(H_star,t,m_max)
    }

    # update prevalence vector {p_m}:
    if (is.null(model_options$p0)){
      p    <- update_prevalence(H_star,alpha,m_max)
    }

    #
    # Record MCMC results:
    #
    if (iter==keepers[keep_index+1]){
      keep_index   <- keep_index +1
      t_samp[keep_index] <- t
      for (j in 1:t){N_samp[mylist[j],keep_index] <- N[mylist[j]]}
      if (is.null(model_options$p0)){
        p_samp[,keep_index]    <- p
      }
      H_star_samp[mylist[1:t],,keep_index] <- H_star # <-- okay, insert H_star at places indicated by mylist.
      for (i in 1:n){z_samp[i,keep_index] <- z[i]}
      mylist_samp[,keep_index] <- mylist
      if (is.null(model_options$theta)){
        theta_samp[,keep_index] <- theta # of length L.
      }
      if (is.null(model_options$psi)){
        psi_samp[,keep_index]   <- psi
      }
      if (is.null(model_options$Q)){
        Q_samp[,,keep_index]       <- Q
      }
      if (is.null(model_options$alpha) && is.null(model_options$p0)){
        alpha_samp[keep_index] <- alpha
      }
      if (!is.null(mcmc_options$predictive_samp) && mcmc_options$predictive_samp){ # <-- slow.
        Y_pred <- (H_star_redun[z,]%*%Q>0.5)+0
        for (i_pred in 1:n){
          for (l_pred in 1:L){
            curr_p_pred <- c(theta[l_pred],psi[l_pred])[2-Y_pred[i_pred,l_pred]]
            Y_pred[i_pred,l_pred] <- rbinom(1,1,curr_p_pred)
          }
        }
        Y_pred_samp[,,keep_index] <- Y_pred
      }
    }

  }# END mcmc iterations.

  res <- list(keepers=keepers,t_samp=t_samp,N_samp=N_samp,z_samp=z_samp,
              H_star_samp = H_star_samp,
              mylist_samp=mylist_samp,keepers=keepers)
  if (is.null(model_options$theta)){
    res$theta_samp <- theta_samp # of length L.
  }
  if (is.null(model_options$psi)){
    res$psi_samp   <- psi_samp
  }
  if (is.null(model_options$Q)){
    res$Q_samp  <- Q_samp
  }
  if (is.null(model_options$p0)){
    res$p_samp <- p_samp
  }
  if (is.null(model_options$alpha) && is.null(model_options$p0)){
    res$alpha_samp <- alpha_samp
  }
  if (!is.null(mcmc_options$predictive_samp) && mcmc_options$predictive_samp){
    res$Y_pred_samp <- Y_pred_samp
  }
  res
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


#' compute the deriviative of the log of density of inactive states'
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
