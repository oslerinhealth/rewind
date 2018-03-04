rm(list=ls())
library(rewind)
library(matrixStats)
library(ars)

# color palette for heatmaps:
YlGnBu5   <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58',"#092d94")
hmcols    <- colorRampPalette(YlGnBu5)(256)

# simulate data:
L0 <- 100
options_sim0  <- list(N = 200,  # sample size.
                      M = 3,   # true number of machines.
                      L = L0,   # number of antibody landmarks.
                      K = 2^3,    # number of true components.,
                      theta = rep(0.8,L0), # true positive rates
                      psi   = rep(0.1,L0), # false positive rates
                      alpha1 = 1 # half of the people have the first machine.
)

image(simulate_data(options_sim0,SETSEED = TRUE)$datmat,hmcols)
simu     <- simulate_data(options_sim0, SETSEED=TRUE)
simu_dat <- simu$datmat

#
# specifying options
#

# check BayesianMixtures.jl for how options were set.

# model options:
model_options0 <- list(
  n   = nrow(simu_dat),
  t_max  = 40,
  m_max  = 10,
  b  = 1, # Dirichlet hyperparameter; in the functions above, we used "b" - also can be called "gamma".
  #Q  = simu$Q,
  a_theta = c(9,1),
  a_psi   = c(1,9),
  #theta = options_sim0$theta,
  #psi   = options_sim0$psi,
  #alpha   = options_sim0$M,
  #p_both      = rep(0.5,3),#,c(0.5,0.5^2,0.5^3,0.5^4,0.5^5)
  log_pk = "function(k) {log(0.1) + (k-1)*log(0.9)}"# Geometric(0.1). Prior for the number of components.
)

# pre-compute the log of coefficients in MFM:
model_options0$log_v <- mfm_coefficients(eval(parse(text=model_options0$log_pk)),
                                         model_options0$gamma,
                                         model_options0$n,
                                         model_options0$t_max+1)

# mcmc options:
mcmc_options0 <- list(
  n_total = 200,
  n_keep  = 200,
  n_split = 5,
  print_mod = 10,
  constrained = TRUE,
  block_update_H = TRUE,
  block_update_Q = !TRUE,
  hmcols =hmcols
)


#
# run posterior algorithm for simulated data:
#
out <- sampler(simu_dat,model_options0,mcmc_options0)


#
# Posterior summaries:
#

# comparing true Q to sampled Q:
ind_post_mode <- mcmc_options0$n_total#which.max(table(pats))
#pdf("diagnosticsQ.pdf",width=12,height=6)
par(mfrow=c(1,4))
image(simu$datmat,main="Data",col=hmcols)
image(simu$xi,main="True presence/absence of proteins)",col=hmcols)
image(order_mat_byrow(simu$Q)$res,main="True Q (ordered)",col=hmcols)
Q_merged <- out$Q_merge_samp[,,ind_post_mode]
image(order_mat_byrow(Q_merged[rowSums(Q_merged)!=0,,drop=FALSE])$res,main="Sampled Q (merged & ordered)",col=hmcols)
#dev.off()

# co-clustering:
comat <- coclust_mat(nrow(simu_dat),out$z_samp,mcmc_options0$n_keep)
image(1:options_sim0$N,1:options_sim0$N, comat,
      xlab="Subjects",ylab="Subjects",col=hmcols,main="co-clustering")
for (k in 1:options_sim0$K){
  abline(h=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
  abline(v=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
}

