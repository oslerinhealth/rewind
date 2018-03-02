rm(list=ls())
library(rewind)
library(matrixStats)
library(ars)
# ---------------------------------------------
# SIMULATIONS
#

# simulate data:
L0 <- 100
options_sim0  <- list(N = 200,  # sample size.
                      M = 3,   # true number of machines.
                      L = L0,   # number of antibody landmarks.
                      K = 8,    # number of true components.,
                      theta = rep(0.99,L0), # true positive rates
                      psi   = rep(0.01,L0), # false positive rates
                      alpha1 = 1 # half of the people have the first machine.
)

image(simulate_data(options_sim0,SETSEED = TRUE)$datmat)
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
  m_max  = 100,
  b  = 1, # Dirichlet hyperparameter; in the functions above, we used "b" - also can be called "gamma".
  #Q  = simu$Q,
  a_theta = c(99,1),
  a_psi   = c(1,99),
  a_alpha = 1, # hyperparameter for IBP alpha.
  b_alpha = 1,
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
  print_mode = 10,
  block_update_H = !TRUE,# NULL is also okay - because truncation level m_both can
  # be very large, using block updating is very expensive (2^m_both
  # multivariate binary patterns)
  block_update_Q = !TRUE # NULL is also okay - slice_sampler() do not use block updating
  # of Q for the same reason as mentioned above.
)


#
# run posterior algorithm for simulated data:
#

out <- slice_sampler(simu_dat,model_options0,mcmc_options0)

