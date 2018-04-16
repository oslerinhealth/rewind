#\dontrun{
rm(list=ls())
library(rewind)
library(matrixStats)
library(ars)

# color palette for heatmaps:
YlGnBu5   <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58',"#092d94")
hmcols    <- colorRampPalette(YlGnBu5)(256)

# simulate data:
L0 <- 100
options_sim0  <- list(N = 100,  # sample size.
                      M = 3,   # true number of machines.
                      L = L0,   # number of antibody landmarks.
                      K = 2^3,    # number of true components.,
                      theta = rep(0.8,L0), # true positive rates
                      psi   = rep(0.1,L0), # false positive rates
                      alpha1 = 1 # half of the people have the first machine.
)

image(simulate_data(options_sim0,SETSEED = TRUE)$datmat,col=hmcols)
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
  m_max  = 5,
  b  = 1, # Dirichlet hyperparameter; in the functions above,
  # we used "b" - also can be called "gamma"!.
  #Q  = simu$Q,
  a_theta = c(9,1),
  a_psi   = c(1,9),
  #theta = options_sim0$theta,
  #psi   = options_sim0$psi,
  #alpha   = options_sim0$M,
  #p_both      = rep(0.5,3),#,c(0.5,0.5^2,0.5^3,0.5^4,0.5^5)
  log_pk = "function(k) {log(0.1) + (k-1)*log(0.9)}"# Geometric(0.1).
  #Prior for the number of components.
)



# pre-compute the log of coefficients in MFM:
model_options0$log_v<-mfm_coefficients(eval(parse(text=model_options0$log_pk)),
                                       model_options0$b,
                                       model_options0$n,
                                       model_options0$t_max+1)

# mcmc options:
mcmc_options0 <- list(
  n_total = 200,
  n_keep  = 200,
  n_split = 5,
  print_mod = 10,
  constrained = TRUE, # <-- need to write a manual about when these options are okay.
  block_update_H = TRUE,
  block_update_Q = !TRUE,
  hmcols = hmcols
)


#
# run posterior algorithm for simulated data:
#
out <- sampler(simu_dat,model_options0,mcmc_options0)

#
# Posterior summaries: (Question, how to reconcile differences of the minimizing
#    elements for, H, Q, Z?)
#

# co-clustering:
comat <- coclust_mat(nrow(simu_dat),out$z_samp,mcmc_options0$n_keep)
image(1:options_sim0$N,1:options_sim0$N, comat,
      xlab="Subjects",ylab="Subjects",col=hmcols,main="co-clustering")
for (k in 1:options_sim0$K){
  abline(h=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
  abline(v=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
}

trans_squared_error <- function(A,B){
  norm(A%*%t(A) - B%*%t(B),"F")
}

myentropy <- function(z){
  n <- length(z)
  tb <- table(z)
  log(n) - sum(tb*log(tb))/n
  }

z2comat <- function(z){
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

## summarize partition C:
## Approach 0: the maximum a posterior: issue - the space is huge.
## Approach 1: using the entropy to get posterior median and the credible interval;
## Approach 2: minimizing the least squared error between the co-clustering indicators
##            and the posterior co-clustering probabilities;
## Approach 3: Wade - estimate (greedy) the best partition using posterior expected loss (VI).

nsamp_C <- dim(out$z_samp)[2]
z_Em <- rep(NA,nsamp_C)
## approach 2:
for (iter in 1:nsamp_C){
  z_Em[iter] <- norm(z2comat(out$z_samp[,iter])-comat,"F")
}
plot(z_Em,type="o",main="|S-E[S|Y]|")
image(1:options_sim0$N,1:options_sim0$N,
      z2comat(out$z_samp[,which.min(z_Em)]),col=hmcols,
      main="clustering at the minimized Z")
for (k in 1:options_sim0$K){
  abline(h=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
  abline(v=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
}

#3 approach 1:
for (iter in 1:nsamp_C){
  z_Em[iter] <- myentropy(out$z_samp[,iter])
  }

bd <- HDInterval::hdi(z_Em)
z_Em[which(z_Em==bd[1])]
z_Em[which(z_Em==bd[2])]
which(z_Em == median(z_Em))

# lower:
image(1:options_sim0$N,1:options_sim0$N,
      z2comat(out$z_samp[,which(z_Em==bd[1])[1]]),col=hmcols,
      main="clustering at the lower entropy quantile")

# upper:
image(1:options_sim0$N,1:options_sim0$N,
      z2comat(out$z_samp[,which(z_Em==bd[2])]),col=hmcols,
      main="clustering at the upper entropy quantile")

# median:
image(1:options_sim0$N,1:options_sim0$N,
      z2comat(out$z_samp[,which(z_Em == median(z_Em))[3]]),col=hmcols,
      main="clustering at the minimized Z")

for (k in 1:options_sim0$K){
  abline(h=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
  abline(v=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
}

## summarize H*:
##
## If we compare it to H, then the clustering matters too, not
## just H^*.
##
##  compute the error |HH' - E[HH'|Y]|_Frobneious
##
nsamp_H <- dim(out$z_samp)[2]
H_Em <- rep(NA,nsamp_H)
for (iter in 1:nsamp_H){
  H_Em[iter] <- trans_squared_error(out$H_star_merge_samp[out$z_samp[,iter],,iter],
                                    simu$H_star[simu$Z,])
}
plot(H_Em,type="o",main="|HH'-E[H0H0'|Y]|")

min_H_Em <- which.min(H_Em) # select the one with minimal squared error.
image(1:options_sim0$N,1:model_options0$m_max,out$H_star_merge_samp[out$z_samp[,min_H_Em],,min_H_Em],
      col=hmcols,xlab="subjects",ylab="latent states")
for (k in 1:options_sim0$K){
  abline(v=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
}

## summarize Q (first need to transpose it)
##
##
## Approach 1:  compute the error |QQ' - E[QQ'|Y]|_Frobneious
## Approach 2: compute the marginal co-positive probability
##            P(\sum_m Q_ml >= 1, \sum_m Q_ml' >= 1) -
##             This is invariant to the relabeling of latent states, or cluster labels.
nsamp_Q <- dim(out$Q_merge_samp)[3]
Q_Em <- rep(NA,nsamp_Q)
for (iter in 1:nsamp_Q){
  Q_Em[iter] <- trans_squared_error(t(out$Q_merge_samp[,,iter]), t(simu$Q))
}
plot(Q_Em,type="o",main="|QQ'-E[Q0Q0'|Y]|")

# comparing true Q to sampled Q:
ind_post_mode <- mcmc_options0$n_total#which.max(table(pats))
#pdf("diagnosticsQ.pdf",width=12,height=6)
par(mfrow=c(1,4))
image(simu$datmat,main="Data",col=hmcols)
image(simu$xi,main="True presence/absence of proteins)",col=hmcols)
image(order_mat_byrow(simu$Q)$res,main="True Q (ordered)",col=hmcols)
Q_merged <- out$Q_merge_samp[,,ind_post_mode]
image(order_mat_byrow(Q_merged[rowSums(Q_merged)!=0,,drop=FALSE])$res,
      main="Sampled Q (merged & ordered)",col=hmcols)
#dev.off()



#}
