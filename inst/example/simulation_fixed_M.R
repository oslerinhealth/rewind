\dontrun{
rm(list=ls())
library(rewind)
library(matrixStats)
library(ars)

#optional packages:
library(mcclust.ext) # http://wrap.warwick.ac.uk/71934/1/mcclust.ext_1.0.tar.gz

# color palette for heatmaps:
YlGnBu5   <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58',"#092d94")
hmcols    <- colorRampPalette(YlGnBu5)(256)

# simulate data:
L0 <- 100
options_sim0  <- list(N = 200,  # sample size.
                      M = 3,    # true number of machines.
                      L = L0,   # number of antibody landmarks.
                      K = 2^3,    # number of true components.
                      theta = rep(0.8,L0), # true positive rates.
                      psi   = rep(0.2,L0), # false positive rates.
                      alpha1 = 1 # half of the people have the first machine.
)

image(simulate_data(options_sim0,SETSEED = TRUE)$datmat,col=hmcols)
simu     <- simulate_data(options_sim0, SETSEED=TRUE)
simu_dat <- simu$datmat

#
# specifying options:
#

# model options:
model_options0 <- list(
  n   = nrow(simu_dat),
  t_max  = 40,
  m_max  = 5,
  b  = 1, # Dirichlet hyperparameter; in the functions above,
  # we used "b" - also can be called "gamma"!.
  #Q  = simu$Q,
  a_theta = replicate(L0, c(9,1)),
  a_psi   = replicate(L0, c(1,9)),
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
  n_keep  = 100,
  n_split = 5,
  print_mod = 10,
  constrained = TRUE, # <-- need to write a manual about when these options are okay.
  block_update_H = TRUE,
  block_update_Q = !TRUE,
  ALL_IN_ONE_START =!TRUE,  # <--- TRUE for putting all subjects in one cluster,
                            # FALSE by starting from a hierechical clustering
                            # (complete linkage) and cut
                    # to produce floor(t_max/4). Consider this as a warm start.
  hmcols = hmcols
)

#
# run posterior algorithm for simulated data:
#
out <- sampler(simu_dat,model_options0,mcmc_options0)


###############################################################
## Posterior summaries:
###############################################################
out0 <- out
out <- postprocess_H_Q(out0)

#Z_SAMP_FOR_PLOT <- out$z_samp  # <---- use pseudo-indicators for co-clustering.
                                # tend to be more granular.
Z_SAMP_FOR_PLOT <- out$z_sci_samp # <--- use scientific-cluster indicators.

# posterior co-clustering probabilities (N by N):
comat <- coclust_mat(nrow(simu_dat),Z_SAMP_FOR_PLOT)
image(1:options_sim0$N,1:options_sim0$N, comat,
      xlab="Subjects",ylab="Subjects",col=hmcols,main="co-clustering")
for (k in 1:options_sim0$K){
  abline(h=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
  abline(v=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
}

#
# Summarize partition C:
## Approach 0: The maximum a posterior: issue - the space is huge (Bell number).
##             Not implemented.
## Approach 1: Minimizing the least squared error between the co-clustering indicators
##            and the posterior co-clustering probabilities; (Dahl 2006)
## Approach 2: Wade - estimate the best partition using posterior expected loss
##             by variation of information (VI) metric.

nsamp_C <- dim(Z_SAMP_FOR_PLOT)[2]
z_Em    <- rep(NA,nsamp_C)

par(mfrow=c(2,3))

## Approach 1:
for (iter in 1:nsamp_C){
  z_Em[iter] <- norm(z2comat(Z_SAMP_FOR_PLOT[,iter])-comat,"F")
}
ind_dahl <- which.min(z_Em)
# plot 1:
image(1:options_sim0$N,1:options_sim0$N,
      z2comat(Z_SAMP_FOR_PLOT[,ind_dahl]),col=hmcols,
      main="Dahl (2006) - least squares to hat{pi}_ij")
for (k in 1:options_sim0$K){
  abline(h=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
  abline(v=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
}
#heatmap(z2comat(out$z_samp[,ind_dahl]),col=hmcols)

## Approach 2 - use "mcclust.ext" pacakge:
# process the posterior samples of cluster indicators:
psm    <- comp.psm(t(Z_SAMP_FOR_PLOT))
# point estimate using all methods:
bmf.VI <- minVI(psm,t(Z_SAMP_FOR_PLOT),method="all",include.greedy=TRUE)
summary(bmf.VI)

# plot 2:
image(1:options_sim0$N,1:options_sim0$N,
      z2comat(bmf.VI$cl["avg",]),col=hmcols,
      main="Wade clustering")

for (k in 1:options_sim0$K){
  abline(h=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
  abline(v=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
}
#heatmap(z2comat(bmf.VI$cl["avg",]),col=hmcols)

# credible sets as defined in Wade and Ghahramani 2017.
bmf.cb <- credibleball(bmf.VI$cl[1,],t(Z_SAMP_FOR_PLOT))

# plot 3:
image(1:options_sim0$N,1:options_sim0$N,
      z2comat(bmf.cb$c.horiz[1,]),col=hmcols,
      main="Wade credible ball - horizontal")

for (k in 1:options_sim0$K){
  abline(h=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
  abline(v=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
}
#heatmap(z2comat(bmf.cb$c.horiz),col=hmcols)

# plot 4:
image(1:options_sim0$N,1:options_sim0$N,
      z2comat(bmf.cb$c.uppervert[1,]),col=hmcols,
      main="Wade credible ball - upper vertical")

for (k in 1:options_sim0$K){
  abline(h=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
  abline(v=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
}
#heatmap(z2comat(bmf.cb$c.uppervert),col=hmcols)

# plot 5:
image(1:options_sim0$N,1:options_sim0$N,
      z2comat(bmf.cb$c.lowervert),col=hmcols,
      main="Wade credible ball - lower vertical")

for (k in 1:options_sim0$K){
  abline(h=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
  abline(v=cumsum(rle(simu$Z)$lengths)[k]+0.5,lty=2)
}
#heatmap(z2comat(bmf.cb$c.lowervert),col=hmcols)

#
# summarize Q (first need to transpose it)
#
# Approach 1: compute the error |QQ' - E[QQ'|Y]|_Frobneious
# Approach 2: compute the marginal co-activation probability
#             P(\sum_m Q_ml >= 1, \sum_m Q_ml' >= 1) -
#             This is invariant to the relabeling of latent states, or cluster labels.

nsamp_Q <- dim(out$Q_merge_samp)[3]
EQQ     <- matrix(0,nrow=dim(out$Q_merge_samp)[2],
                  ncol=dim(out$Q_merge_samp)[2]) # L by L.
EQQ_binary <- EQQ
# Approach 1:
for (iter in 1:nsamp_Q){
  A   <- t(out$Q_merge_samp[,,iter])
  EQQ <- (A%*%t(A)+EQQ*(iter-1))/iter # invariant to the row reordering of Q.
  EQQ_binary <- 0+(A%*%t(A)>0)+EQQ_binary
}

# Approach 2:
EQQ_percent <- EQQ_binary/nsamp_Q
image(EQQ_percent,col=hmcols,main="co-activation patterns across L dimensions")

# Approach 1:
Q_Em <- rep(NA,nsamp_Q)
for (iter in 1:nsamp_Q){
  A <- t(out$Q_merge_samp[,,iter])
  Q_Em[iter] <- norm(A%*%t(A) - EQQ,"F")
}
plot(Q_Em,type="o",main="||QQ'-E[QQ'|Y]||")

# choose the indices minimizing the errors:
ind_of_Q       <- which(Q_Em==min(Q_Em))

#
# visualize truth:
#
pdf(file.path("inst/example_figure/","bmf_truth.pdf"),width=12,height=6)
plot_truth(simu,options_sim0)
dev.off()

#
# visualize the individual latent states depending on whether Q is known or not.
#
if (!is.null(model_options0$Q)){ # <--- Q known.
  H_res     <- matrix(0,nrow=nrow(simu$datmat),ncol=sum(rowSums(model_options0$Q)!=0))
  H_pat_res <- matrix(0,nrow=nrow(simu$datmat),ncol=dim(out$H_star_merge_samp)[3])
  for (l in 1:dim(out$H_star_merge_samp)[3]){
    tmp_mylist <- out$mylist_samp[,l]
    tmp0 <- out$H_star_samp[tmp_mylist[match(out$z_samp[,l],tmp_mylist)],,l]
    #<------ could be simpler!!!
    tmp1 <- tmp0[,colSums(tmp0)!=0,drop=FALSE]
    tmp <- tmp1[,order_mat_byrow(model_options0$Q[rowSums(model_options0$Q)!=0,
                                                  ,drop=FALSE])$ord,drop=FALSE]
    H_pat_res[,l] <- bin2dec_vec(tmp,LOG=FALSE)
    H_res <- (tmp + H_res*(l-1))/l
  }
  apply(H_pat_res,1,table)
  image(f(H_res))
} else{   # <---- Q unknown.
  #
  # This is the best estimated Q:
  #
  Q_merged <- out$Q_merge_samp[,,ind_of_Q[1]] # just picked one.
  NROW_Q_PLOT <- nrow(Q_merged) # sum(rowSums(Q_merged)!=0)
  Q_PLOT <- f(order_mat_byrow(Q_merged)$res) # t(Q_merged)
  #f(order_mat_byrow(Q_merged[rowSums(Q_merged)!=0,,drop=FALSE])$res)
  image(1:ncol(simu$datmat),1:NROW_Q_PLOT,
        Q_PLOT,
        main="Best Q (merged & ordered)",col=hmcols,
        xlab="Dimension (1:L)",
        ylab="Latent State (1:M)",yaxt="n",cex.lab=1.2)
  for (k in 1:NROW_Q_PLOT){
    abline(h=NROW_Q_PLOT-k+0.5,
           lty=2,col="grey",lwd=2)
  }

  #
  # Summarize H* and H*(Z) for individual predictions:
  #
  # If we compare it to H, then the clustering matters too, not
  # just H^*.
  #
  # Approach: just obtain the H samples from those iterations with the best Q
  #            obtained above.
  #
  H_res     <- matrix(0,nrow=nrow(simu$datmat),ncol=sum(rowSums(Q_merged)!=0))
  H_pat_res <- matrix(0,nrow=nrow(simu$datmat),ncol=length(ind_of_Q))
  # columns for the best indices.
  # here ind_of_Q are those that minimized the Q loss.
  for (l in seq_along(ind_of_Q)){
    tmp_mylist <- out$mylist_samp[,ind_of_Q[l]]
    tmp0 <- out$col_merged_H_star_samp[out$z_samp[,ind_of_Q[l]],,ind_of_Q[l]]
    tmp1 <- tmp0[,colSums(tmp0)!=0,drop=FALSE]
    tmp <- tmp1[,order_mat_byrow(Q_merged[rowSums(Q_merged)!=0,,
                                          drop=FALSE])$ord,drop=FALSE]
    H_pat_res[,l] <- bin2dec_vec(tmp,LOG=FALSE)
    H_res <- (tmp + H_res*(l-1))/l
  }
  apply(H_pat_res,1,table)
  pdf(file.path("inst/example_figure/","H_estimated_marginal_prob_vs_truth.pdf"),width=12,height=6)
  par(mfrow=c(1,2))
  image(f(H_res),col=hmcols,main="estimated marginal prob") # <--- get marg prob.
  image(f(simu$Eta[,order_mat_byrow(simu$Q)$ord]),col=hmcols,main="truth")
  dev.off()
  # issues: the order of the rows of Q at ind_of_Q might be different, so need to order them.
}

# posterior distribution over the number of pseudo-clusters T: <-- scientific clusters?
plot(out$t_samp,type="l",ylab="T: #pseudo-clusters")

## individual predictions:
pdf(file.path("inst/example_figure/","individual_pred_simulation.pdf"),height=15,width=12)
par(mar=c(2,8,8,0),mfrow=c(2,1),oma=c(5,5,5,5))
for (i in 1:nrow(simu_dat)){
  plot_individual_pred(apply(H_pat_res,1,table)[[i]]/sum(apply(H_pat_res,1,table)[[i]]),
                       1:model_options0$m_max,paste0("Obs ", i),asp=0.5)
}
dev.off()

# check positive rate estimates:
pdf(file.path("inst/example_figure/","check_PR_post_vs_prior.pdf"),width=12,height=6)
par(mfrow=c(5,2))
for (l in 1:L0){
  par(mar=c(1,1,1,2))
  baker::beta_plot(model_options0$a_psi[1,l],model_options0$a_psi[2,l])
  hist(out$psi_samp[l,],add=TRUE,freq=FALSE)
  legend("topright",legend = l)
  baker::beta_plot(model_options0$a_theta[1,l],model_options0$a_theta[2,l])
  hist(out$theta_samp[l,],add=TRUE,freq=FALSE)
}
dev.off()


#
# posterior of scientific clusters:
#
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

pdf(file.path("inst/example_figure/","posterior_sci_cluster_number.pdf"),
    width=10,height=6)
scatterhist(1:ncol(out$z_sci_samp),apply(out$z_sci_samp,2,
                                         function(v){length(unique(v))}),
            "Index","tilde{T}: #scientific clusters")
dev.off()

# posterior distribution of active/effective dimensions (machines):
pdf(file.path("inst/example_figure/","posterior_effective_machines_number.pdf"),
    width=10,height=6)
plot(table(apply(out$col_merged_H_star_samp,c(3),
                 function(x) sum(colSums(x)>0)))/mcmc_options0$n_keep,
     xlab = "Effect Number of Machines",xlim=c(1,model_options0$m_max),
     ylab ="Posterior Probability",
     main="",lwd=4)
dev.off()

}




