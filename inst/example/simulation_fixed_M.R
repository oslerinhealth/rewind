#\dontrun{
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
options_sim0  <- list(N = 100,  # sample size.
                      M = 3,    # true number of machines.
                      L = L0,   # number of antibody landmarks.
                      K = 2^3,    # number of true components.
                      theta = rep(0.9,L0), # true positive rates.
                      psi   = rep(0.1,L0), # false positive rates.
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




################################################################
## Posterior summaries:
###############################################################

##
## Post-processing H an Q:
##
n_kept <- dim(out$H_star_samp)[3]
Q_merge_samp <- array(0,c(model_options0$m_max,dim(out$Q_samp)[2],n_kept))
H_star_merge_samp <- array(0,c(dim(out$H_star_samp)[1],model_options0$m_max,n_kept))
col_merged_H_star_samp <- H_star_merge_samp
z_sci_samp <- matrix(0,nrow=dim(out$z_samp)[1],ncol=n_kept)
for (iter in 1:n_kept){
  merged_res <- merge_H_Q(out$H_star_samp[,,iter],
                          out$mylist[,iter],
                          out$t_samp[iter],
                          out$Q_samp[,,iter],
                          FALSE,
                          out$z_samp[,iter])
  merged_H <- merged_res$H_star_merge
  merged_Q <- merged_res$Q_merge
  H_star_merge_samp[1:nrow(merged_H),1:ncol(merged_H),iter] <- merged_H
  Q_merge_samp[1:nrow(merged_Q),1:ncol(merged_Q),iter]      <- merged_Q
  z_sci_samp[,iter]    <- merged_res$z_sci
  col_merged_H_star_samp[,1:ncol(merged_H),iter] <- merge_H_col(out$H_star_samp[,,iter])$H_star_merge # <-- may still have extra zeros.
}
out$H_star_merge_samp <- H_star_merge_samp # <-- could have zero columns.
out$Q_merge_samp      <- Q_merge_samp # <-- could have zero rows.
out$z_sci_samp        <- z_sci_samp # <-- add scientific indicators.
out$col_merged_H_star_samp <- col_merged_H_star_samp

#Z_SAMP_FOR_PLOT <- out$z_samp  # <---- use pseudo-indicators for co-clustering.
                                # tend to be more granular.
Z_SAMP_FOR_PLOT <- out$z_sci_samp # <--- use scientific-cluster indicators.

# posterior co-clustering probabilities (N by N):
comat <- coclust_mat(nrow(simu_dat),Z_SAMP_FOR_PLOT,mcmc_options0$n_keep)
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
      z2comat(bmf.cb$c.horiz),col=hmcols,
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

f <- function(m) t(m)[,nrow(m):1]

#
# visualize truth: <------------------ IMPORTANT!
#
pdf("bmf_truth.pdf",width=12,height=6)
par(mfcol=c(4,1),tcl=-0.5,
    mai=c(0.7,0.7,0.3,0.3))

m <- cbind(c(1, 1), c(2, 2),c(3,4))
layout(m, widths=c(6,3,6), heights=c(2,4))
#layout.show(4)
par(mar = c(5, 5, 5, 1))

# image(1:ncol(simu$datmat),1:nrow(simu$datmat),
#       f(simu$datmat),main="Data",col=hmcols,
#       xlab="Dimension (1:L)",
#       ylab="Subject (1:N)",cex.lab=1.2)
# for (k in 1:options_sim0$K){
#   abline(h=100-cumsum(rle(simu$Z)$lengths)[k]+0.5,
#          lty=2,col="grey",lwd=2)
# }

image(1:ncol(simu$datmat),1:nrow(simu$datmat),
      f(simu$xi),#main="Design Matrix",
      col=hmcols,
      xlab="Dimension (1:L)",yaxt="n",
      ylab="Subject (1:N)",cex.lab=1.2)
axis(side=2,at=1:options_sim0$N,labels=rev(1:options_sim0$N),las=2,cex.axis=0.5)
mtext(expression(paste("(",Gamma[il],"): Design Matrix",collapse="")),3,1)
for (k in 1:options_sim0$K){
  abline(h=100-cumsum(rle(simu$Z)$lengths)[k]+0.5,
         lty=2,col="grey",lwd=2)
}

image(1:ncol(simu$Eta),1:nrow(simu$datmat),
      f(simu$Eta[,order_mat_byrow(simu$Q)$ord]),#main="Latent State Profile",
      col=hmcols,
      xlab="Latent State (1:M)",ylab="Subject (1:N)",
      yaxt="n",cex.lab=1.2,xaxt="n")
axis(side=2,at=1:options_sim0$N,labels=rev(1:options_sim0$N),las=2,cex.axis=0.5)
mtext(expression(paste("(",eta[im],"): Latent State",collapse="")),3,1)
axis(side = 1,at=1:options_sim0$M,labels=1:options_sim0$M)

for (k in 1:options_sim0$K){
  abline(h=100-cumsum(rle(simu$Z)$lengths)[k]+0.5,
         lty=2,col="grey",lwd=2)
}
for (k in 1:options_sim0$K){
  abline(v=options_sim0$M-k+0.5,
         lty=1,col="black",lwd=2)
}

image(1:ncol(simu$datmat),1:options_sim0$M,
      f(order_mat_byrow(simu$Q)$res),
      #main="True Q (ordered)",
      col=hmcols,
      xlab="Dimension (1:L)",
      ylab="Latent State (1:M)",yaxt="n",cex.lab=1.2)
mtext(expression(paste("(",Q[ml],"): True Q (ordered)",collapse="")),3,1)

axis(side = 2,at=1:options_sim0$M,labels=1:options_sim0$M)

for (k in 1:options_sim0$M){
  abline(h=options_sim0$M-k+0.5,
         lty=1,col="black",lwd=2)
}

plot.new()
legend("center",legend = c(1,0),col=c("#092D94","#FFFFD9"),
       pch=c(15,15),cex=3,bty="n")
dev.off()


#
# visualize the individual latent states depending on whether Q is known or not.
#
if (!is.null(model_options0$Q)){
  H_res     <- matrix(0,nrow=nrow(simu$datmat),ncol=sum(rowSums(model_options0$Q)!=0))
  H_pat_res <- matrix(0,nrow=nrow(simu$datmat),ncol=dim(out$H_star_merge_samp)[3])
  for (l in 1:dim(out$H_star_merge_samp)[3]){
    tmp_mylist <- out$mylist_samp[,l]
    tmp0 <- out$H_star_samp[tmp_mylist[match(out$z_samp[,l],tmp_mylist)],,l] #<------ wrong!!!
    tmp1 <- tmp0[,colSums(tmp0)!=0,drop=FALSE]
    tmp <- tmp1[,order_mat_byrow(model_options0$Q[rowSums(model_options0$Q)!=0,,drop=FALSE])$ord,drop=FALSE]
    H_pat_res[,l] <- bin2dec_vec(tmp,LOG=FALSE)
    H_res <- (tmp + H_res*(l-1))/l
  }
  apply(H_pat_res,1,table)
  image(H_res)
} else{
  #
  # This is the best estimated Q:
  #
  Q_merged <- out$Q_merge_samp[,,ind_of_Q[1]] # just picked one.
  NROW_Q_PLOT <- nrow(Q_merged) # sum(rowSums(Q_merged)!=0)
  Q_PLOT <- f(order_mat_byrow(Q_merged)$res) #f(order_mat_byrow(Q_merged[rowSums(Q_merged)!=0,,drop=FALSE])$res)
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
  H_pat_res <- matrix(0,nrow=nrow(simu$datmat),ncol=length(ind_of_Q)) # columns for the best indices.
  # here ind_of_Q are those that minimized the Q loss.
  for (l in seq_along(ind_of_Q)){
    tmp_mylist <- out$mylist_samp[,ind_of_Q[l]]
    tmp0 <- out$col_merged_H_star_samp[tmp_mylist[match(out$z_samp[,ind_of_Q[l]],tmp_mylist)],,ind_of_Q[l]]
    tmp1 <- tmp0[,colSums(tmp0)!=0,drop=FALSE]
    tmp <- tmp1[,order_mat_byrow(Q_merged[rowSums(Q_merged)!=0,,drop=FALSE])$ord,drop=FALSE]
    H_pat_res[,l] <- bin2dec_vec(tmp,LOG=FALSE)
    H_res <- (tmp + H_res*(l-1))/l
  }
  apply(H_pat_res,1,table)
  image(f(H_res),col=hmcols)

  # issues: the order of the rows of Q at ind_of_Q might be different, so need to order them.
}


# population quantities: <---------------- important in childhood pneumonia examples.
compute_table <- function(p){
  M <- length(p)
  pat <-  as.matrix(expand.grid(rep(list(0:1),M)),ncol=M)
  res <- rowSums(pat%*%diag(log(p))+(1-pat)%*%diag(log(1-p)))
  names(res) <- apply(pat,1,paste,collapse="")
  res
}
log_eti_samp <- apply(out$p_samp,2,compute_table)
eti_samp_mean <- rowMeans(exp(log_eti_samp))
plot(eti_samp_mean,type="h")

pat <-  as.matrix(expand.grid(rep(list(0:1),model_options0$m_max)),
                  ncol=model_options0$m_max)

pdf("population_fraction.pdf",height=6,width=8)
plot(rep(1:nrow(pat), each = ncol(pat)),
     rep(-rev(1:ncol(pat)), nrow(pat)),
     axes = FALSE, ann = FALSE,
     pch = ifelse(t(pat), 19, 1), cex = 2.5, asp = 1, xpd = NA,
     col = rep(c("grey","black")[1+(eti_samp_mean>quantile(eti_samp_mean,0.75))],
               each=ncol(pat)))

multiplier <- 4
segments(1:nrow(pat),exp(apply(log_eti_samp,1,quantile,0.025))*ncol(pat)*multiplier,
         1:nrow(pat),exp(apply(log_eti_samp,1,quantile,0.975))*ncol(pat)*multiplier)
segments(1:nrow(pat),exp(apply(log_eti_samp,1,quantile,0.25))*ncol(pat)*multiplier,
         1:nrow(pat),exp(apply(log_eti_samp,1,quantile,0.75))*ncol(pat)*multiplier,lwd=4,col="dodgerblue2")
points(1:nrow(pat),rowMeans(exp(log_eti_samp))*ncol(pat)*multiplier,pch=18,xpd = NA)
points(1:nrow(pat),exp(apply(log_eti_samp,1,quantile,0.025))*ncol(pat)*multiplier,lwd=6,pch="-")
points(1:nrow(pat),exp(apply(log_eti_samp,1,quantile,0.975))*ncol(pat)*multiplier,lwd=6,pch="-")
axis(2,at=ncol(pat)*c(0,0.25,0.5)*multiplier,labels=c(0,0.25,0.5),las=2,xpd = NA)

mtext(text = expression(paste("", pi[eta],
                              sep="")),line = 2,side=2,cex=2,las=2,padj=-4)
mtext(text = expression(paste(eta[1],",..., ", eta[M],collapse="")),line = 1,side=2,
      cex=1,adj = 0.5,las=0)

legend("topright",legend = c("95% CI","50% CI"),col=c("black","dodgerblue2"),lwd=c(1,4),
       bty="n")

#
# This is the best estimated Q:
#
Q_merged <- out$Q_merge_samp[,,ind_of_Q[1]] # just picked one.
NROW_Q_PLOT <- nrow(Q_merged) #sum(rowSums(Q_merged)!=0)
Q_PLOT <- t(order_mat_byrow(Q_merged)$res)
#<- f(order_mat_byrow(Q_merged[rowSums(Q_merged)!=0,,drop=FALSE])$res)
image(seq(1,nrow(pat),length=ncol(simu$datmat)),-rev(1:NROW_Q_PLOT)-NROW_Q_PLOT-1,
      Q_PLOT,
      main="Best Q (merged & ordered)",col=hmcols,
      xlab="Dimension (1:L)",
      ylab="Latent State (1:M)",yaxt="n",cex.lab=1.2,add=TRUE)
for (k in 1:(NROW_Q_PLOT+1)){
  abline(h=-k+0.5-NROW_Q_PLOT-1,
         lty=2,col="grey",lwd=2)
}

mtext(text = expression(paste(Q[1.],",..., ", Q[M.],collapse="")),line = 1,side=2,
      cex=1,adj = 0.15,las=0)
dev.off()

# posterior distribution over the number of pseudo-clusters T: <-- scientific clusters?
plot(out$t_samp,type="l",ylab="T: #pseudo-clusters")

## individual predictions:
pdf("individual_pred_simulation.pdf",height=15,width=12)
par(mar=c(2,8,8,0),mfrow=c(2,1),oma=c(5,5,5,5))
for (i in 1:nrow(simu_dat)){
  plot_individual_pred(apply(H_pat_res,1,table)[[i]]/sum(apply(H_pat_res,1,table)[[i]]),
                       1:model_options0$m_max,
                       simu_dat[i,,drop=FALSE],asp=0.5)
}
dev.off()






#}




