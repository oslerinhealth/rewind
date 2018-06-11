#' simulate data from boolean matrix factorization model
#'
#' This function takes in options_sim and simulate as specified by the options.
#'
#' @param options_sim options for simulating data: N - sample size; M - number of machines;
#' L - number of protein landmarks; K - number of ture components; alpha1 - hyperparameter
#' for determining the machine prevalence {p_m}; theta - true positive rates for
#' all protein landmarks; psi - false positive rates for all protein landmarks.
#'
#' @param SETSEED FALSE by default; Set to TRUE if to fix random number generation.
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
#' image(simulate_data(options_sim0,SETSEED = TRUE)$datmat)
#' simu     <- simulate_data(options_sim0, SETSEED=TRUE)
#' simu_dat <- simu$datmat
#' @return \itemize{
#' \item \code{datmat} a matrix for multivariate binary observations from the assumed boolean
#' matrix factorization,
#' \item \code{Q} A Q-matrix,
#' \item \code{H_star} A binary matrix of dimension K by M, for K clusters, and M factors,
#' \item \code{Z} A vector of length N for individual cluster indicators,
#' \item \code{xi} a binary matrix of dimension N by L for true presence or absence of proteins
#' \item \code{Eta} A binary matrix of dimension N by M for presence or absence of factors among
#' all subejcts.
#' }
#' @export
simulate_data <- function(options_sim, SETSEED=FALSE){
  if(SETSEED){set.seed(2018-01-29)}
  N <- options_sim$N
  M <- options_sim$M
  L <- options_sim$L
  K <- options_sim$K
  alpha1 <- options_sim$alpha1
  theta  <- options_sim$theta
  psi    <- options_sim$psi

  # Q matrix: Check Chen 2015 JASA for sufficient condition for identification:
  Q <- simulate_Q(M,L)
  # # This is for random eta_atom: for testing, we now use eta_atom as prespecified.
  # eta_atom <- matrix(0,nrow=K,ncol=M)
  # for (k in 1:K){
  #   V <- rep(0,M)
  #   p <- rep(0,M)
  #   for (m in 1:M){
  #     V[m] <- rbeta(1,alpha1,1)
  #     p[m] <- exp(sum(log(V[1:m])))
  #     eta_atom[k,m] <- rbinom(1,1,p[m])
  #   }
  # }

  if (K==1){
    eta_atom <- matrix(c(1,0,1),nrow=K,ncol=M,byrow=TRUE)
  } else{

    if (K==3 && M==2){eta_atom <- matrix(c(0,0,
                                           1,1,
                                           1,0),nrow=K,ncol=M,byrow=TRUE)
    }
    # if (K==3 && M== 5){eta_atom <- matrix(c(0,0,0,0,0,
    #                                         1,1,0,0,0,
    #                                         1,1,1,1,1),nrow=K,ncol=M,byrow=TRUE)
    # }
    if (K==5 && M== 5){eta_atom <- matrix(c(0,0,0,0,0,
                                            1,1,0,1,0,
                                            1,0,0,1,0,
                                            1,1,1,0,0,
                                            0,0,1,1,1),nrow=K,ncol=M,byrow=TRUE)
    }
    #
    # eta_atom <- matrix(c(0,0,0,1,0,0,0,0,1,0,
    #                      1,1,0,1,0,1,1,0,1,0,
    #                      1,0,0,1,0,1,0,0,1,0,
    #                      1,1,1,0,0,1,1,1,0,0,
    #                      0,0,1,1,1,0,0,1,1,1),nrow=K,ncol=M,byrow=TRUE)

    if (K==2^M){eta_atom <- as.matrix(expand.grid(rep(list(0:1), M)),ncol=M)}
  }

  Z  <- sample(1:K,N,prob=rep(1/K,K),replace=TRUE)
  res_ord_tmp <- order_mat_byrow(eta_atom[Z,,drop=FALSE])
  H  <- res_ord_tmp$res
  Z  <- Z[res_ord_tmp$ord]

  xi <- (H%*%Q>0.5)+0

  res <- matrix(NA,nrow=N,ncol=L)
  for (n in 1:N){
    for (l in 1:L){
      res[n,l] <- stats::rbinom(1,1,prob=xi[n,l]*theta[l]+(1-xi[n,l])*psi[l])
    }
  }
  list(datmat = res,Q=Q,H_star= eta_atom,Z=Z,xi=xi,Eta=H)
}



#' Plot truth for simulation studies
#'
#' @param simu output from \code{\link{simulate_data}}
#' @param options_sim0 options for simulation, see \code{\link{simulate_data}}
#'
#' @return a figure with three parts, design matrix, latent profiles and
#' Q matrix
#' @export
#' @examples
#' # simulate data:
#' L0 <- 100
#' options_sim0  <- list(N = 100,  # sample size.
#'                       M = 3,    # true number of machines.
#'                       L = L0,   # number of antibody landmarks.
#'                       K = 2^3,    # number of true components.
#'                       theta = rep(0.9,L0), # true positive rates.
#'                       psi   = rep(0.1,L0), # false positive rates.
#'                       alpha1 = 1 # half of the people have the first machine.
#' )
#' YlGnBu5   <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58',"#092d94")
#' hmcols    <- colorRampPalette(YlGnBu5)(256)
#' image(simulate_data(options_sim0,SETSEED = TRUE)$datmat,col=hmcols)
#' simu     <- simulate_data(options_sim0, SETSEED=TRUE)
#' plot_truth(simu,options_sim0)
#' @importFrom grDevices colorRampPalette
plot_truth <- function(simu,options_sim0){
  YlGnBu5   <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58',"#092d94")
  hmcols    <- colorRampPalette(YlGnBu5)(256)
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
}
