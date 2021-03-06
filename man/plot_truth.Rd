% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{plot_truth}
\alias{plot_truth}
\title{Plot truth for simulation studies}
\usage{
plot_truth(simu, options_sim0)
}
\arguments{
\item{simu}{output from \code{\link{simulate_data}}}

\item{options_sim0}{options for simulation, see \code{\link{simulate_data}}}
}
\value{
a figure with three parts, design matrix, latent profiles and
Q matrix
}
\description{
Plot truth for simulation studies
}
\examples{
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
YlGnBu5   <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58',"#092d94")
hmcols    <- colorRampPalette(YlGnBu5)(256)
image(simulate_data(options_sim0,SETSEED = TRUE)$datmat,col=hmcols)
simu     <- simulate_data(options_sim0, SETSEED=TRUE)
plot_truth(simu,options_sim0)
}
