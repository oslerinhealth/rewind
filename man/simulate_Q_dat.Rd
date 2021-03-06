% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_mcmc.R
\name{simulate_Q_dat}
\alias{simulate_Q_dat}
\title{Randomly generate a Q matrix within the identifying condition (data-driven)}
\usage{
simulate_Q_dat(M, dat, p = 0.1, frac = 1/4)
}
\arguments{
\item{M}{latent state dimension}

\item{dat}{binary data matrix (rows for observations, columns for dimensions)}

\item{p}{Bernoulli probability of 1 in the Q matrix (except two diagonal matrices)}

\item{frac}{A threshold - this function only initializes the dimensions with at least
\code{frac}*100 percent observed frequencies.}
}
\value{
a binary matrix of dimension M by L
}
\description{
This function initializes Q (if unknown) during MCMC sampling chain within identifiability
constraints. It is a warm start - because it will not assign an "1" to dimension "l" with
few ones in the data. NB: harder to get 1 to zero? easy to get zero to one?
}
\examples{
# simulate data:
L0 <- 100
options_sim0  <- list(N = 200,  # sample size.
                     M = 3,   # true number of machines.
                     L = L0,   # number of antibody landmarks.
                     K = 8,    # number of true components.,
                     theta = rep(0.8,L0), # true positive rates
                     psi   = rep(0.01,L0), # false positive rates
                     alpha1 = 1 # half of the people have the first machine.
)

#image(simulate_data(options_sim0,SETSEED = TRUE)$datmat)
simu     <- simulate_data(options_sim0, SETSEED=TRUE)
simu_dat <- simu$datmat
simulate_Q_dat(5,simu_dat)

}
