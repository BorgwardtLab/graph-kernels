\name{CalculateExponentialRandomWalkKernel}
\alias{CalculateExponentialRandomWalkKernel}
\title{
  Exponential random walk kernel
}
\description{
  This function calculates a kernel matrix of the exponential random walk
  kernel \eqn{K_{ER}}{K_ER}.
}
\usage{
CalculateExponentialRandomWalkKernel(G, par)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
  \item{par}{
    a coefficient \eqn{\beta}{beta}, with which the weight
    \eqn{\lambda_k}{lambda_k} for each step \eqn{k}{k} is given as \eqn{\lambda_k = \beta^k
    / k!}{lambda_k = \beta^k / k!}
  }
}
\value{
  a kernel matrix of the exponential random walk kernel \eqn{K_{ER}}{K_ER}
}
\references{
  Sugiyama, M., Borgwardt, K. M.: \strong{Halting in Random Walk Kernels},
  \emph{Advances in Neural Information Processing Systems (NIPS 2015)}, 2015.
}
\author{
  Mahito Sugiyama
}
\examples{
data(mutag)
K <- CalculateExponentialRandomWalkKernel(mutag[1:5], .1)
}
\keyword{ graph kernel }
\keyword{ random walk kernel }
\keyword{ exponential }

