\name{CalculateKStepRandomWalkKernel}
\alias{CalculateKStepRandomWalkKernel}
\title{
  k-step random walk kernel
}
\description{
  This function calculates a kernel matrix of the \eqn{k}{k}-step random
  walk kernel \eqn{K_{\times}^{k}}{K_x^k}.
}
\usage{
CalculateKStepRandomWalkKernel(G, par)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
  \item{par}{
    a vector of coefficients \eqn{\lambda_0, \lambda_1, \dots,
      \lambda_k}{lambda_0, lambda_1, ..., lambda_k}
  }
}
\value{
  a kernel matrix of the k-step random walk kernel \eqn{K_{\times}^{k}}{K_x^k}
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
K <- CalculateKStepRandomWalkKernel(mutag, rep(1, 2))
}
\keyword{ graph kernel }
\keyword{ k-step walk kernel }

