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
  Gartner, T., Flach, P., Wrobel, S.: \strong{On graph kernels:
    Hardness results and efficient alternatives}, \emph{Learning Theory
    and Kernel Machines (LNCS 2777)}, 129-143 (2003)
  \url{https://link.springer.com/chapter/10.1007/978-3-540-45167-9_11}.

  Sugiyama, M., Borgwardt, K. M.: \strong{Halting in Random Walk
    Kernels}, \emph{Advances in Neural Information Processing Systems
    (NIPS 2015)}, 28, 1630-1638 (2015)
  \url{https://papers.nips.cc/paper/5688-halting-in-random-walk-kernels.pdf}.
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

