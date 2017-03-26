\name{CalculateGeometricRandomWalkKernel}
\alias{CalculateGeometricRandomWalkKernel}
\title{
  Geometric random walk kernel
}
\description{
  This function calculates a kernel matrix of the geometric random walk
  kernel \eqn{K_{GR}}{K_GR}.
}
\usage{
CalculateGeometricRandomWalkKernel(G, par)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
  \item{par}{
    a coefficient \eqn{\lambda}{lambda}, with which the weight
    \eqn{\lambda_k}{lambda_k} for each step \eqn{k}{k} is given as
    \eqn{\lambda_k = \lambda^k}{lambda_k = lambda^k}
  }
}
\value{
  a kernel matrix of the geometric random walk kernel \eqn{K_{GR}}{K_GR}
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
K <- CalculateGeometricRandomWalkKernel(mutag, .1)
}
\keyword{ graph kernel }
\keyword{ random walk kernel }
\keyword{ geometric }

