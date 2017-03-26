\name{CalculateVertexVertexEdgeHistKernel}
\alias{CalculateVertexVertexEdgeHistKernel}
\title{
  Linear kernel combination of vertex label histograms and vertex-edge
  label histograms
}
\description{
  This function calculates a kernel matrix of the linear kernel
  combination \eqn{K_{H}}{K_H} of vertex label histograms
  \eqn{K_{VH}}{K_VH} and vertex-edge label histograms \eqn{K_{VEH}}{K_VEH}.
}
\usage{
CalculateVertexVertexEdgeHistKernel(G, par)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
  \item{par}{
    a coefficient \eqn{\lambda}{lambda}, with which the resulting kernel
    is given as \eqn{K_{VH} + \lambda K_{VEH}}{K_VH + lambda K_VEH}
  }
}
\value{
  a kernel matrix that is equivalent to \eqn{K_{VH} + \lambda K_{VEH}}{K_VH + lambda K_VEH}
}
\references{
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
K <- CalculateVertexVertexEdgeHistKernel(mutag, .1)
}
\keyword{ graph kernel }
\keyword{ vertex label histogram }
\keyword{ vertex-edge label histogram }
\keyword{ linear combination }

