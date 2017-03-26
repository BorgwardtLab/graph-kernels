\name{CalculateVertexEdgeHistGaussKernel}
\alias{CalculateVertexEdgeHistGaussKernel}
\title{
  Gaussian RBF kernel between vertex-edge label histograms
}
\description{
  This function calculates a kernel matrix of the Gaussian RBF kernel
  \eqn{K_{VEH,G}}{K_VEH,G} between vertex-edge label histograms.
}
\usage{
CalculateVertexEdgeHistGaussKernel(G, par)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
  \item{par}{
    \eqn{\sigma}{sigma} in the Gaussian RBF kernel
  }
}
\value{
  a kernel matrix of the Gaussian RBF kernel \eqn{K_{VEH,G}}{K_VEH,G}
  between vertex-edge label histograms
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
K <- CalculateVertexEdgeHistGaussKernel(mutag, .1)
}
\keyword{ graph kernel }
\keyword{ vertex-edge label histogram }
\keyword{ Gaussian kernel }
\keyword{ RBF kernel }
