\name{CalculateEdgeHistKernel}
\alias{CalculateEdgeHistKernel}
\title{
  Linear kernel between edge label histograms
}
\description{
  This function calculates a kernel matrix of the linear kernel
  \eqn{K_{EH}}{K_EH} between edge label histograms.
}
\usage{
CalculateEdgeHistKernel(G)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
}
\value{
  a kernel matrix of the linear kernel \eqn{K_{EH}}{K_EH} between edge
  label histograms
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
K <- CalculateEdgeHistKernel(mutag)
}
\keyword{ graph kernel }
\keyword{ edge label histogram }
