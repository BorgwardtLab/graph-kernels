\name{CalculateShortestPathKernel}
\alias{CalculateShortestPathKernel}
\title{
  Shortest-path kernel
}
\description{
  This function calculates a kernel matrix of the shortest-path kernel
  \eqn{K_{SP}}{K_SP}.
}
\usage{
CalculateShortestPathKernel(G)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
}
\value{
  a kernel matrix of the shortest-path kernel \eqn{K_{SP}}{K_SP}
}
\references{
  Borgwardt, K. M., Kriegel, H.-P.: \strong{Shortest-Path Kernels on
    Graphs}, \emph{Proceedings of the 5th IEEE International Conference on
    Data Mining (ICDM'05)}, 74-81 (2005)
  \url{http://ieeexplore.ieee.org/document/1565664/}.
}
\author{
  Mahito Sugiyama
}
\examples{
data(mutag)
K <- CalculateShortestPathKernel(mutag)
}
\keyword{ graph kernel }
\keyword{ shortest path }
