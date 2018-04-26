\name{CalculateConnectedGraphletKernel}
\alias{CalculateConnectedGraphletKernel}
\title{
  Connected graphlet kernel
}
\description{
  This function calculates a kernel matrix of the graphlet kernel with
  connected graphlets \eqn{K_{CGL}}{K_CGL} between unlabeled graphs.
}
\usage{
CalculateConnectedGraphletKernel(G, par)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
  \item{par}{
    the number \eqn{k}{k} of graphlet nodes (\eqn{k = 3}{k = 3}, 4, or 5 is supported)
  }
}
\value{
  a kernel matrix of the connected graphlet kernel \eqn{K_{CGL}}{K_CGL}
}
\references{
  Shervashidze, N., Vishwanathan, S. V. N., Petri, T., Mehlhorn, K.,
  Borgwardt, K. M.: \strong{Efficient Graphlet Kernels for Large Graph
    Comparison}, \emph{Proceedings of the 12th International Conference
    on Artificial Intelligence and Statistics (AISTATS)}, 5, 488-495
  (2009)
  \url{http://proceedings.mlr.press/v5/shervashidze09a.html}.
}
\author{
  Mahito Sugiyama
}
\examples{
data(mutag)
K <- CalculateConnectedGraphletKernel(mutag, 4)
}
\keyword{ graph kernel }
\keyword{ graphlet }
