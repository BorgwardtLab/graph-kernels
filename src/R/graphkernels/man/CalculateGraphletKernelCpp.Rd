\name{CalculateGraphletKernelCpp}
\alias{CalculateGraphletKernelCpp}
\title{
  An C++ implementation of graphlet kernels
}
\description{
  This function calculates a graphlet kernel matrix.
}
\usage{
CalculateGraphletKernelCpp(graph_adj_all, graph_adjlist_all, k, connected)
}
\arguments{
  \item{graph_adj_all}{
    a list of adjacency matrices
  }
  \item{graph_adjlist_all}{
    a list of adjacency lists
  }
  \item{k}{
    the number \eqn{k}{k} of graphlet nodes
  }
  \item{connected}{
    whether or not graphlets are conneceted
  }
}
\value{
  a kernel matrix of the respective graphlet kernel
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
al.list <- as.list(rep(NA, length(mutag)))
for (i in 1:length(mutag)) { al.list[[i]] <- as_adj_list(mutag[[i]]) }
K <- CalculateGraphletKernelCpp(list(), al.list, 4, 0)
}
\keyword{ graph kernel }
\keyword{ graphlet }
