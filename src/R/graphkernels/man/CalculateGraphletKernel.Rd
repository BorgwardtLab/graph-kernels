\name{CalculateGraphletKernel}
\alias{CalculateGraphletKernel}
\title{
  Graphlet kernel
}
\description{
  This function calculates a kernel matrix of the graphlet kernel
  \eqn{K_{GL}}{K_GL} between unlabeled graphs.
}
\usage{
CalculateGraphletKernel(G, par, store.features = FALSE)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
  \item{par}{
    the number \eqn{k}{k} of graphlet nodes (\eqn{k = 3}{k = 3} or 4 is supported)
  }
  \item{store.features}{
    feature vector representations of graphs are returned if TRUE
  }
}
\value{
  \code{store.features} is \code{FALSE}
  \itemize{
    \item a kernel matrix of the graphlet kernel \eqn{K_{GL}}{K_GL}
  }

  \code{store.features} is \code{TRUE}
  \itemize{
    \item a named list with the following members:
    \describe{
      \item{\code{kernel}}{
	a kernel matrix of the graphlet kernel \eqn{K_{GL}}{K_GL}
      }
      \item{\code{feature}}{
	feature vector representations of graphs
      }
    }
  }
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
K <- CalculateGraphletKernel(mutag, 4)
}
\keyword{ graph kernel }
\keyword{ graphlet }
