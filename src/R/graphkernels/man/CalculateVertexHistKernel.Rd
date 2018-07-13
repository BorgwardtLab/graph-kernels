\name{CalculateVertexHistKernel}
\alias{CalculateVertexHistKernel}
\title{
  Linear kernel between vertex label histograms
}
\description{
  This function calculates a kernel matrix of the linear kernel
  \eqn{K_{VH}}{K_VH} between vertex label histograms.
}
\usage{
CalculateVertexHistKernel(G, store.features = FALSE)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
  \item{store.features}{
    feature vector representations of graphs are returned if TRUE
  }
}
\value{
  \code{store.features} is \code{FALSE}
  \itemize{
    \item a kernel matrix of the linear kernel \eqn{K_{VH}}{K_VH}
    between vertex label histograms
  }

  \code{store.features} is \code{TRUE}
  \itemize{
    \item a named list with the following members:
    \describe{
      \item{\code{kernel}}{
	a kernel matrix of the linear kernel \eqn{K_{VH}}{K_VH}
	between vertex label histograms
      }
      \item{\code{feature}}{
	feature vector representations of graphs
      }
    }
  }
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
K <- CalculateVertexHistKernel(mutag)
}
\keyword{ graph kernel }
\keyword{ vertex label histogram }
