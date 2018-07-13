\name{CalculateKernelCpp}
\alias{CalculateKernelCpp}
\title{
  An C++ implementation of graph kernels
}
\description{
  This function calculates a kernel matrix.
}
\usage{
CalculateKernelCpp(graph_info_list, par_r, kernel_type, store_features)
}
\arguments{
  \item{graph_info_list}{
    a list of \code{igraph} graphs
  }
  \item{par_r}{
    parameters of kernels
  }
  \item{kernel_type}{
    type of kernel
  }
  \item{store_features}{
    feature vector representations of graphs are returned if TRUE
  }
}
\value{
  \code{store_features} is \code{FALSE}
  \itemize{
    \item a kernel matrix of the respective kernel
  }

  \code{store_features} is \code{TRUE}
  \itemize{
    \item a named list with the following members:
    \describe{
      \item{\code{kernel}}{
	a kernel matrix of the respective kernel
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
graph.info.list <- vector("list", length(mutag))
for (i in 1:length(mutag)) { graph.info.list[[i]] <- GetGraphInfo(mutag[[i]]) }
K <- CalculateKernelCpp(graph.info.list, 5, 11, FALSE)
}
\keyword{ graph kernel }

