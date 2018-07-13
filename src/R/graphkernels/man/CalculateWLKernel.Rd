\name{CalculateWLKernel}
\alias{CalculateWLKernel}
\title{
  Weisfeiler-Lehman subtree kernel
}
\description{
  This function calculates a kernel matrix of the Weisfeiler-Lehman
  subtree kernel \eqn{K_{WL}}{K_WL}.
}
\usage{
CalculateWLKernel(G, par, store.features = FALSE)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
  \item{par}{
    the number \eqn{h}{h} of iterations
  }
  \item{store.features}{
    feature vector representations of graphs are returned if TRUE
  }
}
\value{
  \code{store.features} is \code{FALSE}
  \itemize{
    \item a kernel matrix of the Weisfeiler-Lehman subtree kernel \eqn{K_{WL}}{K_WL}
  }

  \code{store.features} is \code{TRUE}
  \itemize{
    \item a named list with the following members:
    \describe{
      \item{\code{kernel}}{
	a kernel matrix of the Weisfeiler-Lehman subtree kernel
	\eqn{K_{WL}}{K_WL}
      }
      \item{\code{feature}}{
	feature vector representations of graphs
      }
      \item{\code{hash}}{
	a (compressed) vertex label and neighbor (compressed) vertex labels corresponding to each feature
      }
    }
  }
}
\references{
  Shervashidze, N., Schweitzer, P., van Leeuwen, E. J., Mehlhorn, K.,
  Borgwardt, K. M.: \strong{Weisfeiler-Lehman Graph Kernels},
  \emph{Journal of Machine Learning Research}, 12, 2359-2561 (2011)
  \url{http://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf}.
}
\author{
  Mahito Sugiyama
}
\examples{
data(mutag)
K <- CalculateWLKernel(mutag, 5)
L <- CalculateWLKernel(mutag, 2, TRUE)
}
\keyword{ graph kernel }
\keyword{ Weisfeiler-Lehman subtree kernel }

