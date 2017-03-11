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
CalculateWLKernel(G, par)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
  \item{par}{
    the number \eqn{h}{h} of iterations
  }
}
\value{
  a kernel matrix of the Weisfeiler-Lehman subtree kernel \eqn{K_{WL}}{K_WL}
}
\references{
  Sugiyama, M., Borgwardt, K. M.: \strong{Halting in Random Walk Kernels},
  \emph{Advances in Neural Information Processing Systems (NIPS 2015)}, 2015.
}
\author{
  Mahito Sugiyama
}
\examples{
data(mutag)
K <- CalculateWLKernel(mutag, 5)
}
\keyword{ graph kernel }
\keyword{ Weisfeiler-Lehman subtree kernel }

