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
}
\keyword{ graph kernel }
\keyword{ Weisfeiler-Lehman subtree kernel }

