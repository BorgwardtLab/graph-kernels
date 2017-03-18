\name{graphkernels-package}
\alias{graphkernels-package}
\alias{graphkernels}
\docType{package}
\title{
  Graph Kernels
}
\description{
  A fast C++ implementation of various graph kernels.
}
\details{
  This library provides a fast C++ implementation of various graph kernels including:
  \itemize{
    \item simple kernels between vertex and/or edge label histograms,
    \item random walk kernels (popular baselines), and
    \item the Weisfeiler-Lehman graph kernel (state-of-the-art).
  }
}
\author{
Mahito Sugiyama

Maintainer: Mahito Sugiyama \code{<mahito@ar.sanken.osaka-u.ac.jp>}
}
\references{
  Sugiyama, M., Borgwardt, K. M.: \strong{Halting in Random Walk Kernels},
  \emph{Advances in Neural Information Processing Systems (NIPS 2015)}, 2015.
}
\keyword{ graph kernel }
\examples{
data(mutag)
KEH <- CalculateEdgeHistKernel(mutag)
  ## compute linear kernel between edge histograms
KkR <- CalculateKStepRandomWalkKernel(mutag, c(.1, .1))
  ## compute k-step random walk kernel (k = 2 in this case)
KWL <- CalculateWLKernel(mutag, 5)
  ## compute Weisfeiler-Lehman subtree kernel
}
