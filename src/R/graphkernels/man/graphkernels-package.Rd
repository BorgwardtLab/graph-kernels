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
  Debnath, A. K., Lopez de Compadre, R. L., Debnath, G., Shusterman,
  A. J., Hansch, C.: \strong{Structure-activity relationship of
    mutagenic aromatic and heteroaromatic nitro compounds. correlation
    with molecular orbital energies and hydrophobicity}, \emph{Journal
    of Medicinal Chemistry}, 34(2), 786-797 (1991)
  \url{http://pubs.acs.org/doi/abs/10.1021/jm00106a046}.

  Gartner, T., Flach, P., Wrobel, S.: \strong{On graph kernels:
    Hardness results and efficient alternatives}, \emph{Learning Theory
    and Kernel Machines (LNCS 2777)}, 129-143 (2003)
  \url{https://link.springer.com/chapter/10.1007/978-3-540-45167-9_11}.

  Shervashidze, N. and Schweitzer, P. and van Leeuwen, E. J. and
  Mehlhorn, K. and Borgwardt, K. M.: \strong{Weisfeiler-Lehman Graph
    Kernels}, \emph{Journal of Machine Learning Research}, 12,
  2359-2561 (2011)
  \url{http://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf}.

  Sugiyama, M., Borgwardt, K. M.: \strong{Halting in Random Walk
    Kernels}, \emph{Advances in Neural Information Processing Systems
    (NIPS 2015)}, 28, 1630-1638 (2015)
  \url{https://papers.nips.cc/paper/5688-halting-in-random-walk-kernels.pdf}.
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
