\name{graphkernels-package}
\alias{graphkernels-package}
\alias{graphkernels}
\docType{package}
\title{
  Graph Kernels
}
\description{
  A fast C++ implementation for computing various graph kernels
  including (1) simple kernels between vertex and/or edge label
  histograms, (2) graphlet kernels, (3) random walk kernels (popular
  baselines), and (4) the Weisfeiler-Lehman graph kernel
  (state-of-the-art).
  
}
\details{
  This library provides the following graph kernels:
  \itemize{
    \item the linear kernel between vertex label histograms
    \item the linear kernel between edge label histograms
    \item the linear kernel between vertex-edge label histograms
    \item the linear kernel combination vertex label histograms and
    vertex-edge label histograms
    \item the Gaussian RBF kernel between vertex label histograms
    \item the Gaussian RBF kernel between edge label histograms
    \item the Gaussian RBF kernel between vertex-edge label histograms
    \item the graphlet kernel
    \item the \eqn{k}{k}-step random walk kernel
    \item the geometric random walk kernel
    \item the exponential random walk kernel
    \item the shortest-path kernel
    \item the Weisfeiler-Lehman subtree kernel
  }
  Given a list of \pkg{igraph} graphs, each function calculates the
  corresponding kernel (Gram) matrix.
}
\author{
Mahito Sugiyama

Maintainer: Mahito Sugiyama \code{<mahito@nii.ac.jp>}
}
\references{
  Borgwardt, K. M., Kriegel, H.-P.: \strong{Shortest-Path Kernels on
    Graphs}, \emph{Proceedings of the 5th IEEE International Conference on
    Data Mining (ICDM'05)}, 74-81 (2005)
  \url{http://ieeexplore.ieee.org/document/1565664/}.

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

  Shervashidze, N., Schweitzer, P., van Leeuwen, E. J., Mehlhorn, K.,
  Borgwardt, K. M.: \strong{Weisfeiler-Lehman Graph Kernels},
  \emph{Journal of Machine Learning Research}, 12, 2359-2561 (2011)
  \url{http://www.jmlr.org/papers/volume12/shervashidze11a/shervashidze11a.pdf}.

  Shervashidze, N., Vishwanathan, S. V. N., Petri, T., Mehlhorn, K.,
  Borgwardt, K. M.: \strong{Efficient Graphlet Kernels for Large Graph
    Comparison}, \emph{Proceedings of the 12th International Conference
    on Artificial Intelligence and Statistics (AISTATS)}, 5, 488-495
  (2009)
  \url{http://proceedings.mlr.press/v5/shervashidze09a.html}.

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
KWL <- CalculateWLKernel(mutag, 5)
  ## compute Weisfeiler-Lehman subtree kernel
}
