\name{mutag}
\alias{mutag}
\title{
  The mutag dataset
}
\description{
  This is the mutag dataset, a well known benchmark dataset for graph
  processing algorithms.
}
\usage{
data(mutag)
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
\keyword{ graph }
\keyword{ mutag }

