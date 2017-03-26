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
  Debnath, A. K., Lopez de Compadre, R. L., Debnath, G., Shusterman,
  A. J., Hansch, C.: \strong{Structure-activity relationship of
    mutagenic aromatic and heteroaromatic nitro compounds. correlation
    with molecular orbital energies and hydrophobicity}, \emph{Journal
    of Medicinal Chemistry}, 34(2), 786-797 (1991)
  \url{http://pubs.acs.org/doi/abs/10.1021/jm00106a046}.
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

