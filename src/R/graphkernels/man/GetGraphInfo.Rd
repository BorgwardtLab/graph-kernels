\name{GetGraphInfo}
\alias{GetGraphInfo}
\title{
  Necessary information of graphs for kernel computation
}
\description{
  This function extracts necessary information of graphs for kernel
  computation.
}
\usage{
GetGraphInfo(g)
}
\arguments{
  \item{g}{
    an \code{igraph} graph
  }
}
\value{
  a list of graph information with the following elements:
  \item{edge}{
    a matrix of edges with their labels
  }
  \item{vlabel}{
    a vector of vertex labels
  }
  \item{vsize}{
    the number of vertices
  }
  \item{esize}{
    the number of edges
  }
  \item{maxdegree}{
    the maximum degree
  }
}
\author{
  Mahito Sugiyama
}
\examples{
data(mutag)
ginfo <- GetGraphInfo(mutag[[1]])
}
\keyword{ graph }

