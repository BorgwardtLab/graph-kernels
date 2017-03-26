\name{CalculateExponentialRandomWalkKernel}
\alias{CalculateExponentialRandomWalkKernel}
\title{
  Exponential random walk kernel
}
\description{
  This function calculates a kernel matrix of the exponential random walk
  kernel \eqn{K_{ER}}{K_ER}.
}
\usage{
CalculateExponentialRandomWalkKernel(G, par)
}
\arguments{
  \item{G}{
    a list of \code{igraph} graphs
  }
  \item{par}{
    a coefficient \eqn{\beta}{beta}, with which the weight
    \eqn{\lambda_k}{lambda_k} for each step \eqn{k}{k} is given as \eqn{\lambda_k = \beta^k
    / k!}{lambda_k = \beta^k / k!}
  }
}
\value{
  a kernel matrix of the exponential random walk kernel \eqn{K_{ER}}{K_ER}
}
\references{
  Gartner, T., Flach, P., Wrobel, S.: \strong{On graph kernels:
    Hardness results and efficient alternatives}, \emph{Learning Theory
    and Kernel Machines (LNCS 2777)}, 129-143 (2003)
  \url{https://link.springer.com/chapter/10.1007/978-3-540-45167-9_11}.
}
\author{
  Mahito Sugiyama
}
\examples{
data(mutag)
K <- CalculateExponentialRandomWalkKernel(mutag[1:5], .1)
}
\keyword{ graph kernel }
\keyword{ random walk kernel }
\keyword{ exponential }

