## library(graphkernels); data(mutag); res <- CalculateVertexEdgeHistKernel(mutag, TRUE)
## data(mutag)
## res <- CalculateEdgeHistKernel(mutag)
## head(res[[2]])

## res <- CalculateWLKernel(mutag, 0)

## K <- res[[2]] %*% t(res[[2]])
## sum(abs(res[[1]] - K))
