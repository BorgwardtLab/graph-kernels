CalculateEdgeHistKernel <- function(G) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, 0, 1)
}

CalculateVertexHistKernel <- function(G) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, 0, 2)
}

CalculateVertexEdgeHistKernel <- function(G) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, 0, 3)
}

CalculateVertexVertexEdgeHistKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 4)
}

CalculateEdgeHistGaussKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 5)
}

CalculateVertexHistGaussKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 6)
}

CalculateVertexEdgeHistGaussKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 7)
}

CalculateGeometricRandomWalkKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 8)
}

CalculateExponentialRandomWalkKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 9)
}

CalculateKStepRandomWalkKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 10)
}

CalculateWLKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 11)
}

GetGraphInfo <- function(g) {
  ## an edge matrix
  E <- as_edgelist(g)
  E <- matrix(as.integer(E), ncol = ncol(E))
  ## there are multiple edge attributes
  if (length(edge_attr(g)) > 1) {
    warning(paste0("There are multiple edge attributes! The first attribute \"",
                   names(edge_attr(g))[1],  "\" is used."))
  }
  ## an edge attribute is missing
  if (length(edge_attr(g)) == 0)
    g <- set_edge_attr(g, "label", value = rep(1, ecount(g)))
  E <- cbind(E, edge_attr(g)[[1]])

  ## a vector of a vertex attribute
  ## there are multiple vertex attributes
  if (length(vertex_attr(g)) > 1) {
    warning(paste0("There are multiple vertex attributes! The first attribute \"",
                   names(vertex_attr(g))[1],  "\" is used."))
  }
  if (length(vertex_attr(g)) == 0)
    g <- set_vertex_attr(g, "label", value = rep(1, vcount(g)))
  v.label <- as.integer(vertex_attr(g)[[1]])
  res <- list(edge = E, vlabel = v.label,
              vsize = vcount(g), esize = ecount(g), maxdegree = max(degree(g)))
  res
}
