CalculateEdgeHistKernel <- function(G, store.features = FALSE) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  res <- CalculateKernelCpp(graph.info.list, 0, 1, store.features)[1:2]
  if (!store.features) res <- res$kernel
  res
}

CalculateVertexHistKernel <- function(G, store.features = FALSE) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  res <- CalculateKernelCpp(graph.info.list, 0, 2, store.features)[1:2]
  if (!store.features) res <- res$kernel
  res
}

CalculateVertexEdgeHistKernel <- function(G, store.features = FALSE) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  res <- CalculateKernelCpp(graph.info.list, 0, 3, store.features)[1:2]
  if (!store.features) res <- res$kernel
  res
}

CalculateVertexVertexEdgeHistKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 4, FALSE)$kernel
}

CalculateEdgeHistGaussKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 5, FALSE)$kernel
}

CalculateVertexHistGaussKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 6, FALSE)$kernel
}

CalculateVertexEdgeHistGaussKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 7, FALSE)$kernel
}

CalculateGeometricRandomWalkKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 8, FALSE)$kernel
}

CalculateExponentialRandomWalkKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 9, FALSE)$kernel
}

CalculateKStepRandomWalkKernel <- function(G, par) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  CalculateKernelCpp(graph.info.list, par, 10, FALSE)$kernel
}

CalculateWLKernel <- function(G, par, store.features = FALSE) {
  graph.info.list <- vector("list", length(G))
  for (i in 1:length(G))
    graph.info.list[[i]] <- GetGraphInfo(G[[i]])
  res <- CalculateKernelCpp(graph.info.list, par, 11, store.features)
  if (!store.features) res <- res$kernel
  res
}

CalculateConnectedGraphletKernel <- function(G, par, store.features = FALSE) {
  if (par < 3) {
    par <- 3
    warning("k = 3 is used (k = 3, 4, or 5 is supported).")
  }
  if (par > 5) {
    par <- 5
    warning("k = 5 is used (k = 3, 4, or 5 is supported).")
  }
  al.list <- as.list(rep(NA, length(G)))
  for (i in 1:length(G)) {
    al.list[[i]] <- as_adj_list(G[[i]])
  }
  am.list <- as.list(rep(NA, length(G)))
  for (i in 1:length(G)) {
    am.list[[i]] <- as_adj(G[[i]])
  }
  res <- CalculateGraphletKernelCpp(am.list, al.list, par, 1, store.features)
  if (!store.features) res <- res$kernel
  res
}

CalculateGraphletKernel <- function(G, par, store.features = FALSE) {
  if (par < 3) {
    par <- 3
    warning("k = 3 is used (k = 3 or 4 is supported).")
  }
  if (par > 4) {
    par <- 4
    warning("k = 4 is used (k = 3 or 4 is supported).")
  }
  al.list <- as.list(rep(NA, length(G)))
  for (i in 1:length(G)) {
    al.list[[i]] <- as_adj_list(G[[i]])
  }
  res <- CalculateGraphletKernelCpp(list(), al.list, par, 0, store.features)
  if (!store.features) res <- res$kernel
  res
}

CalculateShortestPathKernel <- function(G) {
  G.floyd <- as.list(rep(NA, length(G)))
  for (i in 1:length(G)) {
    D <- distances(G[[i]])
    G.floyd[[i]] <- make_full_graph(vcount(G[[i]])) %>% set_edge_attr("weight", value = D[lower.tri(D)])
  }
  CalculateKStepRandomWalkKernel(G, c(0, 1))
}

GetGraphInfo <- function(g) {
  ## a vertex attribute is missing
  if (length(vertex_attr(g)) == 0)
      g <- set_vertex_attr(g, "label", value = rep(1, vcount(g)))
  ## there are multiple vertex attributes
  if (length(vertex_attr(g)) > 1) {
    warning(paste0("There are multiple vertex attributes! The first attribute \"",
                   names(vertex_attr(g))[1],  "\" is used."))
  }
  ## change name of labels to "label"
  names(vertex_attr(g))[1] <- "label"
  ## remove non-integer labels
  l <- as.integer(vertex_attr(g)$label)
  if (sum(is.na(l)) > 0) {
    warning("Non-integer labels are converted to 0.")
    l[is.na(l)] <- 0
  }
  g <- g %>% set_vertex_attr("label", value = l)

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
  if (length(E(g)) > 0)
    E <- cbind(E, edge_attr(g)[[1]])

  ## a vector of a vertex attribute
  v.label <- as.integer(vertex_attr(g)[[1]])
  res <- list(edge = E, vlabel = v.label,
              vsize = vcount(g), esize = ecount(g), maxdegree = max(degree(g)))
  res
}
