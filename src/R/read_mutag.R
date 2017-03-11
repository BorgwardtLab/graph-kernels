library(igraph)
read.mutag <- function() {
  lb <- unlist(read.table("../../data/mutag.label"))
  names(lb) <- NULL
  n <- length(lb)
  G <- vector("list", n)
  for (i in 1:n) {
    G[[i]] <- read_graph(paste0("../../data/mutag/mutag_", i, ".graphml"), format = "graphml")
    graph_attr(G[[i]], "label") <- lb[i]
  }
  G
}
