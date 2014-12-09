get.adjacency.norm <- function(
    g,
    edge.attr
) {
    # compute a column-normalised adjacency matrix to be used to compute the RWR distance
    # function outputs a row for each vertex in g and a column for each vertex in g
    # a sparse matrix is always output
    # g should be an undirected igraph object and edge.attr a edge attribute containing weights
    # if edge.attr is NULL, then it is assumed that all edge weights equal 1
    # in order to speed up the function, numbers rather than names are used to refer to vertices
    # this function is designed to be as fast as possible, and therefore inputs are not checked
    weights <- if (is.null(edge.attr)) 1 else get.edge.attribute(g, edge.attr)
    el <- cbind(get.edgelist(g, names=F), weights) # el: edgelist with weights
    el <- rbind(el, el[el[, 1] != el[, 2], c(2,1,3)])
    weights.colsum <- sapply(split(el[, 3], el[, 1]), sum)
    el[, 3] <- el[, 3] / weights.colsum[el[, 1]]
    sparseMatrix(dims=rep(vcount(g), 2), i=el[, 2], j=el[, 1], x=el[, 3], symmetric=FALSE)
}
