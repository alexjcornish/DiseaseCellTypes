gsc.sub <- function(
    expression,
    genes,
    g,
    rwr.r, 
    rwr.cutoff
) {
    # run the compactness function on a cellular network created using expression and genes
    # this function will be run a large number of times and therefore is important to reduce the time required
    # because of this, no checks are run within the function
    # all checks should be run within the compactness.sig function
    # a single compactness score is output
    
    edge.attr <- "edge.attr"
    genes.vnum <- get.vertex.numbers(g, genes)
    genes.vnum.size <- length(genes.vnum)
    
    # score edges in g
    el <- get.edgelist(g, names=TRUE)
    expression[expression == 0] <- min(expression[expression != 0], na.rm=T)
    edge.weights <- (expression[el[, 1]] * expression[el[, 2]])
    
    # extract column-normalised adjacency matrix
    g <- set.edge.attribute(g, edge.attr, value=edge.weights)
    W <- get.adjacency.norm(g, edge.attr)
    
    # p: a matrix of random walk probabilities
    # plast: p in the previous time point
    # p0: p at t=0
    plast <- p <- p0 <- as(sparseMatrix(dims=c(vcount(g), genes.vnum.size), i=genes.vnum, j=1:genes.vnum.size, x=rep(1, genes.vnum.size)), "dgeMatrix")
    
    # run first RWR iteration
    p <- (1 - rwr.r) * W %*% plast
    p[cbind(genes.vnum, 1:genes.vnum.size)] <- p[cbind(genes.vnum, 1:genes.vnum.size)] + rwr.r
    
    # run remaining RWR iterations
    rwr.cutoff.total <- rwr.cutoff * genes.vnum.size
    while (sum(abs(p - plast)) > rwr.cutoff.total) {
        plast <- p
        p <- (1 - rwr.r) * W %*% plast
        p[cbind(genes.vnum, 1:genes.vnum.size)] <- p[cbind(genes.vnum, 1:genes.vnum.size)] + rwr.r
    }
    
    # rank distances
    d <- apply(-p, 2, rank)[genes.vnum, ]
    
    # compute average distance and output
    mean(d)
}
