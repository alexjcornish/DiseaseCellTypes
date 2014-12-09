distance.rwr <- function(
    g,
    v=V(g),
    edge.attr=NULL,
    rwr.r=0.7, 
    rwr.cutoff=1e-5,
    correct.inf=TRUE
) {
    # compute a RWR-based distance matrix for g using an iterative approach
    # implementation is based on the random walk with restart algorithm from Kohler 2008
    # originaly the function was used to measure the distance between a set and each node, however it can be modified to measure the distance between all node pairs
    # this function contains 
    
    # r is the restart probablility, Kohler2008 don't give details about the value they use or test
    # in Kohler2008, they say values of r between 0.5 and 0.9 are optimal
    # in Glaab2012, they use a value of r = 0.9
    # edge attributes should represent weights, not distances
    # the interations stop when the mean absolute difference for each node falls below rwr.cutoff * v.n.nodes
    # if a node is unconnected, the random walker does not leave this node
    # currently, I convert probabilities to distances using the method from Glaab2012, while additionally setting the diagonal equal to 0
    # currently, this algorithm does not neccessarily produce a symetric matrix

    # check input
    if (!is.igraph(g)) stop("g is not an igraph object")
    if (is.directed(g)) stop("g is directed")
    if (!is.numeric(v) & !is.character(v)) stop("v is not an igraph object, numeric or character vector")
    if (!is.null(edge.attr) & ecount(g) > 0) {
        weights <- get.edge.attribute(g, edge.attr)
        if (is.null(weights)) stop("edge.attr not found on g")
        if (!is.numeric(weights)) stop("edge weights not numeric")
        if (!all(is.finite(weights))) stop("infinite edge weights present")
        if (!all(!is.na(weights))) stop("edge weights cannot be equal to NA")
        if (!all(!is.null(weights))) stop("edge weights cannot be NULL")
        if (!all(weights >= 0)) stop("edge weights should be >= 0")
    }
    if (rwr.r < 0 | rwr.r > 1) stop("rwr.r should be between 0 and 1")
    if (rwr.cutoff <= 0) stop("iteration cutoff should be > 0")
    g.n.nodes <- vcount(g)
    if (!g.n.nodes) return(array(0, dim=c(0,0), dimnames=list(V(g)$name, V(g)$name))) # if no vertices in g, return empty matrix
    
    # setup
    g <- add.loops.to.lonely(g, edge.attr) # add self-loops to unconnected nodes to prevent error
    v.num <- get.vertex.numbers(g, v) # convert input vertices to vertex numbers
    v.n.nodes <- length(v.num)
    W <- get.adjacency.norm(g,  edge.attr) # get W, a column-normalised adjacency matrix

    # p: a matrix of random walk probabilities
    # plast: p in the previous time point
    # p0: p at t=0
    plast <- p <- p0 <- as(sparseMatrix(dims=c(g.n.nodes, v.n.nodes), i=v.num, j=1:v.n.nodes, x=rep(1, v.n.nodes)), "dgeMatrix")
    
    # run first iteration
    p <- (1 - rwr.r) * W %*% plast
    p[cbind(v.num, 1:v.n.nodes)] <- p[cbind(v.num, 1:v.n.nodes)] + rwr.r # add the probability of returning to the start
    
    # run remaining iterations
    ic.total <- rwr.cutoff * v.n.nodes
    while (sum(abs(p - plast)) > ic.total) {
        plast <- p
        p <- (1 - rwr.r) * W %*% plast
        p[cbind(v.num, 1:v.n.nodes)] <- p[cbind(v.num, 1:v.n.nodes)] + rwr.r
    }
    
    # convert probability matrix p to a distance matrix d
    # transposing p ensures that p[i,j] is probability distribution for a walker leaving i
    # Glaab2012 convert probabilities by subtracting the distance from 1 (p <- 1 - p)
    # however, this can produce problems if there are probabilities smaller than 1e-16
    # theres, I will use 1 / p, which can be converted to 1 - p if required
    d <- as.matrix(1 / t(p))
    d[cbind(1:v.n.nodes, v.num)] <- 0
    
    # correct infinite values within the distance matrix 
    # not required when using the Glaab2012 method of converting probabilities to distances
    if (correct.inf) d[!is.finite(d)] <- max(d[is.finite(d)], na.rm=T)
    
    # output distance matrix (the modified p)
    dimnames(d) <- list(V(g)$name[v.num], V(g)$name)
    d
}
