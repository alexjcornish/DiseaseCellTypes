add.loops.to.lonely <- function(
    g, 
    edge.attr=NULL, 
    loop.weight=1
) {
    # add a self loop to each lone vertex in g
    # this function is used in computing the RWR distance
    # if an edge attribute is input, the weight of edge of these loops is set to loop.weight
    
    # check input
    if (class(g) != "igraph") stop("g is not an igraph object")
    if (!is.numeric(loop.weight)) stop("loop weight is not numeric")
    if (!is.null(edge.attr) & ecount(g) > 0) {
        # ecount(g) > 0 requirement present because if no edges are present, no edge weights can be present
        if (!edge.attr %in% list.edge.attributes(g)) stop("edge attr not found on g")
        if (!is.numeric(get.edge.attribute(g, edge.attr))) stop("edge.attr not numeric")
    }
    
    # if there are any lonely nodes, add self loops
    lonely.nodes <- which(degree(g) == 0)
    if (length(lonely.nodes) > 0) {
        g <- add.edges(g, rep(lonely.nodes, each=2))
        
        # add edge weights if required
        if (!is.null(edge.attr)) {
            # if an input g has no edges, then it is neccessary to add the edge weights
            weights <- if (ecount(g) == length(lonely.nodes)) rep(NA, ecount(g)) else get.edge.attribute(g, edge.attr) 
            weights[is.na(weights)] <- loop.weight
            g <- set.edge.attribute(g, edge.attr, value=weights)
        }
    }

    g
}
