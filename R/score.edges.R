score.edges <- function(
    g,
    expression, 
    edge.attr.type=c("weight", "distance"), # whether the edge attribute to add represents weights or distances
    edge.attr="score", # the name of the edge attribute under which distances/weights should be added
    vertex.attr="expression", # the name of the vertex attribute under which expression scores should be added
    expression.miss=median(expression),
    correct.inf=TRUE # if T, then infinite values are not returned
) {
    # reweight the edges of a network using expression data
    # expression should be a named numerical vector of transformed expression values
    # the function doesn't currently take old edge weights into account, but instead assumes that all weights are 1
  
    # setup
    edge.attr.type <- match.arg(edge.attr.type)
    if (!is.igraph(g)) stop("g is not an igraph object")
    if (!is.vector(expression)) stop("expression is not a numeric vector")
    genes.network <- V(g)$name
    genes.expression <- names(expression)
    if (is.null(genes.network)) stop("g vertices should be named")
    if (is.null(genes.expression)) stop("expression should be named")
    if (sum(genes.expression %in% genes.network) == 0) warning("no overlap between g node and expression names")
    
    # add a genes present in the network but missing from the expression to the expression
    genes.missing <- genes.network[!genes.network %in% genes.expression]
    expression <- c(expression, structure(rep(expression.miss, length(genes.missing)), names=genes.missing))
    
    # compute weights (convert to distances if required) and add to network
    el <- get.edgelist(g, names=T)
    if (correct.inf) expression[expression == 0] <- min(expression[expression != 0], na.rm=T)
    weights <- (expression[el[, 1]] * expression[el[, 2]])
    distances <- 1 / weights
    g <- switch(edge.attr.type,
        weight = set.edge.attribute(g, edge.attr, value=as.double(weights)),
        distance = set.edge.attribute(g, edge.attr, value=as.double(distances)) 
    )
    
    # add expression scores to vertices 
    vertex.scores <- expression[V(g)$name]
    g <- set.vertex.attribute(g, name=vertex.attr, value=as.double(vertex.scores))
    
    g
}

