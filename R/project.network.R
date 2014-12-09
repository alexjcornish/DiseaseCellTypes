project.network <- function(
    p.values,
    col.vert=NULL,
    col.other="white",
    method.cor=c("pearson", "kendall", "spearman"), 
    rank.max=4,
    vert.size.max=10,
    vert.label.cex=1,
    vert.p.cutoff=0.1,
    edge.width.max=10,
    edge.col.diff="lightgrey",
    ... # parameters passed to the plot.igraph function
) {
    # project a network using the correlations between the -log10 p-values
    
    # setup
    method.cor <- match.arg(method.cor)
    if (!is.null(col.vert)) {
        if (!all(rownames(p.values) %in% names(col.vert))) stop("not all p.values rownames are in col.vert")
    } 
    
    # compute correlations between -log10 p-values
    adj <- p.values.ml10.corr <- cor(t(-log10(p.values)), method=method.cor)
    
    # set some cells to zero to create adjacency matrix
    diag(adj) <- 0
    adj <- apply(adj, 1, function(row) ifelse(row %in% sort(row, decreasing=T)[1:rank.max], row, 0))
    dimnames(adj) <- dimnames(p.values.ml10.corr)
    g <- graph.adjacency(adj, mode="max", weighted="width", diag=FALSE, add.rownames=TRUE)
    
    # add vertex sizes, which depend upon the number of associations
    v.size <- sqrt(apply(p.values <= vert.p.cutoff, 1, sum) / pi)
    v.size <- v.size * vert.size.max / max(v.size)
    V(g)$size <- v.size
    
    # add vertex colours
    V(g)$color <- if (is.null(col.vert)) rep(col.other, vcount(g)) else col.vert[V(g)$name]
    V(g)$color[is.na(V(g)$color)] <- col.other
    V(g)$color[V(g)$color == "other"] <- col.other
    V(g)$color[V(g)$color == "Other"] <- col.other
    V(g)$frame.color <- V(g)$color # don't include borders
    
    # format vertex labels
    V(g)$label.cex <- rep(vert.label.cex, vcount(g))
    
    # add edge wights
    widths <- E(g)$width
    widths <- widths * edge.width.max / max(widths)
    widths[widths <= 0] <- min(widths[widths > 0])
    E(g)$width <- widths
    
    # add edge colours
    el <- get.edgelist(g)
    el.class <- apply(el, 2, function(diseases) col.vert[diseases])
    el.class.same <- el.class[, 1] == el.class[, 2]
    edge.cols <- rep(edge.col.diff, ecount(g))
    edge.cols[el.class.same] <- el.class[el.class.same, 1]
    edge.cols[edge.cols == col.other] <- edge.col.diff
    E(g)$color <- edge.cols
    
    # plot graph
    plot(g, ...)
    
    invisible(g)
}
