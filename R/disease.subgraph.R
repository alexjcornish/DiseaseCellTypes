disease.subgraph <- function(
    g,
    genes,
    edge.attr="score",
    vert.attr="expression",
    filename.plot=NULL,
    filename.legend.vert=NULL,
    filename.legend.edge=NULL,
    n.bins=500,
    
    # size parameters
    vert.size.disease=6,
    vert.size.not.disease=3,
    edge.width.max=20,
    
    # color parameters
    vert.col.high="black",
    vert.col.low="lightgrey",
    edge.col.high="blue",
    edge.col.low="lightgrey",
    
    ...
) {
    # description
    
    # setup
    if (!vert.attr %in% list.vertex.attributes(g)) stop("vert.attr not found on network")
    if (!edge.attr %in% list.edge.attributes(g)) stop("edge.attr not found on network")
    
    # create a subnetwork containing only the interactors
    genes.interactors <- get.interacting.genes(g, genes)
    g.subgraph <- induced.subgraph(g, sample(c(genes, genes.interactors)))
    genes.subgraph <- V(g.subgraph)$name
    
    # change vertex size and remove labels
    vert.size <- rep(vert.size.not.disease, vcount(g.subgraph))
    vert.size[V(g.subgraph)$name %in% genes] <- vert.size.disease
    V(g.subgraph)$size <- vert.size
    V(g.subgraph)$label <- rep(NA, vcount(g.subgraph))
    
    # shape vertices by disease gene status
    V(g.subgraph)$shape <- rep("circle", vcount(g.subgraph))
    V(g.subgraph)$shape[V(g.subgraph)$name %in% genes] <- "square"
    
    # color vertices by expression scores
    expression.scores <- get.vertex.attribute(g.subgraph, vert.attr)
    vert.breaks <- c(0, seq(0.50, 1, length.out=n.bins))
    vert.bin.number <- cut(expression.scores, breaks=vert.breaks, labels=FALSE, include.lowest=TRUE)
    vert.col.range <- colorRampPalette(c(vert.col.low, vert.col.high))(n.bins)
    V(g.subgraph)$color <-  V(g.subgraph)$frame.color <- vert.col.range[vert.bin.number]
     
    # weight and color edges by score
    edge.col.range <- colorRampPalette(c(edge.col.low, edge.col.high))(n.bins)
    edge.breaks <- c(0, seq(0.25, 1, length.out=n.bins))
    edge.bin.number <- cut(get.edge.attribute(g.subgraph, edge.attr), breaks=edge.breaks, labels=FALSE)
    E(g.subgraph)$width <- edge.bin.number * edge.width.max / n.bins
    E(g.subgraph)$color <- edge.col.range[edge.bin.number]
    
    # produce the plot
    if (!is.null(filename.plot)) {
        pdf(filename.plot, ...)
        plot(g.subgraph, layout=layout.fruchterman.reingold)
        dev.off()
    }
    
    # plot the legends
    if (!is.null(filename.legend.vert)) {
        pdf(filename.legend.vert, height=5, width=10)
        image(z=matrix(vert.breaks, ncol=1), col=vert.col.range, breaks=vert.breaks, yaxt="n")
        dev.off()
    }
    if (!is.null(filename.legend.edge)) {
        pdf(filename.legend.edge, height=5, width=10)
        image(z=matrix(edge.breaks, ncol=1), col=edge.col.range, breaks=edge.breaks, yaxt="n")
        dev.off()
    }
    
    # return the mean number of weights greater than the cutoff
    g.subgraph
}

get.interacting.genes <- function(g, genes) {
    el <- get.edgelist(g)
    el.genes <- el[el[, 1] %in% genes | el[, 2] %in% genes, ]
    inter <- as.vector(el.genes)
    inter <- inter[!inter %in% genes] # remove those interactors in the original genes
    unique(inter)
}
