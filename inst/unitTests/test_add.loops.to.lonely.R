test_add.loops.to.lonely <- function() {
    # TESTS
    # 1) works if g contains no nodes
    # 2) works if g contains no lone nodes
    # 3) works if g contains some lone nodes
    # 4) works if g contains only lonely nodes (adds edge.attr if required)
    # 5) error if g is not an igraph object
    # 6) error if edge.attr specified but not on g
    # 7) error if edge.attr specified but not numeric (say not yet implemented)
    # 8) error if loop.weight is not numeric 
    
    
    # SETUP
    n.genes.base <- 10
    n.genes.lone <- 3
    edge.attr <- "weight"
    edge.attr.char <- "error"
    loop.weight <- 10
    
    gs <- list()
    gs$node0 <- graph.empty(0, directed=F)
    gs$loneF <- barabasi.game(n.genes.base, 2, directed=F)
    gs$loneF <- set.edge.attribute(gs$loneF, edge.attr, value=runif(ecount(gs$loneF)))
    gs$loneF <- set.edge.attribute(gs$loneF, edge.attr.char, value=sample(LETTERS[1:26], ecount(gs$loneF), replace=F))
    gs$loneT <- add.vertices(gs$loneF, n.genes.lone)
    gs$loneALL <- graph.empty(n.genes.base, directed=F)
    gs <- sapply(gs, function(g) if(vcount(g)) set.vertex.attribute(g, "name", value=paste("gene", 1:vcount(g))) else g, simplify=F)
    

    # RUN FUNCTION
    results <- sapply(gs, DiseaseCellTypes:::add.loops.to.lonely, edge.attr=edge.attr, loop.weight=loop.weight, simplify=F)
    

    # RUN TESTS
    # 1, 2, 3, 4
    for (g.name in names(gs)) {      
        checkTrue(is.igraph(results[[g.name]]))
        checkTrue(all(list.vertex.attributes(gs[[g.name]]) %in% list.vertex.attributes(results[[g.name]])))
        checkTrue(all(list.edge.attributes(gs[[g.name]]) %in% list.edge.attributes(results[[g.name]])))
        checkTrue(all(degree(results[[g.name]]) > 0))
        checkEquals(ecount(simplify(results[[g.name]])), ecount(gs[[g.name]]))
        
        if (vcount(gs[[g.name]])) {
            # if there are any vertices in g, check that all new edges are of weight loop.weight
            edge.weights.g <- get.edge.attribute(gs[[g.name]], edge.attr)
            edge.weights.result <- get.edge.attribute(results[[g.name]], edge.attr)
            checkTrue(all(edge.weights.result[!edge.weights.result %in% edge.weights.g] == loop.weight))
        }
    }
    
    checkException(DiseaseCellTypes:::add.loops.to.lonely(5), silent=T) # 5
    checkException(DiseaseCellTypes:::add.loops.to.lonely(gs$loneF, edge.attr="not_present"), silent=T) # 6
    checkException(DiseaseCellTypes:::add.loops.to.lonely(gs$loneF, edge.attr=edge.attr.char), silent=T) # 7
    checkException(DiseaseCellTypes:::add.loops.to.lonely(gs$loneF, loop.weight="1"), silent=T) # 8
}
