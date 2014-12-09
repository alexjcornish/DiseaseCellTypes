test_get.adjacency.norm <- function() {
    # TESTS
    # 1) works when edge.attr present
    # 2) works when edge.attr is NULL
    # 3) works with looped lonely nodes
    # 4) outputs unamed sparse matrix
    
    
    # SETUP
    n.genes <- 5
    gene.names <- paste("gene", 1:n.genes)
    edges <- c(1,2, 1,3, 1,4, 2,3, 2,4, 5,5)
    weights <- c(1, 1, 2, 3, 2, 2)
    edge.attr <- "weights"
    
    g <- graph.empty(n.genes, directed=F)
    g <- add.edges(g, edges)
    g <- set.vertex.attribute(g, "name", value=gene.names)
    g <- set.edge.attribute(g, edge.attr, value=weights)
    
    correct <- list()
    correct$with.edge <- c(0,1/4,1/4,2/4,0, 1/6,0,3/6,2/6,0, 1/4,3/4,0,0,0, 2/4,2/4,0,0,0, 0,0,0,0,2/2)
    correct$without.edge <- c(0,1/3,1/3,1/3,0, 1/3,0,1/3,1/3,0, 1/2,1/2,0,0,0, 1/2,1/2,0,0,0, 0,0,0,0,1/1)
    
    
    # RUN FUNCTION
    results <- list()
    results$with.edge <- DiseaseCellTypes:::get.adjacency.norm(g, edge.attr)
    results$without.edge <- DiseaseCellTypes:::get.adjacency.norm(g, NULL)
    
   
    # RUN TESTS
    tol <- 1e-5
    for (name in names(results)) {
        checkEquals(class(results[[name]])[1], "dgCMatrix")
        checkTrue(all(dim(results[[name]]) == as.integer(n.genes)))   
        checkIdentical(dimnames(results[[name]]), list(NULL, NULL))
        checkTrue(all(as.vector(results[[name]]) - correct[[name]] < tol))
    }
}
