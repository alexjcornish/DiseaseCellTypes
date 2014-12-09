test_score.edges <- function() {
    # TESTS
    # 1) works correctly when answer known
    # 2) produces weights rather than distances when specified
    # 3) corrects missing values when specified
    # 4) ensures that no Inf values are returned when change.zero==T
    # 5) produces error if network contains no gene names
    # 6) produces error if expression contains no gene names
    
    # SETUP
    n.genes <- 4
    gene.names <- paste("gene", 1:n.genes)
    edges <- c(1,2, 2,3, 3,4, 4,1)
    
    g <- graph.empty(n.genes, directed=F)
    g <- add.edges(g, edges)
    g.blank <- g
    V(g)$name <- gene.names
    
    gene.expression <- 0:(n.genes - 1)
    names(gene.expression) <- gene.names
    
    edge.attr <- "score"
    
    
    # RUN FUNCTION
    results <- list()
    results$distance <- score.edges(g, gene.expression, edge.attr.type="distance", edge.attr=edge.attr)
    results$weight <- score.edges(g, gene.expression, edge.attr.type="weight", edge.attr=edge.attr)
    results$missing <- score.edges(g, gene.expression[-1], edge.attr.type="distance", edge.attr=edge.attr, expression.miss=gene.expression[1])
    results$inf <- score.edges(g, gene.expression, edge.attr.type="distance", edge.attr=edge.attr, correct.inf=FALSE)
    
    
    # SETUP TESTS 
    mod.expression <- mod.expression.inf <- gene.expression
    mod.expression[mod.expression == 0] <- min(mod.expression[mod.expression != 0])
    correct.distances <- 1 / (mod.expression * c(mod.expression[-1], mod.expression[1]))
    correct.distances.r.big <- 1 / (mod.expression * c(mod.expression[-1], mod.expression[1]))
    correct.distances.inf <- 1 / (mod.expression.inf * c(mod.expression.inf[-1], mod.expression.inf[1]))
    
    
    # RUN TESTS
    checkIdentical(get.edge.attribute(results$distance, edge.attr), as.numeric(correct.distances)) # 1
    checkIdentical(get.edge.attribute(results$distance, edge.attr), 1 / get.edge.attribute(results$weight, edge.attr)) # 2
    checkTrue(identical(results$distance, results$missing)) # 3
    checkIdentical(get.edge.attribute(results$inf, edge.attr), as.numeric(correct.distances.inf)) # 4
    checkException(res <- score.edges(g.blank, gene.expression, edge.attr), silent=T) # 5
    checkException(res <- score.edges(g, as.numeric(gene.expression), edge.attr), silent=T) # 6
}
