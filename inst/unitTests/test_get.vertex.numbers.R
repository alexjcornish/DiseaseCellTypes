test_get.vertex.numbers <- function() {
    # TESTS 
    # 1) works with input as igraph.vs object, character and numeric vector (output as integer in the same order as input)
    # 2) works when input is NULL (NULL is output)
    # 3) error of g is not an igraph object
    # 4) error if v is not an igraph.vs object, character or numeric vector
    # 5) error if v is igraph.vs object and not all v in g
    # 6) error if v is character vector and not all v in g
    # 7) error if v is numeric vector and not all v in g
    # 8) error if v is character vector and g does not have names
    
    
    # SETUP
    n.genes <- 10
    n.correct.nums <- 4
    correct.nums <- as.integer(sample(n.genes, n.correct.nums))
    
    gene.names <- paste("gene", 1:n.genes)
    g <- g.no.names <- graph.empty(n.genes, directed=F)
    g <- set.vertex.attribute(g, "name", value=gene.names)

    vs <- list(igraph.vs=V(g)[correct.nums], char=gene.names[correct.nums], num=correct.nums)
    v.large <- V(graph.empty(n.genes + 1, directed=F))[n.genes + 1]
    
    
    # RUN FUNCTION
    results <- sapply(vs, function(v) DiseaseCellTypes:::get.vertex.numbers(g, v), simplify=F)
    
    
    # RUN TESTS
    for (result in results) {
        checkTrue(is.integer(result)) # 1
        checkIdentical(result, correct.nums) # 1
    }
    
    checkTrue(is.null(DiseaseCellTypes:::get.vertex.numbers(g, NULL))) # 2 
    checkException(DiseaseCellTypes:::get.vertex.numbers(4, vs[[1]]), silent=T) # 3
    checkException(DiseaseCellTypes:::get.vertex.numbers(g, TRUE), silent=T) # 4
    checkException(DiseaseCellTypes:::get.vertex.numbers(g, v.large), silent=T) # 5
    checkException(DiseaseCellTypes:::get.vertex.numbers(g, "error"), silent=T) # 6
    checkException(DiseaseCellTypes:::get.vertex.numbers(g, n.genes + 1), silent=T) # 7
    checkException(DiseaseCellTypes:::get.vertex.numbers(g.no.names, vs$char), silent=T) # 8
}
