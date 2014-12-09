test_gso.sub <- function() {
    # TESTS
    # 1) works correctly
    
    # SETUP
    n.genes <- 10
    n.genes.sample <- 3
    gene.names <- paste("gene", 1:n.genes)
    gene.set <- sample(gene.names, n.genes.sample)
    expression <- structure(runif(n.genes), names=gene.names)
    
    
    # RUN FUNCTION
    result <- DiseaseCellTypes:::gso.sub(expression, gene.set)
    
    
    # RUN TESTS
    checkEquals(result, mean(expression[gene.set])) # 1)
}
