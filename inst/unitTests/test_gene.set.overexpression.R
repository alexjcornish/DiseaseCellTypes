test_gene.set.overexpression <- function() {
    # TESTS
    # 1) works with toy examples
    # 2) works when expression values not present for all network genes
    # 3) works with expression values equal to 0
    # 4) can run without permutations

    
    # SETUP
    set.seed(1)
    n.genes.random <- 10
    n.genes.select <- 10
    n.genes <- n.genes.random + n.genes.select
    n.contexts <- 50
    n.contexts.high <- 2
    n.perm <- 1000
    
    gene.names <- paste("gene", 1:n.genes)
    genes <- list()
    genes$select <- paste("gene", (n.genes.random + 1):(n.genes))
    genes$select.missing <- paste("gene", (n.genes.random + 1):(n.genes + 3))
    genes$random <- paste("gene", 1:n.genes.random)
    
    context.names <- paste("context", 1:n.contexts)
    contexts.high <- paste("context", 1:n.contexts.high)
    contexts.low <- context.names[!context.names %in% contexts.high]
    
    expression <- array(runif(n.genes * n.contexts), dim=c(n.genes, n.contexts), dimnames=list(gene.names, context.names))
    expression[genes$select, contexts.high] <- expression[genes$select, contexts.high] * 1000
    
    expression <- expression.transform(expression)
    
    
    # RUN FUNCTION
    results <- list()
    results$select <- gene.set.overexpression(expression, genes$select, n.perm) # 1
    results$select.missing <- gene.set.overexpression(expression, genes$select.missing, n.perm) # 2
    results$random <- gene.set.overexpression(expression, genes$random, n.perm) # 1, 3
    results$no.perm <- gene.set.overexpression(expression, genes$select, n.perm=0) # 4
    
    
    # RUN TESTS
    for (result in results) {
        checkTrue(is.list(result))
        checkTrue(all(c("obs", "perm", "pval") %in% names(result)))
    }
    
    # 1)
    checkTrue(all(results$select$obs[contexts.high] > results$select$obs[contexts.low]))
    checkTrue(all(results$select$pval[contexts.high] < results$select$pval[contexts.low]))  
    
    # 2) 
    checkTrue(all(results$select.missing$obs[contexts.high] > results$select.missing$obs[contexts.low]))
    checkTrue(all(results$select.missing$pval[contexts.high] < results$select.missing$pval[contexts.low]))  
    
    # 4)
    checkTrue(is.null(results$no.perm$perm))
    checkTrue(all(is.na(results$no.perm$pval)))
}
