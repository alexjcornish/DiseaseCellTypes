test_gene.set.compactness <- function() {
    # TESTS
    # 1) works on toy example
    # 2) works with lone vertices
    # 3) works when expression values not present for all network genes
    # 4) works with expression values equal to 0
    # 5) can run without permutations
    # 6) works parallel where available 
   
    # SETUP
    n.genes <- 12
    gene.names <- paste("gene", 1:n.genes)
    
    n.contexts <- 10
    context.names <- paste("context", 1:n.contexts)
    
    gene.mixed <- gene.names[1]
    context.high <- context.names[1]
    context.low <- context.names[!context.names %in% context.high]
    
    genes.missing <- gene.names[7:8]
    expression.miss <- 0.001
    
    expression <- list()
    expression$main <- array(runif(n.genes * n.contexts), dim=c(n.genes, n.contexts), dimnames=list(gene.names, context.names))
    expression$main[gene.mixed, context.high] <- expression$main[gene.mixed, context.high] / 100
    expression$main[gene.mixed, context.high] <- 100
    expression$main[genes.missing, ] <- expression.miss
    expression$main[gene.names[5], context.names[5]] <- 0
    expression$some.missing <- expression$main[rownames(expression$main)[!rownames(expression$main) %in% genes.missing], ] 
    
    genes <- gene.names[c(2,3)]
    
    edges <- c(1,2, 1,3, 2,4, 2,5, 2,6, 2,7, 3,8, 3,9, 3,10, 3,11, 4,5, 5,6, 6,7, 8,9, 9,10, 10,11)
    g <- graph.empty(n.genes, directed=F)
    g <- set.vertex.attribute(g, "name", value=gene.names)
    g <- add.edges(g, edges)
    
    n.perm <- 100
    
    
    # RUN FUNCTION
    results <- list()
    results$main <- gene.set.compactness(expression$main, genes, g, n.perm) # 1, 2)
    results$some.missing <- gene.set.compactness(expression$some.missing, genes, g, n.perm, expression.miss=expression.miss) # 3)
    results$no.perm <- gene.set.compactness(expression$main, genes, g, n.perm=0) # 5)
    results$parallel <- gene.set.compactness(expression$main, genes, g, n.perm, parallel=2) # 6)

    
    # RUN TESTS
    for (result in results) {
        checkTrue(is.list(result))
        checkTrue(all(c("obs", "perm", "pval") %in% names(result)))
    }
    
    # 1, 2, 4)
    checkTrue(all(results$main$obs[context.high] <= results$main$obs[context.low]))
    checkTrue(all(results$main$pval[context.high] <= results$main$pval[context.low]))
    
    # 3)
    checkIdentical(results$main$obs, results$some.missing$obs)
    
    # 5) 
    checkTrue(is.null(results$no.perm$perm))
    checkTrue(all(is.na(results$no.perm$pval)))
    
    # 6) 
    checkIdentical(results$main$obs, results$parallel$obs)
}
