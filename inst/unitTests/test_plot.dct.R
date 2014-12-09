test_plot.dct <- function() {
    # TESTS
    # 1) works with both functions
    # 2) works with 0 permutations
    # 3) works with many permutations
    # 4) works when contexts specified
    # 5) includes values below a certain p-value
    # 6) adds p-value names if required
 
    
    # SETUP
    n.genes <- 12
    n.contexts <- 10
    gene.names <- paste("gene", 1:n.genes)
    context.names <- paste("context", 1:n.contexts)
    
    expression <- array(runif(n.genes * n.contexts), dim=c(n.genes, n.contexts), dimnames=list(gene.names, context.names))
    genes <- gene.names[c(2,3)]
    
    edges <- c(1,2, 1,3, 2,4, 2,5, 2,6, 2,7, 3,8, 3,9, 3,10, 3,11, 4,5, 5,6, 6,7, 8,9, 9,10, 10,11)
    g <- graph.empty(n.genes, directed=F)
    g <- set.vertex.attribute(g, "name", value=gene.names)
    g <- add.edges(g, edges)

    n.perms <- as.character(c(0, 100))
    
    results <- list()
    results$expression <- sapply(n.perms, function(n.perm) gene.set.overexpression(expression, genes, n.perm=as.numeric(n.perm)), simplify=F)
    results$compactness <- sapply(n.perms, function(n.perm) gene.set.compactness(expression, genes, g, n.perm=as.numeric(n.perm)), simplify=F)
    
    
    # RUN FUNCTIONS - need to be checked visually  
    # 1, 2)
    plot(results$expression$'0', cutoff=0.5)
    plot(results$compactness$'0', cutoff=0.5)
    
    # 1, 3)
    plot(results$expression$'100', cutoff=0.5)
    plot(results$compactness$'100', cutoff=0.5)
    
    # 4)
    plot(results$expression$'100', cutoff=0.5, context="context 3")
    
    # 5)
    plot(results$expression$'100', cutoff=1)
    
    # 6)
    plot(results$expression$'100', cutoff=1, include.pval=T)
}
