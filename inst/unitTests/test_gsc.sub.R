test_gsc.sub <- function() {
    # TESTS
    # 1) works correctly
    # 2) correctly incorporates rwr.r
    # 3) correctly incorporates rwr.cutoff
    
    
    # SETUP
    n.genes <- 7
    gene.names <- paste("gene", 1:n.genes)
    edges <- c(1,2, 1,3, 1,4, 1,5, 2,3, 2,6, 3,6, 4,5, 4,7, 5,7)
    n.perms <- 50
    
    expression <- list(low1=rep(1, n.genes), high236=rep(1, n.genes))
    expression <- sapply(expression, function(i) structure(i, names=gene.names), simplify=F)
    expression$low1[1] <- 0.01
    
    genes <- list(a=c(2,3,6), a.mixed=c(3,6,2), b=c(4,5,7), c=c(2,4), d=c(1,2))
    genes <- sapply(genes, function(set) gene.names[set], simplify=F)
    
    g <- graph.empty(n.genes, directed=F)
    g <- set.vertex.attribute(g, "name", value=gene.names)
    g <- add.edges(g, edges)
    
    rwr.r <- 0.7
    rwr.cutoff <- 10e-5
    
    
    # RUN FUNCTION
    results <- list()
  
    # 1)
    results$side1 <- DiseaseCellTypes:::gsc.sub(expression$low1, genes$a, g, rwr.r, rwr.cutoff)
    results$side1.mix <- DiseaseCellTypes:::gsc.sub(expression$low1, genes$a.mixed, g, rwr.r, rwr.cutoff)
    results$side2 <- DiseaseCellTypes:::gsc.sub(expression$low1, genes$b, g, rwr.r, rwr.cutoff)
    results$side12 <- DiseaseCellTypes:::gsc.sub(expression$low1, genes$c, g, rwr.r, rwr.cutoff)
    
    # 2)
    results$adj <- DiseaseCellTypes:::gsc.sub(expression$low1, genes$d, g, rwr.r, rwr.cutoff)
    results$adj.high.rwr.r <- DiseaseCellTypes:::gsc.sub(expression$low1, genes$d, g, 1, rwr.cutoff)
    
    # 3)
    results$side12.high.rwr.cutoff <- DiseaseCellTypes:::gsc.sub(expression$low1, genes$c, g, rwr.r, 999999999)
    
 
    # RUN TESTS
    # 1)
    checkEquals(results$side1, results$side1.mix)
    checkEquals(results$side1, results$side2)
    checkTrue(results$side1 < results$side12)
    
    # 2) 
    checkTrue(results$adj.high.rwr.r > results$adj)
    
    # 3)
    checkTrue(results$side12.high.rwr.cutoff > results$side12)
}
