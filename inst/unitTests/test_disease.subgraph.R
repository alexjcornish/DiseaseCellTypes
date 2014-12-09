test_disease.subgraph <- function() {
    # TESTS
    # 1) ensure that the function runs correctly
    # 2) produces error if edge.attr or vert.attr not found
    
    # SETUP
    
    # create cell type-specific interactome
    data(edgelist.string)
    data(expression.fantom5)
    g <- graph.edgelist(as.matrix(edgelist.string[, c("idA", "idB")]), directed=FALSE)
    expression <- expression.transform(expression.fantom5)
    g.myoblast <- score.edges(g, expression[, "myoblast"])
    
    # simulate disease genes
    genes.disease <- sample(V(g.myoblast)$name, 5) 
    
    
    # RUN FUNCTION
    res <- disease.subgraph(g.myoblast, genes.disease)
    
    
    # RUN TESTS
    plot(res) # 1)
    checkException(disease.subgraph(g.myoblast, gene.disease, vert.attr="error"), silent=T) # 2
    checkException(disease.subgraph(g.myoblast, gene.disease, edge.attr="error"), silent=T) # 2   
}
