test_distance.rwr <- function() {
    # TESTS
    # 1) returns numeric matrix
    # 2) dimnames same as vertex names
    # 3) doesn't return negative, infinite, NA or NULL values
    # 4) diagonal always 0
    
    # 5) works on connected and unconnected graphs
    # 6) works on graph with 0 nodes
    # 7) works on graph with 1 node
    
    # 8) correctly incorporates r
    # 9) correctly incorporates iteration.cutoff
    
    # 10) works when v contains a subset of nodes
    # 11) works when v contains a subset of nodes in a random order
    # 12) works when v is of length 1
    
    # 13) works with toy example without edge weights
    # 14) works with toy example with edge weights
    
    # 15) produces error if g is not an undirected igraph object
    # 16) produces error if edge.attr specified, but edge weights not present
    # 17) produces error if edge weights are not numeric, negative, infinite, or NA
    # 18) produces error if r is not between 0 and 1
    # 19) produces error if iteration cutoff <= 0
    
    
    # SETUP 
    n.genes <- list(node0=0, node1=1, toy.connected=6, toy.unconnected=8)
    edge.attr <- "weights"
    toys <- c("toy.connected", "toy.unconnected")
    toy.edges <- c(1,2, 1,3, 1,4, 1,5, 2,3, 2,4, 2,6, 3,4)
    toy.weights <- c(1, 1, 1, 10, 1, 1, 1, 1)
    
    gs <- list()
    gs <- sapply(n.genes, graph.empty, directed=F, simplify=F)
    gs <- sapply(gs, function(g) set.vertex.attribute(g, "name", value=paste("gene", 1:vcount(g), sep="")), simplify=F)
    gs$node0 <- set.vertex.attribute(gs$node0, "name", value=NULL)
    gs[toys] <- sapply(gs[toys], add.edges, toy.edges, simplify=F)
    gs[toys] <- sapply(gs[toys], set.edge.attribute, name=edge.attr, value=toy.weights, simplify=F)
    
    # produce network with error-producing weights
    error.types <- c("char", "neg", "inf", "na")
    weights.error <- sapply(error.types, function(i) toy.weights, simplify=F)
    weights.error$char <- as.character(weights.error$char)
    weights.error$neg[1] <- -1 
    weights.error$inf[2] <- Inf
    weights.error$na[3] <- NA
    g.error <- gs$toy.connected
    for (type in error.types) g.error <- set.edge.attribute(g.error, type, value=weights.error[[type]])
    
    # values of v to test
    vs <- list(subset=c(2,4,6), subset.strange=c(4,6,2), single=c(4))
    
    
    # RUN FUNCTION 
    results <- list()
    results$weights.without <- sapply(gs, distance.rwr, edge.attr=NULL)
    results$weights.with <- sapply(gs[toys], distance.rwr, edge.attr=edge.attr)
    results$other <- list() 
    results$other$big.r <- distance.rwr(gs$toy.unconnected, edge.attr=edge.attr, rwr.r=0.99)
    results$other$big.it.cutoff <- distance.rwr(gs$toy.unconnected, edge.attr=edge.attr, rwr.cutoff=1e+10)
    
    results.v <- sapply(vs, function(v) distance.rwr(gs$toy.connected, v=v, edge.attr=edge.attr), simplify=F)
    
    
    # RUN TESTS
    for (results.outer in results) {
        for (result in results.outer) {
            # 1
            checkTrue(is.matrix(result)) 
            checkTrue(is.numeric(result)) 
            
            # 2
            result.dimnames <- dimnames(result)
            gene.names <- if(nrow(result) == 0) NULL else paste("gene", 1:nrow(result), sep="")
            checkIdentical(gene.names, rownames(result))
            checkIdentical(gene.names, colnames(result))
            
            # 3
            checkTrue(all(result >= 0)) 
            checkTrue(all(is.finite(result)))
            checkTrue(all(!is.na(result)))
            checkTrue(all(!is.null(result)))
            
            # 4
            checkTrue(all(diag(result) == 0))
        }
    }
    
    # 5, 6, 7
    max.connected <- max(results$weights.with$toy.connected)
    checkTrue(all(results$weights.with$toy.unconnected[, 7:8] %in% c(0, max.connected)))
    checkTrue(all(results$weights.with$toy.unconnected[7:8, ] %in% c(0, max.connected)))
    checkIdentical(results$weights.without$node0, array(0, dim=c(0,0), dimnames=list(NULL, NULL)))
    checkIdentical(results$weights.without$node1, array(0, dim=c(1,1), dimnames=list("gene1", "gene1")))
    
    # 8, 9
    checkTrue(all(results$other$big.r >= results$weights.with$toy.unconnected))
    checkEquals(results$other$big.it.cutoff[1, 4], max(results$other$big.it.cutoff)) # i.e. the walker hasn't moved 2 spaces yet

    # 10, 11, 12
    tol <- 0.01
    checkTrue(all(as.vector(results.v$subset) - as.vector(results$weights.with$toy.connected[vs$subset, ]) < tol))
    checkIdentical(results.v$subset[order(vs$subset), ], results.v$subset.strange[order(vs$subset.strange), ])
    checkTrue(all(as.vector(results.v$single) - as.vector(results.v$subset[match(vs$single, vs$subset), ]) < tol))
    
    # 13, 14
    checkTrue(results$weights.with$toy.connected[1,5] < results$weights.with$toy.connected[2,6])
    checkTrue(results$weights.without$toy.connected[1,5] == results$weights.without$toy.connected[2,6])
    
    # 15, 16, 17, 18, 19
    checkException(distance.rwr(5), silent=T)
    checkException(distance.rwr(erdos.renyi.game(5, 0.2, directed=T)), silent=T)
    checkException(distance.rwr(gs$toy.connected, edge.attr="error"), silent=T)
    for (type in error.types) checkException(distance.rwr(g.error, edge.attr=type), silent=T)
    checkException(distance.rwr(gs$toy.connected, edge.attr=edge.attr, r=-1), silent=T)
    checkException(distance.rwr(gs$toy.connected, edge.attr=edge.attr, r=2), silent=T)
    checkException(distance.rwr(gs$toy.connected, edge.attr=edge.attr, iteration.cutoff=-1), silent=T)
}
