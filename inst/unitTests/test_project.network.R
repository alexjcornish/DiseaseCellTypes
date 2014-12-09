test_project.network <- function() {
    # TESTS
    # 1) works with a simple matrix - needs to be inspected visually
    # 2) produces error if not all rownames in class.cols
    
    # SETUP
    # create p-values with rows 1/2 and 4/5 very similar to produce known graph
    n.rows <- 6
    n.cols <- 6
    rows <- paste("row", 1:n.rows)
    cols <- paste("col", 1:n.cols)
    pairs <- cbind(template=c(1,4), change=c(2,5))

    p.values <- array(runif(n.rows * n.cols), dim=c(n.rows, n.cols), dimnames=list(rows, cols))
    for (i in 1:nrow(pairs)) p.values[pairs[i, "change"], ] <- p.values[pairs[i, "template"], ] + (runif(n.cols) - 0.5) / 100000
    p.values[p.values < 0] <- 0
    p.values[p.values > 1] <- 1
    p.values <- 10 ^ (-p.values * 10) # ensure that p-values are not uniformly distributed
    
    col.vertices <- structure(rep(rainbow(n.rows / 2), each=2), names=rows)
    
    
    # RUN FUNCTION
    project.network(p.values, col.vertices, rank.max=1) # 1)
    
    
    # RUN TESTS
    checkException(project.network(p.values, col.vertices[1:3]), silent=T) # 2)
}
