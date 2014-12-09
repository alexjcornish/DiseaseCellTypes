expression.transform <- function(
    x, 
    remove.zero=TRUE
) {
    # transform raw expression values x to relative expression scores and percentile transform
    
    # setup
    if (!is.numeric(x)) stop("x is not numeric")
    if (!is.matrix(x)) stop("x is not a matrix")
    if (ncol(x) < 2) stop("more than once context required")
    if (sum(is.na(x)) > 0) stop("NA values present in x")
    if (remove.zero & all(x == 0)) stop("all x values is 0")
    
    # remove rows where all values equal 0 if required
    if (remove.zero) {
        dimnames.x <- dimnames(x)
        zero.rows <- apply(x, 1, function(row) all(row == 0))
        x <- x[!zero.rows , ]
        
        # if only 1 row contains non-zero values, ensure that it remains a matrix
        if (is.null(dim(x))) x <- x.to.matrix(x, dimnames.x[[1]][!zero.rows], dimnames.x[[2]])
    }
    dimnames.x <- dimnames(x)
    
    # transform the raw values to relative values
    row.mean <- apply(x, 1, mean, na.rm=T)
    row.mean[row.mean == 0] <- 1
    x <- x / row.mean
    
    # percentile-transform expression values
    # genes have a score of 1 if they are the most overexpressed gene in a context and a score of 0 if they are the most underexpressed
    # if x contains only 1 row, then all values are given a percentile score of 1
    x <- apply(x, 2, rank, ties.method="average")
    if (is.null(dim(x))) x <- x.to.matrix(x, dimnames.x[[1]], dimnames.x[[2]])
    if (nrow(x) > 1) x <- (x - 1) / (nrow(x) - 1)
    
    # ensure dimnames are present and output
    dimnames(x) <- dimnames.x
    x
}

x.to.matrix <- function(x, row.names, col.names) array(x, dim=c(1, length(x)), dimnames=list(row.names, col.names)) 
