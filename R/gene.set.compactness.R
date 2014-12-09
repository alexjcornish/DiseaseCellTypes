gene.set.compactness <- function(
    expression,
    genes,
    g,
    n.perm=10000,
    expression.miss=median(expression),
    rwr.r=0.7, 
    rwr.cutoff=1e-5,
    parallel=NULL,
    verbose=TRUE
) {
    # compute the significance of the compactness of a gene set on a number of expression-modified interactomes
    # expression should be a matrix of expression values, containing a row for each gene and a column for each condition
    # g is a network with named vertices
    
    # setup
    gsc.check.input(expression, genes, g, n.perm, expression.miss, rwr.r, rwr.cutoff, parallel) # check input to the function
    cell.types <- colnames(expression)
    genes.required <- unique(V(g)$name, genes)
    expression <- add.missing.expression(expression, genes.required, expression.miss) # expression now contains all genes and network genes
    g <- add.loops.to.lonely(g) # no nodes can have degree 0
    
    # create results list
    res.names <- c("obs", "perm", "pval")
    res <- sapply(res.names, function(i) NULL, simplify=F)
    class(res) <- "dct"
    
    # randomise the expression data to be used for the permutations
    if (n.perm) {
        expression.perm <- sapply(1:n.perm, function(i) expression[cbind(1:nrow(expression), sample(as.integer(1:ncol(expression)), nrow(expression), replace=T))])
        rownames(expression.perm) <- rownames(expression)
    } else {
        expression.perm <- NULL
    }
    
    if (is.null(parallel)) {
        # run sub function without cluster
        
        res$obs <- sapply(cell.types, function(cell.type) gsc.sub(expression[, cell.type], genes, g, rwr.r, rwr.cutoff))
        if (n.perm) res$perm <- sapply(1:n.perm, function(i) gsc.sub(expression.perm[, i], genes, g, rwr.r, rwr.cutoff))
        
    } else {
        # run sub function over cluster if possible
        
        # setup cluster
        if (verbose) message("establishing cluster of size ", parallel, "... ", appendLF=F)
        cl <- makeCluster(parallel, type="SOCK", verbose=FALSE)       
        clusterExport(cl, list("expression", "cell.types", "genes", "g", "rwr.r", "rwr.cutoff"), envir=environment())
        clusterEvalQ(cl, library(DiseaseCellTypes))
        message("done")
        
        # compute compactness for observed expression
        if (length(cell.types) == 1) {
            # parSapply produces error if cell.types is of length 1
            res$obs <- gsc.sub(expression[, cell.types[1]], genes, g, rwr.r, rwr.cutoff)
            names(res$obs) <- cell.types
        } else {
            res$obs <- parSapply(cl, cell.types, function(cell.type) gsc.sub(expression[, cell.type], genes, g, rwr.r, rwr.cutoff))
        }
        
        # compute compactness for permuted expression
        if (n.perm) {
            if (n.perm == 1) {
                # parSapply produces error if cell.types is of length 1
                res$perm <- gsc.sub(expression[, cell.types[1]], genes, g, rwr.r, rwr.cutoff)
            } else {
                res$perm <- parSapply(cl, 1:n.perm, function(i) gsc.sub(expression.perm[, i], genes, g, rwr.r, rwr.cutoff))
            }
        }
       
        # close cluster
        stopCluster(cl)
    }
    
    # compute p-values
    res$pval <- sapply(res$obs, function(obs) if (is.null(res$perm)) NA else max(mean(res$perm <= obs), 1 / n.perm))
   
    res
}



gsc.check.input <- function(
    expression, 
    genes, 
    g, 
    n.perm, 
    expression.miss, 
    rwr.r, 
    rwr.cutoff, 
    parallel
) {
    # check input to the compactness.sig function
    
    # expression
    if (!is.numeric(expression)) stop("expression is not numeric")
    if (!is.matrix(expression)) stop("expression is not a matrix")
    if (is.null(rownames(expression))) stop("expression does not contain rownames")
    if (is.null(colnames(expression))) stop("expression does not contain colnames")
    
    # g
    if (!is.igraph(g)) stop("g is not an igraph object")
    if (is.null(V(g)$name)) stop("g does not contain any node names")
    
    # genes
    if (!is.character(genes)) stop("genes not a character vector")
    if (!all(genes %in% V(g)$name)) stop("not all genes on network")
    
    # n.perm
    if (n.perm < 0) stop("n.perm should be >= 0")
    
    # expression.miss
    if (!is.numeric(expression.miss)) stop("expression.miss should be numeric")
    
    # rwr.r
    if (rwr.r < 0 | rwr.r > 1) stop("rwr.r should be between 0 and 1")
    
    # rwr.cutoff
    if (rwr.cutoff <= 0) stop("rwr.cutoff should be > 0")
    
    # parallel
    if (!is.null(parallel) & !is.numeric(parallel)) stop("parallel should be NULL or numeric")
}
