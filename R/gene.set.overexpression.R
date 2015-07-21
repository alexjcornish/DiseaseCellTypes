gene.set.overexpression <- function(
    expression,
    genes,
    n.perm=10000,
    expression.miss=median(expression)
) {
    # compute the significance of the overexpression of a gene set
    
    # setup
    cell.types <- colnames(expression)
    gso.check.input(expression, genes, n.perm, expression.miss)
    expression <- add.missing.expression(expression, genes, expression.miss) # only the genes are required in the output expression matrix

    # create results list
    res.names <- c("obs", "perm", "pval")
    res <- sapply(res.names, function(i) NULL, simplify=F)
    class(res) <- "dct"
    
    # compute observed and permuted values
    res$obs <- sapply(cell.types, function(cell.type) gso.sub(expression[, cell.type], genes))
    if (n.perm) {
        expression.perm <- sapply(1:n.perm, function(i) expression[cbind(1:nrow(expression), sample(as.integer(1:ncol(expression)), nrow(expression), replace=T))])
        if (is.null(dim(expression.perm))) expression.perm <- array(expression.perm, dim=c(1, n.perm), dimnames=list(genes, NULL))
        rownames(expression.perm) <- rownames(expression)
        res$perm <- sapply(1:n.perm, function(i) gso.sub(expression.perm[, i], genes)) 
    }
    res$pval <- sapply(res$obs, function(obs) if (is.null(res$perm)) NA else max(mean(obs <= res$perm), 1 / n.perm))
    res
}



gso.check.input <- function(
    expression, 
    genes, 
    n.perm, 
    expression.miss
) {
    # check the input to the expression.sig function
    
    # expression
    if (!is.matrix(expression)) stop("expression is not a matrix")
    if (!is.numeric(expression)) stop("expression is not numeric")
    if (is.null(rownames(expression))) stop("expression does not contain gene names")
    if (sum(is.na(expression)) > 0) stop("NA values in expression")
    if (sum(is.null(expression)) > 0) stop("NULL values in expression")
    if (sum(!is.finite(expression)) > 0) stop("Infinite values in expression")
    
    # genes
    if (!is.character(genes)) stop("genes in not a character vector")
    if (sum(is.na(genes)) > 0) stop("NA values in genes")
    if (sum(is.null(genes)) > 0) stop("NULL values in genes")
    if (sum(genes %in% rownames(expression)) == 0) warning("no genes in expression")
    
    # n.perm
    if (n.perm < 0) stop("n.perm is less than 0")
    if (!is.numeric(n.perm)) stop("n.perm is not numeric")
    
    # expression missing
    if (expression.miss < 0) stop("expression.miss is less that 0")
    if (!is.numeric(expression.miss)) stop("expression.miss is not numeric")
}
