gso.sub <- function(
    expression,
    genes
) {
    # compute the mean or median expression of genes
    if (length(expression) == 1) names(expression) <- genes
    if (!all(genes %in% names(expression))) warning("not all genes in expression data")
    mean(expression[genes], na.rm=TRUE)
}
