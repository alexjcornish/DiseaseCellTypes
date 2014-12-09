gso.sub <- function(
    expression,
    genes
) {
    # compute the mean or median expression of genes
    if (!all(genes %in% names(expression))) warning("not all genes in expression data")
    mean(expression[genes], na.rm=TRUE)
}
