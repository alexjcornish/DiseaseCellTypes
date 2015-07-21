add.missing.expression <- function(
    expression,
    genes,
    score.missing
) {
    # add genes in 'genes' missing in 'expression' to 'expression' with score 'score.missing'
    # remove genes in 'expression' not represented in 'genes'
    
    contexts <- colnames(expression)
    n.contexts <- ncol(expression)
    
    if (is.null(rownames(expression))) stop("genes not found as rownames in expression")
    genes.found <- genes[genes %in% rownames(expression)]
    genes.missing <- genes[!genes %in% rownames(expression)]
    
    expression.data <- array(as.numeric(expression[genes.found, ]), dim=c(length(genes.found), n.contexts), dimnames=list(genes.found, contexts))
    expression.no.data <- array(score.missing, dim=c(length(genes.missing), n.contexts), dimnames=list(genes.missing, contexts))
    rbind(expression.data, expression.no.data)
}
