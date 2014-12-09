add.missing.expression <- function(
    expression,
    genes,
    score.missing
) {
    # add genes in 'genes' missing in 'expression' to 'expression' with score 'score.missing'
    # remove genes in 'expression' not represented in 'genes'
    n.contexts <- ncol(expression)
    contexts <- colnames(expression)
    missing.genes <- genes[!genes %in% rownames(expression)]
    expression <- expression[rownames(expression) %in% genes, ] # remove genes not in genes
    if (!is.matrix(expression)) expression <- array(expression, dim=c(length(expression)/length(contexts), length(contexts)), dimnames=list(names(expression), contexts))
    rbind(expression, array(score.missing, dim=c(length(missing.genes), n.contexts), dimnames=list(missing.genes, contexts)))
}
