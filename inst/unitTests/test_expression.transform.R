test_expression.transform <- function() {
    # 1) always returns a matrix of values between 0 and 1
    # 2) returns correct values when zero rows not removed
    # 3) returns correct values when zero rows removed
    # 4) works when only one gene input
    
    
    # SETUP
    n.genes <- 4
    n.contexts <- 3
    gene.names <- paste("gene", 1:n.genes)
    context.names <- paste("context", 1:n.contexts)
    expression <- array(c(0,3,1.5,0, 1,0,1.5,0, 2,0,0,0), dim=c(n.genes, n.contexts), dimnames=list(gene.names, context.names))    
    expression.single <- array(expression[1, ], dim=c(1, n.contexts), dimnames=list(gene.names[1], context.names))
    
    # there are the correct results, they have been worked out by hand
    correct <- list()
    correct$zero.row <- array(c(1,6,4,1, 4,1,6,1, 6,2,2,2)/6, dim=c(n.genes, n.contexts), dimnames=list(gene.names, context.names))
    correct$no.zero.row <- array(c(0,4,2, 2,0,4, 4,1,1)/4, dim=c(n.genes - 1, n.contexts), dimnames=list(gene.names[-4], context.names))
    correct$single <- array(c(1,1,1), dim=c(1, n.contexts), dimnames=list(gene.names[1], context.names))
    
    
    # RUN FUNCTION
    results <- list()
    results$zero.row <- expression.transform(expression, remove.zero=FALSE)
    results$no.zero.row <- expression.transform(expression, remove.zero=TRUE)
    results$single <- expression.transform(expression.single, remove.zero=FALSE)
    
    
    # RUN TESTS
    # 1
    for (result in results) {
        checkTrue(max(result) <= 1)
        checkTrue(min(result) >= 0)
        checkTrue(is.matrix(result))
    }
    
    # 2
    checkTrue(identical(dimnames(results$zero.row), dimnames(expression)))
    checkIdentical(results$zero.row, correct$zero.row)
    
    # 3
    checkTrue(!identical(dimnames(results$no.zero.row), dimnames(expression)))
    checkIdentical(results$no.zero.row, correct$no.zero.row)
    
    # 4
    checkIdentical(results$single, correct$single)
}
