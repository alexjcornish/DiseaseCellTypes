test_add.missing.expression <- function() {
    # TESTS
    # 1) add genes in genes missing from expression
    # 2) removes genes from expression missing from genes
    # 3) added genes are given scores of expression.missing
    # 4) works when only a single context input
    # 5) works when all genes from genes are in expression
    # 6) works when no genes from genes are in expression
    # 7) works when contexts not named
    # 8) works when only a single gene is input and found
    # 9) works when only a single gene is intput and not found
    # 10) outputs numeric matrix with gene names
    
    
    # SETUP
    expression.names <- c("all", "single", "no.context")
    n.genes <- c(4, 4, 4)
    n.contexts <- c(3, 1, 3)
    context.names <- c(T, T, F)
    names(n.genes) <- names(n.contexts) <- names(context.names) <- expression.names
    expression <- sapply(expression.names, function(name) array(runif(n.genes[name] * n.contexts[name]), dim=c(n.genes[name], n.contexts[name]), dimnames=list(paste("gene", 1:n.genes[name]), if (context.names[name]) paste("context", 1:n.contexts[name]) else NULL)), simplify=F)
    gene0 <- "gene 1"
    context0 <- "context 1"
    expression$all[gene0, context0] <- 0
    
    genes.numbers <- list(mixed=c(1, 2, 3, 11, 12, 13), all=c(1, 2, 3), none=c(11, 12, 13))
    genes <- sapply(genes.numbers, function(i) paste("gene", i), simplify=F)

    score.missing <- 1
    
    
    
    # RUN FUNCTION
    results <- list()
    results$all.mixed <- DiseaseCellTypes:::add.missing.expression(expression$all, genes$mixed, score.missing) # 1, 2, 3
    results$single.mixed <- DiseaseCellTypes:::add.missing.expression(expression$single, genes$mixed, score.missing) # 4
    results$all.all <- DiseaseCellTypes:::add.missing.expression(expression$all, genes$all, score.missing) # 5
    results$all.none <- DiseaseCellTypes:::add.missing.expression(expression$all, genes$none, score.missing) # 6
    results$no.context.mixed <- DiseaseCellTypes:::add.missing.expression(expression$no.context, genes$mixed, score.missing) # 7
    results$gene.single.found <- DiseaseCellTypes:::add.missing.expression(expression$all, "gene 1", score.missing) # 8
    results$gene.single.not.found <- DiseaseCellTypes:::add.missing.expression(expression$all, "gene 11", score.missing) # 9
    

    
    # RUN TESTS
    # 1, 2, 3
    checkTrue(all(genes$mixed %in% rownames(results$all.mixed)))
    checkTrue(all(rownames(results$all.mixed %in% genes$mixed)))
    checkTrue(all(as.vector(results$all.mixed[!rownames(results$all.mixed) %in% genes$mixed, ]) == score.missing))
    
    # 4
    checkEquals(rownames(results$single.mixed), genes$mixed)
    checkEquals(colnames(results$single.mixed), colnames(expression$single))

    # 5
    checkEquals(dimnames(results$all.all), list(genes$all, colnames(expression$all)))
    
    # 6 
    checkTrue(all(as.vector(results$all.none) == score.missing))
    checkEquals(dimnames(results$all.none), list(genes$none, colnames(expression$all)))
    
    # 7 
    checkEquals(rownames(results$no.context.mixed), genes$mixed)
    checkEquals(ncol(results$no.context.mixed), ncol(expression$no.context), checkNames=F)
    checkTrue(is.null(colnames(results$no.context.mixed)))
    
    # 8 
    checkEquals(rownames(results$gene.single.found), "gene 1")
    checkIdentical(colnames(results$gene.single.found), colnames(expression$all))
    checkTrue(all(as.vector(results$gene.single.found) == expression$all["gene 1", ]))
    
    # 9
    checkEquals(rownames(results$gene.single.not.found), "gene 11")
    checkIdentical(colnames(results$gene.single.not.found), colnames(expression$all))
    checkTrue(all(as.vector(results$gene.single.not.found) == score.missing))
    
    # 10
    for (result in results) {
        checkTrue(is.matrix(result))
        checkTrue(is.numeric(result))
        checkTrue(!is.null(rownames(result)))
    }
}
