plot.dct <- function(
    x,
    contexts=NULL,
    cutoff=0.05, 
    include.pval=FALSE, # if TRUE, the p-values are added to the end of the context names
    ...
) {
    # plot the output of the gene.set.compactness or gene.set.overexpression function
    # may want to change the class name when I've got the package name finalised
    # if contexts is NULL, only the contexts with p-value <= cutoff are plotted
    # if contexts is not NULL, only the contexts specified are plotted and the cutoff is ignored
    # these contexts should be used to name x$obs
    
    # plotting details
    xlab <- "Expression sum"
    main <- "Expression: observed v. permuted"
    border <- "grey"
    col <- "grey"
    
    # setup
    if (class(x) != "dct") stop("plot.dct should be applied to objects of class dct")
    if (!"obs" %in% names(x)) stop("observed values not found in x")
    if (!"perm" %in% names(x)) stop("permuted values not found in x")
    if (!"pval" %in% names(x)) stop("pvalues not found in x")
    contexts.x <- names(x$obs)
    if (!is.null(contexts)) {
        if (is.null(contexts.x)) stop("observed values in x do not have context names")
        if (!all(contexts %in% contexts.x)) stop("not all contexts in x")
    }
    
    # identify the contexts to plot
    contexts.plot <- if (is.null(contexts) & sum(is.na(x$pval)) == 0) contexts.x[x$pval <= cutoff] else contexts
    
    # identify the values to plot    
    obs <- x$obs[contexts.plot]
    perm <- x$perm
    
    # compute range to plot
    if (length(c(obs, perm)) == 0) {
        # if there are no values to plot, plot empty graph
        xlim <- c(0, 1)
    } else {
        value.min <- min(c(obs, perm))
        value.max <- max(c(obs, perm))
        value.d <- value.max - value.min
        xlim <- if (value.d == 0) c(0, max(1, value.max + 0.5)) else c(max(0, value.min - value.d * 0.1), value.max + value.d * 0.1)
    }
    
    # create text to plot for each context
    if (length(contexts.plot)) {
        contexts.text <- sapply(contexts.plot, function(context) paste(context, ifelse(include.pval, paste(", p=", x$pval[context], sep=""), ""), sep=""), simplify=T)
        text.max.length <- max(str_length(contexts.text))
    } else {
        contexts.text <- NULL
        text.max.length <- 0
    }
    
    # plotting parameters
    par(mfrow=c(1, 1), mar=c(5,4,text.max.length * 0.45,2) + 0.1, oma=c(0,0,2,0))
    
    if (!is.null(perm)) {
        # if permutations are included
        hist.output <- hist(perm, xlim=xlim, col=col, border=border, xlab=xlab, main="", ...)
    } else {
        # if permutations are not included, don't plot histogram
        neg <- min(xlim) - 1000
        hist.output <- hist(neg, xlim=xlim, breaks=c(neg - 1, neg + 1), col=col, border=border, xlab=xlab, main="", ...)
    }
    title(main, outer=T)
    
    # plot lines 
    par(xpd=FALSE)
    for (context in contexts.plot) abline(v=obs[context], lwd=2, col="red")
    
    # plot text above lines
    par(xpd=NA)
    count.max <- max(pretty(c(0, max(hist.output$counts))))
    for (context in contexts.plot) text(x=obs[context], y=count.max * 1.02, labels=contexts.text[context], srt=90, pos=4, offset=0)
}
