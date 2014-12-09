project.network.legend <- function(
    col.classes,
    col.other="white", 
    pt.cex=3,
    ... # passed to the legend function
) {
    # plot a legend to a series of colours
    # if multiple catagories have the col.other colour, then they are combined under an 'Other' label
    
    if (!is.null(col.other)) {
        catagories.to.combine <- col.classes[col.classes == col.other]
        col.classes <- col.classes[!col.classes %in% catagories.to.combine]
    }
    col.classes <- col.classes[order(names(col.classes))]
    if (length(catagories.to.combine)) col.classes <- c(col.classes, Other=col.other)
    
    # plot legend on blank plot
    par(oma=rep(0, 4), mar=rep(0, 4))
    plot.new()
    legend("center", legend=names(col.classes), col=col.classes, lty=0, pch=20, bty="n", pt.cex=pt.cex, ...)
    par(mar=c(5,4,4,2) + 0.1)
}
