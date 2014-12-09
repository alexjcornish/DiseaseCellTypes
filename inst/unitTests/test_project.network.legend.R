test_project.network.legend <- function() {
    # TESTS
    # 1) works a set of colours
    # 2) works with 'other' colours
    
    # SETUP
    n.classes <- 9
    n.other <- 3
    col.other <- "white"
    rows <- paste("Class", 1:n.classes)
    col.classes <- structure(rainbow(n.classes), names=rows)
    col.classes[sample(length(col.classes), n.other)] <- col.other
    
    # RUN FUNCTION
    # the results need to be checked visually
    project.network.legend(col.classes, col.other) # 1, 2)
}
