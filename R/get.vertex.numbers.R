get.vertex.numbers <- function(
    g, 
    v
) {
    # given a graph g and a set of vertex IDs v, return the number of the vertices in g present in v
    # v can be a igraph vertex (igraph.vs) object, a character vector or a numeric vector
    # igraph.vs objects are of type numeric
    # v can only be a character vector if vertex names are stored under the vertex attribute "name"
    # if not all v in g, then an error is produced
    
    if (class(g) != "igraph") stop("g is not an igraph object")
    if (is.null(v)) return(NULL)
    if (!is.numeric(v) & !is.character(v)) stop("v is not an igraph object, character or numeric vector")
    
    # v is igraph.vs object
    if (class(v) == "igraph.vs") {
        if (!all(v %in% V(g))) stop("not all v in g")
        v.num <- match(v, V(g))
    }
    
    # v is character vector of vertex names
    if (is.character(v)) {
        if (is.null(V(g)$name)) stop("g does not have any vertex names")
        if (!all(v %in% V(g)$name)) stop("not all v in g")
        v.num <- match(v, V(g)$name)
    }
    
    # v is numeric vector of vertex numbers
    if (class(v) != "igraph.vs" & is.numeric(v)) {
        if (max(v) > vcount(g)) stop("not all v in g")   
        v.num <- v
    }
    
    as.integer(v.num)
}
