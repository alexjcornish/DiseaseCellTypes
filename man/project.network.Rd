\name{project.network}
\alias{project.network}
\title{
Plot a network built from p-value correlations.
}
\description{
Produce and plot a correlation network from a matrix of p-values, such as the diseasome in the associated paper.
}
\usage{
project.network(p.values, col.vert=NULL, col.other="white", 
    method.cor=c("pearson", "kendall", "spearman"), rank.max=4,
    vert.size.max=10, vert.label.cex=1, vert.p.cutoff=0.1, 
    edge.width.max=10, edge.col.diff="lightgrey", \dots)
}
\arguments{
  \item{p.values}{Numeric matrix. The p-values. Each row becomes a vertex in the network. 
}
  \item{col.vert}{Named character vector. The color of each vertex in any format accepted by \code{xspline} (colors in RGB, numeric color IDs or symbolic color names). Names must be the same as the rownames of \code{p.values}. If \code{NULL}, then all vertices colored using \code{col.other}. 
}
  \item{col.other}{Character scalar. If \code{col.vert} is \code{NULL}, or the color of a vertex is \code{NA} or \code{'other'}, then the vertex is colored this color. 
}
  \item{method.cor}{Character scalar. The method used to compute the correlations between the rows.
}
  \item{rank.max}{Numeric scalar. The number of edges to add to each vertex.
}
  \item{vert.size.max}{Numeric scalar. The maximum size of the vertices.
}
  \item{vert.p.cutoff}{Numeric scaler. The p-value cutoff used identify associations for vertex sizes.
}
  \item{vert.label.cex}{Numeric scalar. The size of the vertex labels.
}
  \item{edge.width.max}{Numeric scalar. The maximum width of the edges.
}
  \item{edge.col.diff}{Character scalar. The color of edges between vertices of different colors.
}
  \item{\dots}{Additional arguments to be passed to \code{plot.igraph}.
}
}
\details{
Plot a correlation network computed from a matrix of p-values. The -log10 of p-values is first computed, then the correlation between rows computed using the method specified in \code{method.cor}. A vertex is created for each of the rows in \code{p.values}. Each vertex is connected to the \code{rank.max} vertices it correlates most strongly with. Vertices are sized by the number of significant associations they are involved in. Edge widths represent the strength of the correlation.

Vertices are colored according to \code{col.vert}. If an edge connects two vertices of different color, then it is colored \code{edge.col.diff}. If an edge connects two vertices of the same color, then it is also colored this color. 

The correlation network is returned by the function (invisibly).

This function is used to create the diseasome in the associated paper.
}
\value{
The computed correlation network in the form of an \code{igraph} object.
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\references{
Paper in preparation.
}
\seealso{
\code{\link[DiseaseCellTypes:project.network.legend]{project.network.legend}}
}
\examples{
# simulate p-values and create correlation network
n.diseases <- 5
n.cells <- 6
disease.names <- paste("disease", 1:n.diseases)
p.values <- array(runif(n.diseases * n.cells), dim=c(n.diseases, n.cells), 
    dimnames=list(disease.names, NULL))
col.vert <- structure(rainbow(n.diseases), names=disease.names)
project.network(p.values, col.vert, rank.max=2) 
}
