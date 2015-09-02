\name{disease.subgraph}
\alias{disease.subgraph}
\title{
Extract and plot a disease subgraph.
}
\description{
Plot a subgraph containing a set of disease-associated genes and their interacting partners.
}
\usage{
disease.subgraph(g, genes, edge.attr="score", vert.attr="expression", 
    filename.plot=NULL, filename.legend.vert=NULL, filename.legend.edge=NULL, 
    n.bins=500, vert.size.disease=6, vert.size.not.disease=3, edge.width.max=20,
    vert.col.high="black", vert.col.low="lightgrey", edge.col.high="blue",
    edge.col.low="lightgrey", \dots)
}
\arguments{
\item{g}{\code{igraph} object. The network from which the subgraph is extracted.
}
\item{genes}{Character vector. The names of the disease-associated genes.
}
\item{edge.attr}{Character scalar. The edge attribute under which the edge scores are saved.
}
\item{vert.attr}{Character scalar. The vertex attribute under which the expression scores are saved.
}
\item{filename.plot}{Character scalar. The name of the file within which the subgraph is plotted. If \code{NULL} the subgraph is not plotted.
}
\item{filename.legend.vert}{Character scalar. The name of the file within which the vertex color legend is plotted. If \code{NULL} the legend is not plotted.
}
\item{filename.legend.edge}{Character scalar. The name of the file within which the edge color legend is plotted. If \code{NULL} the legend is not plotted.
}
\item{n.bins}{Numeric scalar. The number of bins the colors are split into.
}
\item{vert.size.disease}{Numeric scaler. The size of the disease gene vertices. 
}
\item{vert.size.not.disease}{Numeric scaler. The size of the non-disease gene vertices. 
}
\item{edge.width.max}{Numeric scaler. The maximum width of the edges. 
}
\item{vert.col.high}{Named character vector. The color of high-scoring vertices in any format accepted by \code{xspline} (colors in RGB, numeric color IDs or symbolic color names). Names must be the same as the rownames of \code{p.values}. If \code{NULL} then all vertices colored using \code{col.other}.  
}
\item{vert.col.low}{Named character vector. The color of low-scoring vertices.
}
\item{edge.col.high}{Named character vector. The color of high-scoring edges.
}
\item{edge.col.low}{Named character vector. The color of low-scoring edges.
}
\item{\dots}{Additional arguments to be passed to \code{pdf}. 
}
}
\details{
Extract a subgraph from a context-specific interactome containing a set of disease genes and their interacting partners. The context-specific interactome \code{g} should be created using the \code{score.edges} function.

Disease genes are represented as squares and non-disease genes as circles. Vertices are colored according to the expression of the genes, along a scale from \code{vert.col.low} for the lowest-weight vertices to \code{vert.col.high} for the highest. Edges are colored by their score, along a scale from \code{edge.col.low} for the lowest-weight edges to \code{edge.col.high} for the highest. Edges are also weighted by their score. 
}
\value{
The subgraph in the form of an \code{igraph} object.
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\seealso{
\code{\link[DiseaseCellTypes:score.edges]{score.edges}}
}
\examples{
# create cell type-specific interactome
data(edgelist.string)
data(expression.fantom5)
g <- graph.edgelist(as.matrix(edgelist.string[, c("ID.A", "ID.B")]), directed=FALSE)
expression <- expression.transform(expression.fantom5)
g.myoblast <- score.edges(g, expression[, "myoblast"])

# simulate disease genes
genes.disease <- sample(V(g.myoblast)$name, 5) 

# produce cell type-specific interactome
subgraph <- disease.subgraph(g.myoblast, genes.disease)
plot(subgraph)
}
