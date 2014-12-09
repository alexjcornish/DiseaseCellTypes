\name{score.edges}
\alias{score.edges}
\title{
Build a context-specific interactome 
}
\description{
Integrate gene expression data with a protein-protein interaction (PPI) network to build a context-specific interactome.
}
\usage{
score.edges(g, expression, edge.attr.type=c("weight", "distance"), 
    edge.attr="score", vertex.attr="expression", 
    expression.miss=median(expression), correct.inf=TRUE)
}
\arguments{
\item{g}{\code{igraph} object. The protein-protein interaction (PPI) network. Vertex names should be saved under the vertex attribute named \code{name}.
}
\item{expression}{Named numeric vector. The percentile-normlized relative gene expression scores, possible produced using the \code{expression.transform} function.
}
\item{edge.attr.type}{Character scalar. What the edge scores should represent. If \code{weight}, then larger scores indicate stronger interactions. If \code{distance}, then smaller scores indicate stronger interactions.
}
\item{edge.attr}{Character scalar. The edge attribute under which the edge scores are saved.
}
\item{vertex.attr}{Character scalar. The vertex attribute under which the expression scores are saved.
}
\item{expression.miss}{Numeric scalar. The expression score to give to genes not present in \code{expression}.
}
\item{correct.inf}{Logical. If \code{TRUE}, then expression scores equal to zero are replaced with the smallest non-zero value (to avoid infinite distances).
}
}
\details{
Create a context-specific interactome by reweighting the edges of a network. These are the interactomes used by the \code{gene.set.compactness} function and described in the accompanying paper.

If \code{edge.attr.type == "weight"}, then edges with highly-expressed interacting partners are given larger edge scores, indicating stronger interactions. Weights are computed using:

\code{w(i,j) = x(i) * x(j)}

where \code{w(i,j)} is the weight of the edge connecting protein \code{i} and \code{j} and \code{x(i)} is the gene expression score for protein \code{i}. Pre-existing edge scores are not considered when edges are reweighted.

If \code{edge.attr.type == "distance"}, then edge scores are computed using: 

\code{d(i,j) = 1 / (x(i) * x(j))}

where \code{d(i,j)} is the distance along the edge connecting protein \code{i} and protein \code{j}.

The names of the genes in \code{expression} should correspond to the gene names in \code{g}. I an gene in \code{g} is not represented in \code{expression}, then the gene is given an expression score equal to \code{expression.miss}.
}
\value{
\code{igraph} object. The new edge scores will be saved under the edge attribute specified by \code{edge.attr}
}
\references{
Paper in preparation.
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\examples{
# create myoblast-specific interactome using FANTOM5 expression data
data(edgelist.string)
data(expression.fantom5)
g <- graph.edgelist(as.matrix(edgelist.string[, c("idA", "idB")]), directed=FALSE)
expression <- expression.transform(expression.fantom5)
g.myoblast <- score.edges(g, expression[, "myoblast"])
}
