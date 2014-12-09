\name{distance.rwr}
\alias{distance.rwr}
\title{
Compute graph distance using random walks with restart (RWR)
}
\description{
Use the random walk with restart (RWR) method to compute the distance between vertex pairs on a graph.
}
\usage{
distance.rwr(g, v=V(g), edge.attr=NULL, rwr.r=0.7, rwr.cutoff=1e-5, correct.inf=TRUE)
}
\arguments{
\item{g}{\code{igraph} object. The graph on which to work.
}
\item{v}{\code{igraph} object or numeric vector. The vertices from which each distance is computed.
}
\item{edge.attr}{Character scalar. The name of the edge attribute under which edge weights are found. If \code{NULL}, then all edges are assumed to have equal weight.
}
\item{rwr.r}{Numeric scalar. The restart probability.
}
\item{rwr.cutoff}{Numeric scalar. Value used to compute the change cutoff that controls when iterations are terminated.
}
\item{correct.inf}{Logical. If \code{TRUE}, then infinite distances produced by the RWR method are converted to the largest finite distance.
}
}
\details{
Use the iterative random walk with restart (RWR) method described by Kohler et al. to compute the distances between pairs of vertices in a graph. This method is used by the \code{gene.set.compactness} function to compute the distances between the vertices in the vertex set. Distances are computed between each vertex in \code{v} and every vertex in \code{g}. 

Let \code{A} be a column-normalized adjacency matrix for \code{g} using the edge weights specified in \code{edge.attr}. Larger edges weights should represent stronger interactions. \code{p} is a probability matrix with dimensions equal to the number of vertices in \code{g}. The element \code{p^{t}(i,j)} is the probability that at time \code{t}, a random walker starting from vertex \code{i} is located at vertex \code{j}. \code{p^0} is the initial probability distribution and an identity matrix. \code{r} is the restart probability specified by \code{rwr.r}. The final distances are computed iteratively using:

\code{p^(t+1) = (1-r)Ap^t + rp^0}

To save computational time, distances are only computed between vertices in \code{v} and vertices in \code{g}. Iterations are stopped when the distance (Manhattan) between \code{p^(t-1)} and \code{p^t} falls below \code{rwr.cutoff} multipled by the number of vertices in \code{v}.

Distances between vertex pairs are computed by taking the reciprocal of the final probabilities.
}
\value{
Numeric matrix. The distances between the vertices in \code{v} and the vertices in \code{g}.
}
\references{
Kohler, S., Bauer, S., Horn, D. et al. (2008) \emph{Walking the interactome for prioritization of candidate disease genes.} The American Journal of Human Genetics. 82, 949-58.
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\examples{
# create a graph and compute the vertex pair distances
g <- barabasi.game(6, directed=FALSE)
distance.rwr(g)
plot(g, layout=layout.fruchterman.reingold)
}
