\name{gene.set.compactness}
\alias{gene.set.compactness}
\title{
Run the gene set compactness method
}
\description{
Identify disease-cell-type associations using the gene set compactness (GSC) method.
}
\usage{
gene.set.compactness(expression, genes, g, n.perm=10000, expression.miss=median(expression), 
    rwr.r=0.7, rwr.cutoff=1e-5, parallel=NULL, verbose=TRUE)
}
\arguments{
  \item{expression}{Numeric matrix. Matrix of percentile-normalised relative gene expression scores. Rows should represent genes and be named with the gene name. Columns should represent contexts (i.e. cell types) and be named with the context name. \code{expression} should contain at least 2 columns and can be produced using the \code{expression.transform} function. 
} 
  \item{genes}{Character vector. The names of the disease-associated genes.
}
  \item{g}{\code{igraph} object. The PPI network used to created the context-specific interactomes. Gene names should be stored under a vertex attribute named \code{name}.
}
  \item{n.perm}{Numeric scalar. Number of permutations to complete.
}
  \item{expression.miss}{Numeric scalar. The expression score given to genes not present in \code{expression}.
}
  \item{rwr.r}{Numeric scalar. The RWR restart probability.
}
  \item{rwr.cutoff}{Numeric scalar. The value used to determine when the iterative algorithm used to compute the RWR distances terminates.
}
  \item{parallel}{Numeric scalar or \code{NULL}. If a numeric scalar and parallel computing is available then this value determines the number of cores that the computation will be split over. See the \code{snow} package for further details.
}
  \item{verbose}{Logical. If \code{TRUE} then messages about the progress of the function are displayed. 
}
}
\details{
Use a permutation-based approach to identify context-specific interactomes in which a set of genes are significantly more compact than expected by chance. The compactness score is defined by Glaab et al. as the mean distance between pairs of vertices in a set of vertices. Here, we use the random walk with restart (RWR) method to compute the distances between vertex pairs, as described by Kohler et al. \code{rwr.r} determines the probability that a random walker will return to the start node, while \code{rwr.cutoff} determines when the iterative process will terminate (as described for the \code{distance.rwr} function).

Context-specific interactomes are generated within the function by integrating expression and protein-protein interaction (PPI) data (\code{g}). Edges are re-scored using the product of the expression scores of the interactors. Expression data should represent percentile-normalised relative gene expression scores, created using the \code{expression.transform} function. The context-specific interactomes are not output by this function, but can be created separately using the \code{score.edges} function. 

Observed interactomes are created for each context in \code{expression}. \code{n.perm} permuted interactomes are created by permuting the expression scores, then using these permuted profiles to create permuted interactomes. Empirical p-values, describing whether the gene set is significantly more compact on any of the observed interactomes than expected by chance given the PPI data used to created the interactomes, are produced for each context by counting the number of permuted interactomes in which the compactness score is smaller than the compactness score in each observed interactome.

Contexts can represent a range of biological contexts, including cell types, tissues and disease states. This method is refered to as the `gene set compactness' method in Cornish et al. (2015).
}
\value{
A object of class \code{dct}. This object is a list and contains the following elements
\item{obs}{The compactness score of the set of genes in \code{genes} in each of the contexts in \code{expression}.
}
\item{perm}{The compactness score of the set of genes in \code{genes} in each of the permuted interactomes. \code{NULL} if no permutations are completed.
}
\item{pval}{The compactness-based p-value for each of the contexts. \code{NA} if no permutations are completed.
}
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\references{
Cornish, A.J., Filippis, I., David, A. and Sternberg, M.J.E. (2015) \emph{Exploring the cellular basis of human disease through a large-scale mapping of deleterious genes to cell types.} Genome Medicine. 7:95.

Cornish, A.J. and Markowetz, F. (2014) \emph{SANTA: Quantifying the Functional Content of Molecular Networks.}. PLOS Computational Biology. 10:9, e1003808.

Glaab, E., Baudot A., Krasnogor N. and Valencia A. (2010) \emph{Extending pathways and processes using molecular interaction networks to analyse cancer genome data.} BMC Bioinformatics. 11(1): 597:607.

Kohler, S., Bauer, S., Horn, D. et al. (2008) \emph{Walking the interactome for prioritization of candidate disease genes.} The American Journal of Human Genetics. 82, 949-958.
}
\seealso{
\code{\link[DiseaseCellTypes:expression.transform]{expression.transform}},
\code{\link[DiseaseCellTypes:plot.dct]{plot.dct}},
\code{\link[DiseaseCellTypes:score.edges]{score.edges}}
}
\examples{
# parameters used in the example
n.genes <- 50
n.contexts <- 5
gene.names <- paste("gene", 1:n.genes)
context.names <- paste("context", 1:n.contexts)

# create expression data
expression <- array(runif(n.genes * n.contexts) * 10, dim=c(n.genes, n.contexts), 
    dimnames=list(gene.names, context.names))
expression <- expression.transform(expression)

# create the gene set
genes <- sample(gene.names, 4)

# create the PPI data
g <- barabasi.game(n.genes, directed=FALSE)
V(g)$name <- gene.names

# run the function
res <- gene.set.compactness(expression, genes, g, n.perm=20)
res
}
