\name{gene.set.overexpression}
\alias{gene.set.overexpression}
\title{
Run the gene set overexpression method
}
\description{
Identify disease-cell-type associations using the gene set overexpression (GSO) method.
}
\usage{
gene.set.overexpression(expression, genes, n.perm=10000, expression.miss=median(expression))
}
\arguments{
  \item{expression}{Numeric matrix. Matrix of percentile-normalised relative gene expression scores. Rows should represent genes and be named with the gene name. Columns should represent contexts (i.e. cell types) and be named with the context name. \code{expression} should contain at least 2 columns and can be produced using the \code{expression.transform} function. 
}
  \item{genes}{Character vector. The names of the disease-associated genes.
}
  \item{n.perm}{Numeric scalar. Number of permutations to complete.
}
  \item{expression.miss}{Numeric scalar. The expression score to give to genes not present in \code{expression}.
}
}
\details{
Use a permutation-based appraoch to identify contexts (i.e. cell types) in which sets of genes are overexpressed. Permuted expression profiles are created by randomly redistributing the expression scores of each gene between contexts. The mean expression score of the genes in \code{genes} is computed for each observed and permuted expression profile. An empirical p-value is produced describing the probability that a set of disease-associated genes is as overexpressed as obsered in a profile by chance is produced by counting the number of permuted profiles in which the mean expression score of the genes is greater than the mean expression score of the genes in the observed profile.

Expression data should represent percentile-normalised relative expression scores, possibly created using the \code{expression.transform} function. Contexts can represent are range of biological entities, including cell types, tissues and disease states. This method is refered to as the `gene set overexpression' method in Cornish et al. (2015).
}
\value{
A object of class \code{dct}. This object is a list and contains the following elements
\item{obs}{The observed mean expression score for the set of genes in \code{genes} in each of the contexts in \code{expression}.
}
\item{perm}{The permuted mean expression scores. \code{NULL} if no permutations are completed.
}
\item{pval}{The overexpression-based p-value for each of the contexts. \code{NA} if no permutations are completed.
}
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\references{
Cornish, A.J., Filippis, I., David, A. and Sternberg, M.J.E. (2015) \emph{Exploring the cellular basis of human disease through a large-scale mapping of deleterious genes to cell types.} Genome Medicine. 7:95.
}
\seealso{
\code{\link[DiseaseCellTypes:expression.transform]{expression.transform}},
\code{\link[DiseaseCellTypes:plot.dct]{plot.dct}}
}
\examples{
# create data and run expression sig
n.genes <- 50
n.contexts <- 5
gene.names <- paste("gene", 1:n.genes)
context.names <- paste("context", 1:n.contexts)
expression <- array(runif(n.genes * n.contexts) * 10, dim=c(n.genes, n.contexts), 
    dimnames=list(gene.names, context.names))
expression <- expression.transform(expression)
genes <- sample(gene.names, 4)
res <- gene.set.overexpression(expression, genes, n.perm=20)
res
}
