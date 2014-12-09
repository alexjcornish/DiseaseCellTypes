\name{gene.set.overexpression}
\alias{gene.set.overexpression}
\title{
Compute gene set overexpression significance
}
\description{
Identify associations between sets of disease-associated genes and context (i.e. cell types) using the gene set overexpression method.
}
\usage{
gene.set.overexpression(expression, genes, n.perm=10000, expression.miss=median(expression))
}
\arguments{
  \item{expression}{Numeric matrix. Matrix of percentile-transformed relative expression scores. Rows should represent genes and be named with the gene name. Columns should represent contexts (i.e. cell types) and be named with the context name. \code{expression} should contain at least 2 columns and can be produced using the \code{expression.transform} function. 
}
  \item{genes}{Character vector. The names of the genes in the gene set to test. 
}
  \item{n.perm}{Numeric scalar. Number of permutations to complete.
}
  \item{expression.miss}{Numeric scalar. The expression score to give to genes not present in \code{expression}.
}
}
\details{
Use a permutation-based appraoch to identify contexts (i.e. cell types) in which sets of genes are overexpressed. Permuted expression profiles are created by randomly redistributing the expression scores of each gene between contexts. The mean expression score of the genes in \code{genes} is computed for each observed expression profile and the permuted expression profiles. An empirical p-value, based on the number of permuted mean scores greater than each observed mean scores, is produced for each context. 

Expression data should represent percentile-normalised relative expression scores, possibly created using the \code{expression.transform} function. Contexts can represent are range of biological entities, including cell types, tissues and disease states. This method is refered to as the `gene set overexpression' method in the associated paper.
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
Paper in preparation
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
