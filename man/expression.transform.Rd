\name{expression.transform}
\alias{expression.transform}
\title{
Transform raw gene expression values
}
\description{
Transform a matrix of raw gene expression values to a matrix of percentile-normalized relative gene expression scores.
}
\usage{
expression.transform(x, remove.zero=TRUE)
}
\arguments{
\item{x}{Numeric matrix. Raw expression values. Each row should represent a gene and each column a context (i.e. a cell type).
}
\item{remove.zero}{Logical. If \code{TRUE} rows in which all values equal zero are removed.
}
}
\details{
Convert raw gene-wise expression values to percentile-normalised relative expression scores for use within the \code{gene.set.compactness} and \code{gene.set.overexpression} functions. First, each gene-wise expression value is divided by the mean expression values across all contexts to produce relative values. These relative values are then percentile-normalized so that the expression scores for each context (each column in \code{x}) range uniformly between 0 and 1. A score of 1 indicates that the gene is the most relatively-overexpressed gene in the context, while a score of 0 indicates that it is the most relatively-underexpressed.
}
\value{
Numeric matrix. Percentile-normalized relative expression scores. 
}
\references{
Cornish, A.J., Filippis, I., David, A. and Sternberg, M.J.E. (2015) \emph{Exploring the cellular basis of human disease through a large-scale mapping of deleterious genes to cell types.} Genome Medicine. 7:95.
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\examples{
data(expression.fantom5)
expression.fantom5[1:3, 1:3]

# transform gene expression using the mean method
expression.rel <- expression.transform(expression.fantom5)
expression.rel[1:3, 1:3]
}
