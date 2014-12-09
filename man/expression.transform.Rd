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
\item{x}{Numeric matrix. Absolute expression values. Each row should represent a gene and each column a context (i.e. a cell type).
}
\item{remove.zero}{Logical. If \code{TRUE}, rows in which all values equal zero are removed.
}
}
\details{
Convert raw gene-wise expression values to relative expression scores for use within the \code{gene.set.compactness} and \code{gene.set.overexpression}. First, each gene-wise expression value is divided by the mean expression values across all contexts to produce relative values. These relative values are then percentile-normalized, so that the expression scores of each profile (each column in \code{x}) range uniformly between 0 and 1. A score of 1 indicates that the gene is the most relatively-overexpressed gene in a profile, while a score of 0 indicates that it is the most relatively-underexpressed.
}
\value{
Numeric matrix. Percentile-normalized relative expression scores. 
}
\references{
Paper in preparation.
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
