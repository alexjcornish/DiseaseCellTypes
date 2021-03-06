\name{plot.dct}
\alias{plot.dct}
\title{
Plot the results of the \code{gene.set.overexpression} or \code{gene.set.compactness} function
}
\description{
Plot the observed scores and the distribution of permuted scores output by the \code{gene.set.overexpression} or \code{gene.set.compactness} function
}
\usage{
\S3method{plot}{dct}(x, contexts=NULL, cutoff=0.05, include.pval=FALSE, \dots)
}
\arguments{
  \item{x}{Object of class \code{dct} output by the \code{gene.set.overexpression} or \code{gene.set.compactness} function.
}
  \item{contexts}{Character vector. If not \code{NULL} then these contexts are plotted. This overrides the specified \code{cutoff}.
}
  \item{cutoff}{Numeric scalar. The observed scores of context with p-values less than this are plotted.
}
  \item{include.pval}{Logical. If \code{TRUE} then the computed p-value is plotted next to each observed score.
}
  \item{\dots}{Additional arguments to be passed to \code{plot}.      
}
}
\details{
Plot the output of the \code{gene.set.overexpression} or \code{gene.set.compactness} function. If more than 1 permutation is completed then a histogram showing the distribution of permutation scores is shown in grey. Observed scores for contexts are shown as red lines. 

If \code{contexts} is \code{NULL} then the observed scores of contexts with p-values less than \code{cutoff} are shown. If \code{contexts} is not \code{NULL} then the observed scores of contexts specified by \code{contexts} are shown.
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\references{
Cornish, A.J., Filippis, I., David, A. and Sternberg, M.J.E. (2015) \emph{Exploring the cellular basis of human disease through a large-scale mapping of deleterious genes to cell types.} Genome Medicine. 7:95.
}
\seealso{
\code{\link[DiseaseCellTypes:gene.set.compactness]{gene.set.compactness}},
\code{\link[DiseaseCellTypes:gene.set.overexpression]{gene.set.overexpression}}
}
\examples{
# run the gene.set.overexpression function on toy data
n.genes <- 50
n.contexts <- 20
gene.names <- paste("gene", 1:n.genes)
context.names <- paste("context", 1:n.contexts)
expression <- array(runif(n.genes * n.contexts), dim=c(n.genes, n.contexts), 
    dimnames=list(gene.names, context.names))
genes <- sample(gene.names, 4)
res <- gene.set.overexpression(expression, genes, n.perm=1000)

# plot the results
plot(res, cutoff=0.2)
}
