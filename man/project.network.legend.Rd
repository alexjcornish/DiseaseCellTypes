\name{project.network.legend}
\alias{project.network.legend}
\title{
Plot legend for the \code{project.network} function.
}
\description{
Plot a legend of classes in a separate plot. To be used with the \code{project.network} function.
}
\usage{
project.network.legend(col.classes, col.other="white", pt.cex=3, \dots)
}
\arguments{
  \item{col.classes}{Named character vector. The color of each class, in any format accepted by \code{xspline} (colors in RGB, numeric color IDs or symbolic color names). Vector should be named with the name of each class.
}
  \item{col.other}{Character scalar. Classes of this color are combined under an 'Other' class. 
}
  \item{pt.cex}{Numeric scalar. The expansion factor for the points.
}
  \item{\dots}{Additional arguments to be passed to \code{legend}.
}
}
\details{
Create a legend of classes. Classes of the same color can be combined under an 'Other' class if they all have the same color as \code{col.other}.
}
\author{
Alex J. Cornish \email{a.cornish12@imperial.ac.uk}
}
\references{
Paper in preparation.
}
\seealso{
\code{\link[DiseaseCellTypes:project.network]{project.network}}
}
\examples{
# create classes with 'Others' and plot legend
n.classes <- 5
n.other <- 2
col.other <- "white"
class.names <- paste("Disease", 1:n.classes)
col.classes <- structure(rainbow(n.classes), names=class.names)
col.classes[sample(length(col.classes), n.other)] <- col.other
project.network.legend(col.classes, col.other)
}
