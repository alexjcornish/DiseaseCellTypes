\name{DiseaseCellTypes.data}
\alias{disease.classes}
\alias{disease.genes}
\alias{edgelist.string}
\alias{expression.fantom5}
\alias{pvalues.gsc}
\alias{pvalues.gso}
\alias{pvalues.text}
\alias{res.gsc}
\title{
Pre-processed data sets for the DiseaseCellTypes vignette
}
\description{
Pre-processed text-mining, gene expression, network and disease data and results from the accompanying paper.
}
\details{
\describe{
\item{disease.classes}{Classes of diseases represented within \code{disease.genes}. Classes were found by mapping each disease to a MeSH term and identifying ancestors at the second level of the MeSH ontology. When diseases map to multiple MeSH terms, the most frequently-occuring MeSH term is selected.}
\item{disease.genes}{Disease-associated genes obtained from the DisGeNET database (version 2.1). These associations have been filtered, as described in Cornish et al. (2015). DisGeNET data is available from the DisGeNET website (http://www.disgenet.org/) under the Open Database License.}
\item{edgelist.string}{Protein-protein interactions from the STRING database (version 9.1, downloaded on 2014-09-08). Only physical interactions with a confidence score greater than 0.8 have been included. Protein identifiers have been mapped to Ensembl gene identifiers. STRING data is freely available from the STRING website (http://string-db.org) under the Creative Commons Attribution 3.0 License.}
\item{expression.fantom5}{Normalized gene-wise expression counts for primary cell types from the FANTOM5 project. Samples have been normalized, grouped and combined, as described in Cornish et al. (2015). FANTOM5 data is made available from the FANTOM5 project website (http://fantom.gsc.riken.jp/) under the Creative Common Attribute 4.0 International license.}
\item{pvalues.gsc}{Disease-cell type association p-values computed using the gene set compactness (GSC) method. 10000 permutations were used.}
\item{pvalues.gso}{Disease-cell type association p-values computed using the gene set overexpression (GSO) method. 10000 permutations were used.}
\item{pvalues.text}{Disease-cell type association p-values computed using the text-mining method. Text mining was completed on 2014-08-31. On some occasions, multiple diseases and cell types are mapped to the same MeSH term, leading to identical rows and columns.}
\item{res.gsc}{Output of the \code{gene.set.compactness} function for use in vignette example.}
}
}
\references{
Our paper in preparation. 

Bauer-Mehren, A., Rautschka, M., Sanz, F. et al. (2010). \emph{DisGeNET: a Cytoscape plugin to visualize, integrate, search and analyze gene-disease networks.} Bioinformatics. 26, 2924-6.

Becker, K., Barnes, K., Bright, T. et al. (2004). \emph{The Genetic Association Database.} Nature Genetics. 36, 431-2.

Cheung, W.A, Ouellette, B.F., and Wasserman, W.W. (2012). \emph{Inferring novel gene-disease associations using Medical Subject Heading Over-representation Profiles.} Genome Medicine. 4:75.

Forrest, A., Kawaji,. H, Rehli, M. et al. (2014). \emph{A promoter-level mamalian expression atlas.} Nature. 507, 462-70.

Franceschini, A., Szklarczyk, D., Frankild, S. et al. (2013) \emph{STRING v9.1: protein-protein interaction networks, with increased coverage and integration.} Nucleic Acids Research. 41, D808-15.

Jackson, D.A., Somers, K.M. and Harvey, H.H. (1989) \emph{Similarity measures: Measures of co-occurrence and association or simply measures of co-occurrence?} The American Naturalist. 133, 436-53.

Robinson, M.D., McCarthy, D.J. and Smyth, G.K. (2010) \emph{edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.} Bioinformatics. 26, 139:40.

The UniProt Consortium (2010). \emph{The Universal Protein Resource (UniProt) in 2010.} Nucleic Acids Research. 38, D142-8.
}
\examples{
data(disease.genes)
}
