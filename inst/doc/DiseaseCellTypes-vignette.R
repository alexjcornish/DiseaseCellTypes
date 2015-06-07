## ----setup, echo=FALSE---------------------------------------------------
options(width = 75)
options(useFancyQuotes=FALSE)

## ----hook-printfun, echo=FALSE-------------------------------------------
library(knitr)
library(formatR)
knit_hooks$set(printfun = function(before, options, envir) {
    if (before) return()
    txt = capture.output(dump(options$printfun, '', envir = envir))
    ## reformat if tidy=TRUE
    if (options$tidy)
        txt = tidy.source(text=txt, output=FALSE, width.cutoff=30L, keep.comment=TRUE, keep.blank.line=FALSE)$text.tidy
    paste(c('\n```r\n', txt, '\n```\n'), collapse="\n")
})

## ----create_network------------------------------------------------------
# for reproducibility
set.seed(1)

# load DiseaseCellTypes package
require(DiseaseCellTypes)

# load PPI data
data(edgelist.string)
edgelist.string[1:3, ]

# create igraph object
g <- graph.edgelist(as.matrix(edgelist.string[, c("ID.A", "ID.B")]), directed=FALSE)

## ----transform_expression------------------------------------------------
# load FANTOM5 gene expression data
data(expression.fantom5)
expression.fantom5[1:3, 1:2]

# transform to percentile-normlized relative gene expression scores
expression <- expression.transform(expression.fantom5)
expression[1:3, 1:2]

## ----disease_genes-------------------------------------------------------
# load disease-gene associations
data(disease.genes)
disease.genes[["uveitis, anterior"]]

## ----gsc-----------------------------------------------------------------
# identify cell types associated with panic disorder-
disease <- "bipolar disorder"
genes <- disease.genes[[disease]]

# remove genes not contained within the network
genes <- genes[genes %in% V(g)$name]

# run the gene set compactness function
#res.gsc <- gene.set.compactness(expression, genes, g, n.perm=100)

## ----gsc_load------------------------------------------------------------
# load the previously-computed results
data(res.gsc)

## ----gsc_investigate, results="asis"-------------------------------------
# look at the p-value of the neuron
res.gsc$pval[["neuron"]]

# plot the results to compare permuted and observed scores
plot(res.gsc)

## ----gso-----------------------------------------------------------------
# run the GSO method
res.gso <- gene.set.overexpression(expression, genes, n.perm=100)

## ----gso_investigate, results="asis"-------------------------------------
# look at the p-value of the neuron
res.gso$pval[["neuron"]]

# plot the results to compare permuted and observed scores
plot(res.gso)

## ----create_neuron_interactome-------------------------------------------
# create global network
data(edgelist.string)
g <- graph.edgelist(as.matrix(edgelist.string[, c("ID.A", "ID.B")]), directed=FALSE)

# load and transform expression data
data(expression.fantom5)
expression <- expression.transform(expression.fantom5)

# create neuron-specific interactome
g.neuron <- score.edges(g, expression[, "neuron"])

## ----monocyte_interactome_attributes-------------------------------------
# edge scores
E(g.neuron)$score[1:5]

# vertex expression scores
V(g.neuron)$expression[1:5]

## ----recreate_gsc_gso----------------------------------------------------
# # parameters 
# min.n.genes <- 6 # minumum number of genes for a disease to be tested
# n.perm <- 10000 # number of permutations
# rwr.r <- 0.7 # the RWR restart probability  
# rwr.cutoff <- 1e-5 # the RWR iteration termination cutoff
# parallel <- 4 # the number of cores to use 
# 
# # load and create global network
# data(edgelist.string)
# g <- graph.edgelist(as.matrix(edgelist.string[, c("ID.A", "ID.B")]), directed=FALSE)
# genes.g <- V(g)$name
# 
# # load and format gene expression data
# data(expression.fantom5)
# expression <- expression.transform(expression.fantom5)
# genes.expression <- rownames(expression)
#  
# # load disease-associated genes
# data(disease.genes)
# disease.genes <- sapply(disease.genes, function(genes) genes[genes %in% genes.g], simplify=F) # remove disease-associated genes not found in network or expression
#
# # select diseases with at least min.n.genes associated genes
# disease.genes <- disease.genes[sapply(disease.genes, length) >= min.n.genes]
# diseases <- names(disease.genes)
# 
# # for each disease, apply the GSC method
# pvalues.gsc <- t(sapply(diseases, function(disease) gene.set.compactness(expression, disease.genes[[disease]], g, n.perm, rwr.r=rwr.r, rwr.cutoff=rwr.cutoff, parallel=parallel)$pval))
# 
# # for each disease, apply the GSO method
# pvalues.gso <- t(sapply(diseases, function(disease) gene.set.overexpression(expression, disease.genes[[disease]], n.perm)$pval))

## ----setup_heatmap, results="asis"---------------------------------------
# load the associations
data(pvalues.gsc)

# compute the -log10 of the p-values
pvalues.ml10 <- -log10(pvalues.gsc)

# plot heatmap
require(gplots)
n.cols <- 500
cols <- colorRampPalette(c("#D3DDDC", "#02401B"))(n.cols)
col.breaks <- c(0, seq(-log10(0.05), max(pvalues.ml10), length.out=n.cols))
heatmap.2(pvalues.ml10[1:10, 1:10], Rowv=TRUE, Colv=TRUE, dendrogram="both", breaks=col.breaks, col=cols, notecol="black", trace="none", margins=c(18,15), key=TRUE, density.info="none")

## ----load_pvalues--------------------------------------------------------
# load the associations
data(pvalues.gsc)
data(pvalues.gso)
data(pvalues.text)

# keep associations for diseases with 6 or more associated genes
data(disease.genes)
diseases <- names(disease.genes)[sapply(disease.genes, length) >= 6]
pvalues.gsc <- pvalues.gsc[diseases, ]
pvalues.gso <- pvalues.gso[diseases, ]
pvalues.text <- pvalues.text[diseases, ]

## ----adjust_pvalues------------------------------------------------------
# adjust p-values for multiple testing
pvalues.adj.gsc <- array(p.adjust(as.vector(pvalues.gsc), method="BH"), dim=dim(pvalues.gsc), dimnames=dimnames(pvalues.gsc))
pvalues.adj.gso <- array(p.adjust(as.vector(pvalues.gso), method="BH"), dim=dim(pvalues.gso), dimnames=dimnames(pvalues.gso))
pvalues.adj.text <- array(p.adjust(as.vector(pvalues.text), method="BH"), dim=dim(pvalues.text), dimnames=dimnames(pvalues.text))

## ----compare_associations------------------------------------------------
# count the number of associations supported by text
res <- matrix(NA, 2,2, dimnames=list(c("supported", "not supported"), c("GSO", "GSC")))
cutoff <- 0.10
res["supported", "GSO"] <- sum(pvalues.adj.gso <= cutoff & pvalues.adj.text <= cutoff)
res["supported", "GSC"] <- sum(pvalues.adj.gsc <= cutoff & pvalues.adj.text <= cutoff)
res["not supported", "GSO"] <- sum(pvalues.adj.gso <= cutoff & !pvalues.adj.text <= cutoff)
res["not supported", "GSC"] <- sum(pvalues.adj.gsc <= cutoff & !pvalues.adj.text <= cutoff)
res

## ----diseasome, results="asis"-------------------------------------------
# load the disease classes
data(disease.classes)

# diseases to set colours for
disease.cols <- c("Cardiovascular Diseases", "Digestive System Diseases", "Immune System Diseases", "Mental Disorders", "Nervous System Diseases", "Nutritional and Metabolic Diseases", "Respiratory Tract Diseases", "Skin and Connective Tissue Diseases", "Urogenital Diseases and Pregnancy Complications")
col.other <- "white"

# set colour for each disease class
classes <- sort(unique(disease.classes))
cols <- rep(col.other, length(classes))
names(cols) <- classes
for (i in 1:length(disease.cols)) cols[disease.cols[i]] <- rainbow(length(disease.cols))[i]

# identify the class color for each disease
disease.cols <- cols[disease.classes]
names(disease.cols) <- names(disease.classes)

# create the diseasome
project.network(pvalues.adj.gsc[names(disease.cols), ], col.vert=disease.cols, col.other=col.other, vert.size.max=5, vert.label.cex=0.6, edge.width.max=5, layout=layout.fruchterman.reingold) 
project.network.legend(cols, col.other)

## ----create_monocyte_interactome-----------------------------------------
# create global network
data(edgelist.string)
g <- graph.edgelist(as.matrix(edgelist.string[, c("ID.A", "ID.B")]), directed=FALSE)

# load expression data
data(expression.fantom5)
expression <- expression.transform(expression.fantom5)

# create monocyte-specific interactome
g.monocyte <- score.edges(g, expression[, "monocyte"])

## ----identify_disease_genes----------------------------------------------
# identify all psoriasis-associated genes
data(disease.genes)
psoriasis.genes <- disease.genes[["psoriasis"]]

# remove genes not present in the network
psoriasis.genes <- psoriasis.genes[psoriasis.genes %in% V(g.monocyte)$name]

# remove genes with more than 15 interacting partners
psoriasis.genes <- psoriasis.genes[degree(g.monocyte)[psoriasis.genes] <= 15]

## ----psoriasis_subgraph--------------------------------------------------
# produce psoriasis subgraph
g.subgraph <- disease.subgraph(g.monocyte, psoriasis.genes, n.bins=500, edge.width.max=10)
plot(g.subgraph)

## ----session_info--------------------------------------------------------
sessionInfo()

