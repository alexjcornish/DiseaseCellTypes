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
        txt = tidy.source(text=txt, output=FALSE,
                          width.cutoff=30L, keep.comment=TRUE,
                          keep.blank.line=FALSE)$text.tidy
    paste(c('\n```r\n', txt, '\n```\n'), collapse="\n")
    })

## ----create_network------------------------------------------------------
# load DiseaseCellTypes package
set.seed(999) # for reproducibility
require(DiseaseCellTypes)

# load STRING PPI data
data(edgelist.string)
edgelist.string[1:2, ]

# create igraph object
g <- graph.edgelist(as.matrix(edgelist.string[, c("idA", "idB")]), directed=FALSE)

## ----transform_expression------------------------------------------------
# load FANTOM5 gene expression data
data(expression.fantom5)
expression.fantom5[1:3, 1:2]

# transform to percentile-normlized relative expression scores
expression <- expression.transform(expression.fantom5)
expression[1:3, 1:2]

## ----disease_genes-------------------------------------------------------
# load disease-gene associations
data(disease.genes)
disease.genes[["stomach ulcer"]]

## ----gsc-----------------------------------------------------------------
# identify cell types associated with panic disorder-
disease <- "panic disorder"
genes <- disease.genes[[disease]]

# remove genes not contained within the network
genes <- genes[genes %in% V(g)$name]

# run the gene set compactness function
#res.gsc <- gene.set.compactness(expression, genes, g, n.perm=100)

## ----gsc_load------------------------------------------------------------
# load the results of the gene.set.compactness function
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
g <- graph.edgelist(as.matrix(edgelist.string[, c("idA", "idB")]), directed=FALSE)

# load and transform expression data
data(expression.fantom5)
expression <- expression.transform(expression.fantom5)

# create neuron-specific interactome
g.neuron <- score.edges(g, expression.fantom5[, "neuron"])

## ----monocyte_interactome_attributes-------------------------------------
# edge scores
E(g.neuron)$score[1:5]

# vertex expression scores
V(g.neuron)$expression[1:5]

## ----recreate_gsc_gso----------------------------------------------------
# # parameters 
# n.diseases <- 100 # number of diseases
# n.perm <- 10000 # number of permutations
# rwr.r <- 0.7 # the RWR restart probability  
# rwr.cutoff <- 1e-5 # the RWR iteration termination cutoff
# parallel <- 10 # the number of cores to use 
# 
# # load and create global network
# data(edgelist.string)
# g <- graph.edgelist(as.matrix(edgelist.string[, c("idA", "idB")]), directed=FALSE)
# genes.g <- V(g)$name
# 
# # load and format gene expression data
# data(expression.fantom5)
# expression <- expression.transform(expression.fantom5)
# expression <- expression[rownames(expression) %in% genes.g, ]
# genes.expression <- rownames(expression)
#  
# # load disease-associated genes
# data(disease.genes)
# 
# # remove disease-associated genes not found in network or expression
# disease.genes <- sapply(disease.genes, function(genes) genes[genes %in% genes.g], simplify=F)
# disease.genes <- sapply(disease.genes, function(genes) genes[genes %in% genes.expression], simplify=F)
# 
# # remove diseases not identified as diseases in MeSH
# diseases.to.remove <- c("dna damage", "oxidative stress", "impaired cognition", "retinal vascular occlusion")
# disease.genes <- disease.genes[!names(disease.genes) %in% diseases.to.remove]
# 
# # select the 100 diseases with the most disease-associated genes for use
# disease.genes <- disease.genes[order(sapply(disease.genes, length), decreasing=T)][1:n.diseases]
# 
# # for each disease, apply the GSC method
# res.full.gsc <- list()
# for (disease in names(disease.genes)) {
#     res.full.gsc[[disease]] <- gene.set.compactness(expression, disease.genes[[disease]], g, n.perm, rwr.r=rwr.r, rwr.cutoff=rwr.cutoff, parallel=parallel)
# }
# 
# # for each disease, apply the GSO method
# res.full.gso <- list()
# for (disease in names(disease.genes)) {
#     res.full.gso[[disease]] <- gene.set.overexpression(expression, disease.genes[[disease]], n.perm)
# }
# 
# # convert to patrix of p-values
# pvalues.gsc <- t(sapply(res.full.gsc, function(disease) disease$pval))
# pvalues.gso <- t(sapply(res.full.gso, function(disease) disease$pval))

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

# ensure that the rows and columns of the matrices are ordered the same
diseases <- rownames(pvalues.gsc)
cell.types <- colnames(pvalues.gsc)
pvalues.gso <- pvalues.gso[diseases, cell.types]
pvalues.text <- pvalues.text[diseases, cell.types]

## ----adjust_pvalues------------------------------------------------------
# adjust p-values for multiple testing
pvalues.adj.gsc <- array(p.adjust(as.vector(pvalues.gsc), method="BH"), dim=dim(pvalues.gsc), dimnames=dimnames(pvalues.gsc))
pvalues.adj.gso <- array(p.adjust(as.vector(pvalues.gso), method="BH"), dim=dim(pvalues.gso), dimnames=dimnames(pvalues.gso))
pvalues.adj.text <- array(p.adjust(as.vector(pvalues.text), method="BH"), dim=dim(pvalues.text), dimnames=dimnames(pvalues.text))

## ----compare_associations------------------------------------------------
# count the number of associations supported by text
res <- matrix(NA, 2,2, dimnames=list(c("supported by text", "not supported by text"), c("GSO", "GSC")))
cutoff <- 0.10
res["supported by text", "GSO"] <- sum(pvalues.adj.gso < cutoff & pvalues.adj.text < cutoff)
res["supported by text", "GSC"] <- sum(pvalues.adj.gsc < cutoff & pvalues.adj.text < cutoff)
res["not supported by text", "GSO"] <- sum(pvalues.adj.gso < cutoff & !pvalues.adj.text < cutoff)
res["not supported by text", "GSC"] <- sum(pvalues.adj.gsc < cutoff & !pvalues.adj.text < cutoff)
res

## ----diseasome, results="asis"-------------------------------------------
# load the disease classes
data(disease.classes)
col.other <- "white"

# set colours for each of the classes
classes <- sort(unique(disease.classes))
classes.cols <- structure(rep(col.other, length(classes)), names=classes)
classes.cols["Cardiovascular Diseases"] <- "red"
classes.cols["Digestive System Diseases"] <- "orange"
classes.cols["Immune System Diseases"] <- "yellow"
classes.cols["Mental Disorders"] <- "lightgreen"
classes.cols["Musculoskeletal Diseases"] <- "darkgreen"
classes.cols["Nutritional and Metabolic Diseases"] <- "blue"
classes.cols["Skin and Connective Tissue Diseases"] <- "purple"
classes.cols["Urogenital Diseases and Pregnancy Complications"] <- "pink"

# identify the class color for each disease
disease.cols <- classes.cols[disease.classes]
names(disease.cols) <- names(disease.classes)

# create cell type diseasome
project.network(pvalues.gsc, col.vert=disease.cols, col.other=col.other, vert.size.max=5, vert.label.cex=0.6, edge.width.max=5, layout=layout.fruchterman.reingold) 
project.network.legend(classes.cols, col.other)

## ----create_monocyte_interactome-----------------------------------------
# create global network
data(edgelist.string)
g <- graph.edgelist(as.matrix(edgelist.string[, c("idA", "idB")]), directed=FALSE)

# load expression data
data(expression.fantom5)
expression <- expression.transform(expression.fantom5)

# create monocyte-specific interactome
g.monocyte <- score.edges(g, expression[, "monocyte"])

## ----identify_disease_genes----------------------------------------------
# identify all psoriasis-associated genes
data(disease.genes)
psoriasis.genes <- disease.genes[["psoriasis"]]

# remove genes not present in network and for which we have no expression data
psoriasis.genes <- psoriasis.genes[psoriasis.genes %in% V(g.monocyte)$name]
psoriasis.genes <- psoriasis.genes[psoriasis.genes %in% rownames(expression)]

# remove genes with more than 15 interacting partners
psoriasis.genes <- psoriasis.genes[degree(g.monocyte)[psoriasis.genes] <= 15]

## ----psoriasis_subgraph--------------------------------------------------
# produce psoriasis subgraph
# in the paper, 500 bins, rather than 20 are used
g.subgraph <- disease.subgraph(g.monocyte, psoriasis.genes, n.bins=20, edge.width.max=10)
plot(g.subgraph)

## ----edge_enrichment-----------------------------------------------------
# setup
edge.attr <- "score"
n.perm <- 20
edge.cutoff <- 0.9 # edges with weight higher than 0.9 are within the top 1% of monocyte edge
proportion.perm <- rep(NA, n.perm)

# create observed subgraph
g.subgraph.obs <- disease.subgraph(g.monocyte, psoriasis.genes, edge.attr=edge.attr, n.bins=20, edge.width.max=10)
proportion.obs <- mean(get.edge.attribute(g.subgraph.obs, edge.attr) > edge.cutoff)

for (n in 1:n.perm) {
    # create permuted subgraph
    expression.cell <- structure(sample(expression[, "monocyte"]), names=rownames(expression))
    g.monocyte.perm <- score.edges(g, expression.cell)
    g.subgraph.perm <- disease.subgraph(g.monocyte.perm, psoriasis.genes, edge.attr=edge.attr, n.bins=20, edge.width.max=10)
    
    # compute proportion of edges with weights created than edge.cutoff
    proportion.perm[n] <- mean(get.edge.attribute(g.subgraph.perm, edge.attr) > edge.cutoff)    
}

# compute empirical p-value
# we set the lower limit of the p-value as 1/n.perm
max(mean(proportion.perm >= proportion.obs), 1/n.perm)

## ----session_info--------------------------------------------------------
sessionInfo()

