DiseaseCellTypes
===

An R-package for the identification of disease-associated cell types using cell type-specific interactomes. This package provides 2 methods for identifying disease-cell type associations: gene set compactness (GSC) and gene set overexpression (GSO). Also provided are functions for building cell type-specific interactomes and cell type-based diseasomes.

Details on how to run the methods and recreate the results of the referenced paper are provided in the package vignette (note, this html file may need to be saved and reopened if it does not initially display correctly): https://github.com/alexjcornish/DiseaseCellTypes/blob/master/inst/doc/DiseaseCellTypes-vignette.html?raw=TRUE

Further details on each of the functions contained within the packages are available in the package manual: https://github.com/alexjcornish/DiseaseCellTypes/blob/master/inst/doc/DiseaseCellTypes-manual.pdf?raw=TRUE


Additional data
===========

Disease-cell type associations: http://alexjcornish.github.io/Disease_Cell_Association_Data/

Cell type-specific interactomes: http://alexjcornish.github.io/Cell_Type_Interactomes/


Download package
===========

Stable package binaries can be downloaded from the following URLs. 

Mac OS/X binary: https://github.com/alexjcornish/DiseaseCellTypes_Binaries/blob/master/DiseaseCellTypes_0.9.0.tgz?raw=TRUE

Windows binary: https://github.com/alexjcornish/DiseaseCellTypes_Binaries/blob/master/DiseaseCellTypes_0.9.0.zip?raw=TRUE

Source package: https://github.com/alexjcornish/DiseaseCellTypes_Binaries/blob/master/DiseaseCellTypes_0.9.0.tar.gz?raw=TRUE


Development
===========

The development version of DiseaseCellTypes can be downloaded and installed using the devtools R-package. DiseaseCellTypes depends on a number of other R-packages, including gplots, igraph, Matrix, psych, snow and stringr.

```
#install.packages("devtools")
require(devtools)
install_github("alexjcornish/DiseaseCellTypes")
```


References
===========

Paper under preparation.
