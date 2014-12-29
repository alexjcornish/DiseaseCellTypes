DiseaseCellTypes
==========

An R-package for identifying disease-associated cell types using cell type-specific interactomes. This package provides 2 methods for identifying disease-cell type associations: gene set compactness (GSC) and gene set overexpression (GSO). Also provided are functions for building cell type-specific interactomes and cell type-based diseasomes.

How to run the methods and recreate the results of the referenced paper is described in the [package vignette][1]. Further details on each of the package functions are available in the [package manual][2]. Both the [package vignette][1] and [package manual][2] are contained within the package. 


Additional data
----------

- [Disease-cell type associations][3]
- [Cell type-specific interactomes][4]


Download package binaries
----------

- [Mac OS/X binary][5]
- [Windows binary][6]
- [Source package][7]


Development
----------

The development version of DiseaseCellTypes can be downloaded and installed using the devtools R-package. DiseaseCellTypes depends on a number of other R-packages, including gplots, igraph, Matrix, psych, snow and stringr.

```
#install.packages("devtools")
require(devtools)
install_github("alexjcornish/DiseaseCellTypes")
```


References
----------

Paper under preparation.


[1]: https://cdn.rawgit.com/alexjcornish/DiseaseCellTypes/master/inst/doc/DiseaseCellTypes-vignette.html?raw=TRUE
[2]: https://github.com/alexjcornish/DiseaseCellTypes/blob/master/inst/doc/DiseaseCellTypes-manual.pdf?raw=TRUE
[3]: http://alexjcornish.github.io/Disease_Cell_Association_Data/
[4]: http://alexjcornish.github.io/Cell_Type_Interactomes/
[5]: https://github.com/alexjcornish/DiseaseCellTypes_Binaries/blob/master/DiseaseCellTypes_0.9.0.tgz?raw=TRUE
[6]: https://github.com/alexjcornish/DiseaseCellTypes_Binaries/blob/master/DiseaseCellTypes_0.9.0.zip?raw=TRUE
[7]: https://github.com/alexjcornish/DiseaseCellTypes_Binaries/blob/master/DiseaseCellTypes_0.9.0.tar.gz?raw=TRUE

