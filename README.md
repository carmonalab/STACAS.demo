# scRNA-seq data integration using STACAS

[STACAS](https://github.com/carmonalab/STACAS) is a method for scRNA-seq integration. It is based on the [Seurat](https://cran.r-project.org/web/packages/Seurat/index.html) integration framework, but adds important innovations:

* **anchors are down-weighted** based on their distance in reciprocal PCA space
* **integration trees** are constructed based on the 'centrality' of datasets, as measured by the sum of their re-weighted anchor scores
* **cell type labels**, if known, can be given as input to the algorithm to perform **semi-supervised integration**

In [this demo](https://carmonalab.github.io/STACAS.demo/STACAS.demo.html) we will show the application of STACAS to integrate a collection of human immune cell datasets.

See also how STACAS compares to other integration methods, and how it avoids overcorrecting batch effects for heterogeneous data sets, at this [demo for T cell data integration](https://carmonalab.github.io/STACAS.demo/Tcell.demo.html)


### How to play the demos

Copy this repository to your local system:
```
git clone https://github.com/carmonalab/STACAS.demo
```
Then move to the `STACAS.demo` folder and open the `STACAS.demo.Rproj` in RStudio. In this enviroment, you can follow step-by-step the commands outlined in `STACAS.demo.Rmd`


For installation and documentation on the `STACAS` package refer to the [STACAS Github repository](https://github.com/carmonalab/STACAS)

### Citation

Massimo Andreatta, Santiago J Carmona, STACAS: Sub-Type Anchor Correction for Alignment in Seurat to integrate single-cell RNA-seq data, Bioinformatics, 2020, btaa755, https://doi.org/10.1093/bioinformatics/btaa755
