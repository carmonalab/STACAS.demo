# scRNA-seq data integration using STACAS


STACAS is a method for scRNA-seq integration. It is based on the [Seurat](https://cran.r-project.org/web/packages/Seurat/index.html) integration framework, but adds important innovations:


* **anchors can be filtered/down-weighted** based on their distance in reciprocal PCA space
* **integration trees** are constructed based on the 'centrality' of datasets, as measured by the sum of their anchor weights
* **cell type labels**, if known, can be given as input to the algorithm to perform **semi-supervised integration**

In [this demo](https://carmonalab.github.io/STACAS.demo/STACAS.demo.html) we will show the application of STACAS to integrate a collection of PBMC datasets.

See also how STACAS integration compares to Seurat, and how it avoids overcorrecting batch effects for heterogeneous data sets, at this [demo for T cell data integration](https://carmonalab.github.io/STACAS.demo/Tcell.demo.html)


### How to play the demos

Copy this repository to your local system:
```
git clone --depth 1 https://gitlab.unil.ch/carmona/STACAS.demo.git
```
Then move to the `STACAS.demo` folder and open the `STACAS.demo.Rproj` in RStudio. In this enviroment, you can follow step-by-step the commands outlined in `STACAS.demo.Rmd`


For installation and documentation on the `STACAS` package refer to the [STACAS Github repository](https://github.com/carmonalab/STACAS)
