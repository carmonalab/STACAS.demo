# scRNA-seq data integration using STACAS

[This demo](https://carmonalab.github.io/STACAS.demo/tutorial.html) will run you through a complete dataset integration using Seurat 3 and STACAS. It reproduces the results and figures in [Andreatta & Carmona, Bioinformatics 2021](https://doi.org/10.1093/bioinformatics/btaa755)

We will be using the four following datasets:

* Cd8+ tumor infiltrating lymphocytes (TILs), from Carmona et al. GEO: [GSE116390](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116390)

* Cd8+/Cd4+  (TILs), from Xiong et al. ArrayExpress: [E-MTAB-7919](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7919/) 

* Cd4+ TILs, from Magen et al. GEO: [GSE124691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124691)

* Cd4+ T cells from draining lymph nodes, from Magen et al. GEO: [GSE124691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124691)


### How to play the demo

Copy this repository to your local system:
```
git clone --depth 1 https://gitlab.unil.ch/carmona/STACAS.demo.git
```
Then move to the `STACAS.demo` folder and open the `STACAS.demo.Rproj` in RStudio. In this enviroment, you can follow step-by-step the commands outlined in `STACAS.vignette.Rmd`

If you have problems with the installation, you may also download a [Docker image](https://hub.docker.com/r/mandrea1/stacas_demo) of the project with all dependencies pre-installed.

See also a precompiled html version of this demo at: [STACAS Tutorial](https://carmonalab.github.io/STACAS.demo/tutorial.html)

For installation and documentation on the `STACAS` package refer to the [STACAS Github repository](https://github.com/carmonalab/STACAS)
