# Endometriosis scRNA-seq

This repository accompanies the manuscript
> Tan Y, Flynn WF, Sivajothi S, Bozal S, Luo D, Luciano AA, Robson P, Luciano
> DE, Courtois ETC. "Cellular Phenotyping of Endometriosis at the single cell
> level highlights complex heterogeneity" (2020)
currently available on [bioRxiv][biorxiv_link].


## For development

1.  Clone this repository somewhere with
    ```{bash}
    cd /path/for/manuscript/code
    git clone https://github.com/TheJacksonLaboratory/endometriosis-scrnaseq.git
    ```
2.  Link for data objects:
    ```{bash}
    cd endometriosis-scrnaseq/data
    ln -s /path/to/global_20201234.h5ad global.h5ad
    ln -s /path/to/fibroblasts_20201234.h5ad fibro.h5ad
    ...
    ```
3.  Make a notebook in the notebooks folder.
4.  In that notebook, do the following after your standard imports:
    ```{python}
    import os, sys
    sys.path.append(f"{os.getcwd()/../code/plotting")
    from plot_settings import *
    from palettes import *
    from scatter import *
    from matrix import *
    ```


[biorxiv_link]: https://biorxiv.org
