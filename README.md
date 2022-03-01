# Endometriosis scRNA-seq

This repository accompanies the manuscript
> Tan Y, Flynn WF, Sivajothi S, Bozal S, Luo D, Luciano AA, Robson P, Luciano
> DE, Courtois ETC. "Cellular Phenotyping of Endometriosis at the single cell
> level highlights complex heterogeneity" (2020)
currently available on [bioRxiv][biorxiv_link].


## Data Availability

### Raw and processed data
Currently this data is under embargo at GEO at [GSE179640][geo_link] to be
released upon publication.

### Analyzed data objects
Additionally upon publication, the analyzed objects (in `.h5ad` format) will be
available for download or exploratory viewing at the 
[JAX Single Cell portal][portal_link].


## Issues with this code or repository
We are committed to making this analysis and its accompanying figures
reproducible.  If you experience problems, a bug, a question, etc., please open
an issue and tag @wflynny or @yulianatan.


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
    sys.path.append(f"{os.getcwd()}/../code")
    from plotting.plot_settings import *
    from plotting.palettes import *
    from plotting.scatter import *
    from plotting.matrix import *
    ```


[biorxiv_link]: https://www.biorxiv.org/content/10.1101/2021.07.28.453839v1
[geo_link]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE179640
[portal_link]: https://singlecell.jax.org
