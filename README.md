
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OneNet

<!-- badges: start -->
<!-- badges: end -->

OneNet is a stability based network inference aggregation method from
whole metagenome sequencing data, which aims at fostering
reproducibility and precision.

## Installation

You can install the latest version of OneNet from the with the following
command :

``` r
remotes::install_github(repo = "metagenopolis/OneNet@*release")
```

## Tutorial

A tutorial on the liver cirrhosis dataset from [Qin et. al
2014](https://pubmed.ncbi.nlm.nih.gov/25079328/) is available. You need
to build package vignettes.

``` r
remotes::install_github(repo = "metagenopolis/OneNet@*release", build_vignettes = TRUE, force=TRUE)
library(OneNet)
vignette(package='OneNet', "Demo_offline")
```

## Citation

If you find OneNet useful, please cite:

Champion C, Momal R, Le Chatelier E, Sola M, Mariadassou M, Berland M
(2024) OneNet—One network to rule them all: Consensus network inference
from microbiome data. PLoS Comput Biol 20(12): e1012627.
<https://doi.org/10.1371/journal.pcbi.1012627>
