
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

## Liver data set

This package includes a data set “liver”, which contains

- Abundances (semi-discrete) and metadata for a trial on healthy and
  liver cirrhosis patients from [Qin et. al
  2014](https://pubmed.ncbi.nlm.nih.gov/25079328/).
- A taxonomy table.
- Network inferences with 6 methods for a minimal species prevalence of
  90% and 50%.

``` r
library(OneNet)
#> Le chargement a nécessité le package : ggplot2
#> Le chargement a nécessité le package : SpiecEasi
#> Registered S3 methods overwritten by 'lava':
#>   method    from
#>   plot.sim  huge
#>   print.sim huge
## basic example code
data("liver")
str(liver, max.level=1)
#> List of 9
#>  $ abundances     :'data.frame': 1990 obs. of  217 variables:
#>  $ taxonomy       :'data.frame': 1990 obs. of  3 variables:
#>   ..- attr(*, "na.action")= 'omit' Named int 1991
#>   .. ..- attr(*, "names")= chr "health_status"
#>  $ meta           :'data.frame': 216 obs. of  1 variable:
#>  $ infer_prev90   :List of 7
#>  $ counts_prev_ill: num [1:114, 1:122] 275 0 4487 11 7 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ infer_prev_ill :List of 7
#>  $ aggreg_ill     :'data.frame': 7381 obs. of  4 variables:
#>  $ msp_set_50     : chr [1:122] "msp_0003" "msp_0005" "msp_0007" "msp_0008" ...
#>  $ msp_set_90     : chr [1:26] "msp_0005" "msp_0007" "msp_0010" "msp_0011" ...
```

## Citation

If you find OneNet useful, please cite:

Champion C, Momal R, Le Chatelier E, Sola M, Mariadassou M, Berland M
(2024) OneNet—One network to rule them all: Consensus network inference
from microbiome data. PLoS Comput Biol 20(12): e1012627.
<https://doi.org/10.1371/journal.pcbi.1012627>
