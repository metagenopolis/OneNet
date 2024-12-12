
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OneNet

<!-- badges: start -->
<!-- badges: end -->

OneNet is a stability based network inference aggregation method from
whole metagenome sequencing data, which aims at fostering
reproducibility and precision.

License: MIT + file LICENSE

## Installation

You can install OneNet from the [forgemia gitlab](metagenopolis/OneNet)
with:

``` r
remotes::install_gitlab(repo = "metagenopolis/OneNet", 
                        host = "forgemia.inra.fr")
```

Or the development version with

``` r
remotes::install_gitlab(repo = "metagenopolis/OneNet@dev", 
                        host = "forgemia.inra.fr")
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
## basic example code
data("liver")
str(liver, max.level=1)
#> List of 5
#>  $ abundances  : chr [1:1990, 1:238] "msp_0001" "msp_0002" "msp_0003" "msp_0004" ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ metadata    : tibble [222 × 27] (S3: tbl_df/tbl/data.frame)
#>  $ infer_prev50:List of 7
#>  $ infer_prev90:List of 7
#>  $ taxonomy    :'data.frame':    1990 obs. of  21 variables:
```
