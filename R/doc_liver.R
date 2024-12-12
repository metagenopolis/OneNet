#' Liver cirrhosis data set
#'
#' This data set contains abundance data, species taxonomy and sample metadata from [Qin et. al 2014](https://pubmed.ncbi.nlm.nih.gov/25079328/).
#' The liver data set also contains all results from 7 inference methods used on a filtered subset consisting of species with 
#' prevalence higher than 90% (`infer_prev90`). Those results have been processed with `adapt_mean_stability()` followed by 
#' `compute_aggreg_measures()` to produce consensus networks (`aggreg_data90`). The same was done for taxa with prevalence higher
#' than 50% (`aggreg_data50`).
#'
#' The dataset has been cleaned to remove duplicates and contaminated samples.
#'
#' @name liver
#' @docType data
#' @source [abundance](https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/9NGGMR) and [taxonomy table](https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/G8QQ0R)
#' @references Qin N and al. (2014) Alterations of the Human Gut Microbiome in Liver Cirrhosis. Nature, 513(7516), 59–64.
#' See a data description at [Qin et. al 2014](https://pubmed.ncbi.nlm.nih.gov/25079328/).
#' @format
#' A list with the following components:
#' * abundances: a 1990 x 216 matrix with (normalized) abundances of 1990 MSPs in the 216 samples of the Qin et al. dataset
#' * metadata: a 216 x 27 data.frame with 27 descriptors for each sample in the dataset
#' * infer_prev90: list of outputs reconstructed for different regularization parameter (lambda) values using 7 different methods (PLNnetwork, gCoda, SpiecEasi, SPRING, Magma, EMtree, ZiLN), using only MSPs with prevalence higher than 90%
#' * aggreg_data90: A 903 x 4 data.frame containing for each potential edge (row) between high prevalence (>90%) species its summary statistics (mean, norm-2, IVW) and the number of methods for which
#'    the edge has high inclusion frequency (> 0.9). Computed from `infer_prev90`.
#' * aggreg_data50: A 12403 x 4 data.frame containing for each potential edge (row) between high prevalence (>50%) species its summary statistics (mean, norm-2, IVW) and the number of methods for which
#'    the edge has high inclusion frequency (> 0.9).Computed from `infer_prev50` (not provided). 
#' * taxonomy: a 1990 x 21 data.frame with descriptors of the MSPs
#' * msp_set_50: character vector of the 158 MSPs with prevalence higher than 50%
#' * msp_set_90: character vector of the 43 MSPs with prevalence higher than 90%
#'
#' @keywords data
"liver"
