% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{clustering_SBM}
\alias{clustering_SBM}
\title{Perform clustering of square matrices according to Stochastic Bloc Modelling}
\usage{
clustering_SBM(
  edge_freq,
  freq.thresh = 0.9,
  selected_taxo,
  palette = palette_default()
)
}
\arguments{
\item{edge_freq}{Edge final frequency, typically the result of an inference aggregation.}

\item{freq.thresh}{Threshold of edge detection.}

\item{selected_taxo}{Taxonomic information about the MGS of the data set.}

\item{palette}{vector of colours defining the cluster palette.}
}
\value{
A list containing
\itemize{
\item{groups: }{vector of group membership.}
\item{species_groups: }{list (or tibble if as.tibble=TRUE) of each group's taxonomic information (name of msp, species,
genus, and phylum).}
\item{graph: }{graphical output, the network with nodes coloured according to their group.}
}
}
\description{
Perform clustering of square matrices according to Stochastic Bloc Modelling
}
\examples{
data("liver")
abund<-liver$abundances
meta<-liver$metadata
taxo=liver$taxonomy \%>\% dplyr::select(msp_name, species)
counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9,verbatim=TRUE)$data
species<-liver$msp_set_90
inf_90<-liver$infer_prev90
adapt<-adapt_mean_stability(inf_90, mean.stability=0.8,plot=TRUE)
selected_taxo=dplyr::left_join(data.frame(msp_name=species), taxo, by="msp_name")
mean<-compute_aggreg_measures(adapt$freqs)$mean
clustering_output<-clustering_SBM(mean, selected_taxo=selected_taxo)
}
