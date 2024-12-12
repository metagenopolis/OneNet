## code to prepare `liver` dataset goes here
counts.file<-"https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/9NGGMR"
taxo.file<- "https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/G8QQ0R"
#inferences on all population
inf_prev90<-readRDS("data-raw/raw_allpop_prev90.rds" )
#aggregated inferences on all population
aggreg_data90=readRDS("data-raw/aggreg_data_90.rds")
aggreg_data50=readRDS("data-raw/aggreg_data_50.rds")
taxonomy<- read.table(file = url(taxo.file), 
                                 header=TRUE, 
                                 stringsAsFactors = FALSE, 
                                 row.names=1)
counts <- read.table(file=url(counts.file), 
                              header=TRUE, 
                              sep=",", 
                              dec=".", 
                              stringsAsFactors = FALSE, 
                              row.names=1) %>% select(-health_status)
liver <- list(
  abundances=counts,
  metadata=meta,
  infer_prev90=inf_prev90,
  aggreg_data90= aggreg_data90,
  aggreg_data50= aggreg_data50,
  taxonomy=taxonomy,
  msp_set_50 = as.character(taxonomy[rowMeans(counts > 0) > 0.5, "msp"]),
  msp_set_90 = as.character(taxonomy[rowMeans(counts > 0) > 0.9, "msp"])
)

usethis::use_data(liver,overwrite = TRUE)



# data("liver")
# abund<-liver$abundances
# meta<-liver$metadata
# taxo=liver$taxonomy
# counts_from_table<-get_count_table(abund.table=abund, mgs=taxo$msp_name, prev.min=0.9, verbatim=FALSE,sample_id = meta$SRA)$data
# some_inferences<-all_inferences(counts_from_table, rep.num=100,edge_thresh = 0.9,n.levels = 100,offset = TRUE,covar=NULL,
#                                 cores=25)
# saveRDS(some_inferences,"data-raw/raw_allpop_prev90.rds")
# adapted<-adapt_mean_stability(some_inferences,mean.stability=0.9, plot=TRUE)
# final_frequences<-adapted$freqs
# aggreg_data=compute_aggreg_measures(final_frequences)
# saveRDS(aggreg_data,"data-raw/aggreg_data_90.rds")

