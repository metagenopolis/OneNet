## code to prepare `liver` dataset goes here
counts.file<-"https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/9NGGMR"
taxo.file<- "https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/G8QQ0R"

#inferences on the entire population with prevalence 90
inf_prev90<-readRDS("data-raw/raw_allpop_inference_prev90.rds" )

#counts data on the ill population with prevalence 50
counts_prev50<-readRDS("data-raw/raw_illpop_counts.rds" )

#inferences on the ill population with prevalence 50
inf_prev50<-readRDS("data-raw/raw_illpop_inference.rds" )

#aggregated inferences on the ill population with prevalence 50
aggreg_data50=readRDS("data-raw/aggreg_illpop.rds")

taxonomy<- read.table(file = url(taxo.file), 
                                 header=TRUE, 
                                 stringsAsFactors = FALSE, 
                                 row.names=1) %>%
                                 mutate(msp_name=rownames(.)) %>%
                                 na.omit()

counts <- read.table(file=url(counts.file), 
                              header=TRUE, 
                              sep=",", 
                              dec=".", 
                              stringsAsFactors = FALSE, 
                              row.names=1) %>% 
                              select(-health_status) %>%
                              t() %>% 
                              data.frame() %>% 
                              mutate(id_mgs=rownames(.))%>% 
                              relocate(id_mgs)

meta <- read.table(file=url(counts.file), 
                      header=TRUE, 
                      sep=",", 
                      dec=".", 
                      stringsAsFactors = FALSE, 
                      row.names=1) %>% 
                      select(health_status)

liver <- list(
  abundances=counts,
  taxonomy=taxonomy,
  meta= meta,
  infer_prev90=inf_prev90,
  counts_prev_ill= counts_prev50,
  infer_prev_ill=inf_prev50,
  aggreg_ill= aggreg_data50,
  msp_set_50 = colnames(counts_prev50),
  msp_set_90 = as.character(taxonomy[rowMeans(counts > 0) > 0.9, "msp_name"])
)

usethis::use_data(liver,overwrite = TRUE, compress = "xz")

## code to prepare `liver` dataset goes here
# counts.file<-"https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/9NGGMR"
# taxo.file<- "https://entrepot.recherche.data.gouv.fr/api/access/datafile/:persistentId?persistentId=doi:10.57745/G8QQ0R"
# taxonomy<- read.table(file = url(taxo.file), 
#                       header=TRUE, 
#                       stringsAsFactors = FALSE, 
#                       row.names=1)%>%
#                        mutate(msp_name=rownames(.)) %>%
#                        na.omit()
# 
# counts <- read.table(file=url(counts.file),
#                      header=TRUE,
#                      sep=",",
#                      dec=".",
#                      stringsAsFactors = FALSE,
#                      row.names=1) %>%
#                      select(-health_status) %>%
#                      t() %>%
#                      data.frame() %>%
#                      mutate(id_mgs=rownames(.))%>%
#                      relocate(id_mgs)
# 
# counts_from_table<-get_count_table(abund.table=counts, prev.min=0.9, verbatim=FALSE)$data
# rep.num=20; n.levels=100; edge_thresh=0.9; cores=2
# some_inferences<-all_inferences(data=counts_from_table,edge_thresh=edge_thresh,
#                                 rep.num = rep.num,n.levels = n.levels, cores = cores,
#                                 Offset=TRUE)
#saveRDS(some_inferences,"data-raw/raw_allpop_inference_prev90.rds")

# counts_ill <- read.table(file=url(counts.file), 
#                          header=TRUE, 
#                          sep=",", 
#                          dec=".", 
#                          stringsAsFactors = FALSE, 
#                          row.names=1) %>% 
#                          filter(health_status=="P") %>% 
#                          select(-health_status) %>% 
#                          t() %>% 
#                          data.frame() %>% 
#                          mutate(id_mgs=rownames(.))%>% 
#                          relocate(id_mgs)
# 
# counts_from_table<-get_count_table(abund.table=counts_ill, prev.min=0.5, sample_id=colnames(counts_ill), verbatim=FALSE)$data
# saveRDS(counts_from_table,"raw_illpop_counts.rds")
# rep.num=20; n.levels=100; edge_thresh=0.9; cores=2
# some_inferences<-all_inferences(data=counts_from_table,edge_thresh=edge_thresh,
#                                 rep.num = rep.num,n.levels = n.levels, cores = cores,
#                                 Offset=TRUE)
# saveRDS(some_inferences,"raw_illpop_inference.rds")
# adapted<-adapt_mean_stability(some_inferences,mean.stability=0.8, plot=TRUE)
# final_frequences<-adapted$freqs
# aggreg_data=compute_aggreg_measures(final_frequences)
# saveRDS(aggreg_data,"aggreg_illpop.rds")
