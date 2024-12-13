#' Count Table
#'
#' Get count table from MGS normalized abundances.
#'
#' Reads data from a provided path to MGP server or the abundance table directly.
#' Extracts species and transforms abundances into counts by setting the minimal positive
#' observed value to 1.
#'
#' @param abund.path Path to the abundance table in MGP server.
#' @param sample_id Ids of samples to keep in the final table.
#' @param prev.min Double between 0 and 1. Minimal prevalence threshold of MGS to keep in final
#' table.
#' @param abund.table Abundance table. By default this table should have the mgs names as first
#' column, and be in tsv or rds format.
#' @param mgs Character vector giving MGS names, if they are not specified in the abundance table
#'  first column.
#' @param verbatim Boolean controlling verbosity.
#'
#' @return
#' A list containing
#' \itemize{
#'  \item{data: }{the final count table (tibble).}
#'  \item{prevalences: }{a tibble gathering the prevalence of each MGS.}
#'}
#' @export
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate filter select
#' @importFrom utils tail read.delim
#' @examples
#' # With table
#' data("liver")
#' abund<-liver$abundances
#' meta=liver$metadata
#' taxo=liver$taxonomy
#' count_from_table<-get_count_table(abund.table=abund,prev.min=0.9, verbatim=FALSE)
#' str(count_from_table)
get_count_table<-function(abund.path=NULL,sample_id=NULL,abund.table=NULL,prev.min,verbatim=TRUE,
                          mgs=NULL){
  if(!is.null(abund.path)){
    metaformat<-tail(strsplit(abund.path,"[.]")[[1]],1)
    if(metaformat=="rds"){
      data_table <- as.matrix(readRDS(abund.path))
    }else{
      data_table <- as.matrix(read.delim(abund.path))
    }
  }else{
    data_table<-abund.table
  }
  if(is.null(mgs)){
    species<-data_table[,1]
    data_table<-data_table[,-1]
  }else{
    species<-mgs
  }

  data_table<-apply(data_table,2,function(x){ as.numeric(x)})
  pct_zero_before=sum(data_table==0)/length(data_table)
  data_table<-as_tibble(data_table)
  #make sure we have samples with metadata
  if(!is.null(sample_id)){
    data_table<-data_table[,colnames(data_table)%in%sample_id]
  }

  #prevalence filtering
  counts<-data_table%>% mutate(prev=rowMeans(.>0),mgs=species)   %>% filter(prev>prev.min)
  prevalences<-counts %>% select(prev, mgs)
  counts=counts %>% select(-mgs, -prev) %>% t(.)
  colnames(counts)<-prevalences$mgs
  #set minimal positive value to 1
  cst=1/min(counts[counts!=0])
  counts=round(counts*cst)
  #Verbatim
  pct_zero_after=sum(counts==0)/length(counts)
  if(verbatim){cat(paste0("Preprocessing step output for species prevalence>",
                          prev.min*100,"% : \n   -from ", length(species),
                          " to ", ncol(counts)," species","\n   -from ",
                          round(pct_zero_before*100,1),"% to ",
                          round(pct_zero_after*100,1),"% zero values."))}

  return(list(data=counts, prevalences=prevalences))
}



#' Perform network inference with several methods
#'
#' @param data Count data set.
#' @param rep.num Number of subsamples.
#' @param n.levels Size of probability grid.
#' @param cores Number of cores for parallel computation.
#' @param edge_thresh Threshold on selection frequencies to set the number of selected edges.
#' @param covar Covariate data set.
#' @param Offset Boolean for the computation of offset terms. If TRUE uses GMPR, or RLE on failure.
#' @param methods Vector of characters indicating with which methods the inference should be
#'  performed among "PLNnetwork","SpiecEasi","gCoda","EMtree","Magma", "ZiLN", and "SPRING".
#' @param seed A numerical seed for subsampling (required by the SPRING package).
#' @param subsample.ratio  The subsampling ratio. The default value is 10*sqrt(n)/n when n>144 and 0.8 when n<=144, where n is the sample size.
#'
#' @return A list of outputs from the desired individual methods.
#' @export
#' @importFrom SpiecEasi sparseiCov
#'
#' @examples
#' data("liver")
#' abund<-liver$abundances
#' meta=liver$metadata
#' taxo <- liver$taxonomy
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9, verbatim=FALSE)$data
#' some_inferences<-all_inferences(counts_from_table, rep.num=4)
#' str(some_inferences, max.level=1)
all_inferences <-function(data,edge_thresh=0.9, rep.num=100, n.levels=100,
                         Offset=TRUE, covar=NULL, cores=1,
                         methods=c("PLNnetwork","SpiecEasi","gCoda","EMtree","Magma","SPRING","ZiLN"),
                         seed=10010,subsample.ratio=NULL){
  if(!is.null(covar)) covar<-data.frame(covar)

  #individual inferences
  res<-list()
  if("PLNnetwork"%in%methods){
    output_PLN=just_PLN(data, rep.num, n.levels, edge_thresh, Offset,covar=covar,subsample.ratio=subsample.ratio)
    res$PLNnetwork<-output_PLN
  }
  if("gCoda"%in%methods){
    if(exists("output_PLN")){
      lambda.seq<-output_PLN$lambda.seq
    }else{lambda.seq<-NULL}
    output_gcoda=just_gcoda(data, rep.num,n.levels, cores, edge_thresh,
                            lambda.seq=,covar=covar,subsample.ratio=subsample.ratio)
    res$gCoda<-output_gcoda
  }

  if("SpiecEasi"%in%methods){
    output_spiec=just_spiec(data, rep.num, n.levels, cores,edge_thresh, covar=covar,subsample.ratio=subsample.ratio)
    res$SpiecEasi<-output_spiec
  }
  if("SPRING"%in%methods){
    output_spring=just_spring(data, rep.num, n.levels, cores,edge_thresh, covar=covar, seed=seed,subsample.ratio=subsample.ratio)
    res$SPRING<-output_spring
  }
  if("Magma"%in%methods){
    output_magma=just_magma(data, rep.num, n.levels, edge_thresh, Offset,covar=covar,subsample.ratio=subsample.ratio)
    res$Magma<-output_magma
  }
  if("EMtree"%in%methods){
    output_EMtree=just_EMtree(data, rep.num,n.levels, edge_thresh,Offset,subsample.ratio=subsample.ratio,
                              covar=covar)
    res$EMtree<-output_EMtree
  }

  if("ZiLN"%in%methods){
    output_ZiLN=just_ZiLN(data, rep.num,n.levels, cores, edge_thresh,subsample.ratio=subsample.ratio)
    res$ZiLN<-output_ZiLN
  }

  return(res)
}


#' Adapt mean stability
#'
#'Adapts each method's stability via the edge density, so that the mean stability is
#' the desired value.
#'Adapting the mean stability of the set of methods aims at adapting each method to
#'the same precision level.
#'
#' @param inference_collection Collection of inferences, result from the all_inferences() function.
#' @param mean.stability Value for the mean stability.
#' @param plot Boolean.
#'
#' @return A list containing
#' \itemize{
#'  \item{freqs: }{a tibble gathering the edge selection frequencies corresponding to
#'  the adapted stability of each method.}
#'  \item{stab_data: }{a tibble giving the adapted stability and density of each method.}
#'}
#'
#' @export
#' @importFrom ggplot2 aes ggplot geom_hline geom_point geom_line geom_vline coord_cartesian labs theme_light
#' @importFrom dplyr group_by filter row_number pull select
#' @importFrom tibble as_tibble
#' @importFrom EMtree freq_selec
#' @importFrom SpiecEasi sparseiCov
#' @examples
#' data("liver")
#' abund<-liver$abundances
#' meta=liver$metadata
#' taxo <- liver$taxonomy
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9, verbatim=FALSE)$data
#' some_inferences<-liver$infer_prev90
#' adapted_stab_data<-adapt_mean_stability(some_inferences, mean.stability=0.8)
adapt_mean_stability<-function(inference_collection,mean.stability=0.8,
                               plot=FALSE){
  ne=inference_collection[[1]]$nedges
  methods<-names(inference_collection)
  #--- number of edges for mean stability
  lambda_stab_density= do.call(rbind,lapply(seq_along(inference_collection), function(index_method){
    cbind(inference_collection[[index_method]]$lambda_stab,method=methods[index_method])
  }))

  evol_adapt=do.call(rbind,lapply(seq(10,ne/4,10), function(ntest){
    filter_stab_density=lambda_stab_density%>%group_by(method) %>%
      filter(abs(nedges-ntest)==min(abs(nedges-ntest)))
    res= data.frame(mstab=mean(filter_stab_density$stability), nedge=ntest)
    return(res)
  }))
  nedge.min=evol_adapt$nedge[which.min(evol_adapt$mstab)]
  evol_adapt=evol_adapt %>% filter(nedge<nedge.min)
  nedge.final = evol_adapt$nedge[which.min(abs(evol_adapt$mstab-mean.stability))]
  if(plot){
    g=evol_adapt %>% as_tibble() %>%
      ggplot(aes(nedge,mstab))+
      geom_hline(yintercept =mean.stability, color="gray70",linetype="dashed")+
      geom_point()+geom_line()+
      geom_vline(xintercept = nedge.final, color="darkorange",linetype="dashed")+
      coord_cartesian(xlim=c(0,min(max(evol_adapt$nedge),nedge.final*3)))+
      labs(x="Density",y="Mean stability")+
      theme_light()
    print(g)
  }
  adapt_stab=lambda_stab_density%>%group_by(method) %>%
    filter(abs(nedges-nedge.final)==min(abs(nedges-nedge.final))) %>% filter(row_number()==1)

  #--- gather method's optimal frequencies in all_optim_freqs

  all_optim_freqs<-do.call(cbind,lapply(seq_along(inference_collection), function(num){
    meth<-methods[num]
    lambda<-adapt_stab$lambda[adapt_stab$method==meth]
    freq<-suppressMessages(get_vec(lambda,inference_collection[[num]],meth))
    return(freq)
  })) %>% 
    `colnames<-`(methods) %>%
    as_tibble()

  adapt_stab<-adapt_stab %>% dplyr::select(-lambda)
  return(list(freqs=all_optim_freqs, stab_data=adapt_stab))
}

#' Perform network inference with several methods
#'
#' @param data Count data set.
#' @param rep.num Number of subsamples.
#' @param methods Vector of characters indicating with which methods the inference should be
#'  performed among "PLNnetwork","SpiecEasi","gCoda","EMtree","Magma", "ZiLN", and "SPRING".
#' @param fast Boolean for the fast computation of OneNet (without PLNnetwork)
#' @param parallel Boolean for the computation of OneNet with all methods parallelized
#' @inheritDotParams just_EMtree 
#' @return A list of outputs from the desired individual methods.
#' 
#' @details
#' For details of the arguments used by the different methods, see the just_*() functions (e.g. [just_PLN()] and [just_EMtree()])
#' 
#' @export
#' @importFrom SpiecEasi sparseiCov
#'
#' @examples
#' data("liver")
#' abund<-liver$abundances
#' meta=liver$metadata
#' taxo <- liver$taxonomy
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9, verbatim=FALSE)$data
#' some_inferences<-all_inferences_new(data=counts_from_table, rep.num=4)
#' str(some_inferences, max.level=1)
all_inferences_new <-function(data, rep.num=100, 
                              methods=c("PLNnetwork","SpiecEasi","gCoda","EMtree","Magma","SPRING","ZiLN"),
                              fast=FALSE, parallel=FALSE, ...){
  if(fast){
    methods=setdiff(methods, "PLNnetwork")
  }
  if(parallel){
    res=mclapply(methods, apply_method, mc.cores=length(methods), data, rep.num,  ...)
    names(res)=methods
  }else{
    res=lapply(methods, apply_method, data, rep.num,...)
    names(res) = methods
  }
  
  return(res)
}

