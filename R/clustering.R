#' ggplot heatmap with personalized ordering
#'
#' @param data A data table.
#' @param no.names Boolean.
#' @param no.xnames Boolean.
#' @param order Order for both rows and columns (only for square matrices).
#' @param row.order Order specific to rows.
#' @param col.order Order specific to columns.
#'
#' @return An ordered ggplot heatmap
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes theme_light geom_tile guides labs theme element_text coord_fixed coord_flip element_blank
#' @noRd
ggimage<-function(data, no.names=FALSE,no.xnames=FALSE, order=NULL, row.order=NULL,
                  col.order=NULL){

  num=FALSE
  p=ncol(data) ; n=nrow(data)
  if(is.null(col.order)) col.order=1:p
  if(is.null(row.order)) row.order=1:n
  melted_data <- melt(as.matrix(data))
  if(!is.null(order)){
    if(n==p) row.order =   col.order = order
  }
  if(is.null(colnames(data))){
    col.level=1:p
  }else{col.level=colnames(data)}

  if(is.null(rownames(data))){
    row.level=1:n
  }else{row.level=rownames(data)}

  melted_data$Var1 <- factor( melted_data$Var1, levels = row.level[row.order])
  melted_data$Var2 <- factor( melted_data$Var2, levels = col.level[col.order])

  g=ggplot(melted_data, aes(x=.data$Var1, y=.data$Var2, fill=.data$value)) + theme_light()+
    labs(x="",y="")+ geom_tile() +guides(fill=FALSE)+
    theme(plot.title = element_text(size=10, hjust=0.5))
  if(no.names){
    g=g+theme(axis.text=element_blank(),axis.ticks =element_blank())+coord_flip()
  }else{
    if(no.xnames){
      g=g+theme(axis.text.y=element_blank(),axis.ticks.y =element_blank())+coord_flip()
    }
  }
  g
}


#' Membership from SBM
#'
#' @param param.sbm Output value from estim_sbm().
#' @param force.groups Number of groups when forced. Default to NULL.
#' @param G Square matrix being clustered.
#' @param plot Boolean.
#'
#' @return Vector of group membership.
#' @noRd
get_cluster_groups<-function(param.sbm,force.groups=NULL, G, plot=FALSE){
  if(is.null(force.groups)){
    nb_groups<-param.sbm$best
  }else{nb_groups=force.groups}
  groups=get_groups(param.sbm$estim.sbm, k=nb_groups)
  sizes=unlist(lapply(groups, length))
  if(plot) ggimage(G,order = unlist(groups[order(sizes)]), no.names = TRUE)
  #assign groups to species
  cluster_groups<-rep(0,ncol(G))
  lapply(seq_along(groups),function(num.group){
    cluster_groups[groups[[num.group]]]<<-num.group
  })
  return(cluster_groups)
}
#' Extracts SBM groups
#'
#' @param estim.sbm Estimated SBM object
#' @param k Number of groups to extract
#'
#' @return List of nodes indices per group
#'
#' @noRd
get_groups<-function(estim.sbm, k){
  paramEstimSBMBern <- extractParamBM(estim.sbm,k)
  #--- extract k groups from inferred probabilities
  clusters= lapply(1:k, function(z){
    (which(paramEstimSBMBern$Z==z))
  })
  return(clusters)
}

#' Graph with group colouring
#'
#' Wrapper of EMtree::draw_network, which uses ggraph and tidygraph.
#'
#' @param G Adjacency matrix, possibly weighted.
#' @param groups Groups on nodes, possibly the result from a clustering.
#' @param stored_layout Optional layout (table of coordinates x and y).
#' @param legend Boolean for legend printing.
#' @param names Nodes labels.
#'
#' @return A graph with nodes coloured according to their group.
#' @importFrom EMtree draw_network
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices rainbow
#' @noRd
coloured_graph<-function(G, groups,stored_layout=NULL, legend=TRUE, names=NULL){
  if(!is.null(names)){
    size=3
  }else{size=4}
  nb_groups=length(unique(groups))
  if(nb_groups<=8){
    pal=brewer.pal(nb_groups, "Dark2")
  }else{pal= rainbow(nb_groups)}

  g<-draw_network(G,layout="nicely",stored_layout = stored_layout, btw_rank =1,label_size = size,
                  nodes_size = rep(2, nb_groups), remove_isolated = FALSE,
                  node_groups = groups, pal_nodes = pal,
                  pal_edges = "#31374f", legend=legend, nodes_label = names)$G
  return(g)
}

#' Extract spcies information in each group
#'
#' @param selected_taxo  Taxonomic information about the MGS of the data set.
#' @param groups Vector of group membership.
#' @param plot Boolean for graphic output.
#'
#' @return list of taxonomic information for each group.
#' @importFrom dplyr select group_by summarise mutate filter
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal labs
#' @importFrom tibble as_tibble
#' @noRd
groups_features<-function(selected_taxo,  groups, plot=TRUE){
  taxo_groups<-selected_taxo %>%as_tibble() %>% mutate(group=groups)
  #groups=vecteur des groupes des espèces
  if(plot){
    g<-taxo_groups %>%
      select(phylum, group) %>% group_by(group,phylum) %>%
      summarise(value=n()) %>%
      ggplot(aes(x=group, y=value, fill=phylum))+
      geom_bar(stat="identity")+theme_minimal()+labs(x="",y="")
    print(g)
  }
  groups_species<-lapply(1:max(groups), function(num){
    taxo_groups %>% filter(group==num) %>% select(msp_name,species,genus, phylum)
  })
  return(groups_species)
}
#' Wrapper for extracting SBM estimated parameters
#'
#' @param BMobject object from the estimation of SBM
#' @param k Number of groups
#'
#' @return SBM parameters
#' @noRd
extractParamBM <- function(BMobject,k){
  model <- BMobject$model_name
  membership_name <-  BMobject$membership_name
  res <- list()
  if (model == 'bernoulli') { res$alpha <- BMobject$model_parameters[k][[1]]$pi}
  if ((membership_name == 'SBM') |  (membership_name == 'SBM_sym')) {
    res$tau <-  BMobject$memberships[[k]]$Z
    res$Z <- apply(res$tau, 1, which.max)
    n <- nrow(BMobject$memberships[[k]]$Z)
    res$pi <-  colSums(BMobject$memberships[[k]]$Z)/n
    res$k <- length(res$pi)
  }
  ########## ordering
  if(k>1){if ((membership_name == 'SBM') |  (membership_name == 'SBM_sym')) {
    o <- switch(model,
                bernoulli  =  order(res$alpha %*% matrix(res$pi,ncol = 1),decreasing = TRUE),
                1:res$k
    )
    res$pi <- res$pi[o]
    res$alpha <- res$alpha[o,o]
    res$tau <- res$tau[,o]
    res$Z <- apply(res$tau, 1, which.max)
  }}
  return(res)
}

#' Estimates the SBM
#'
#' @param Amat Square matrix to cluster.
#' @param cores Number of cores for parallel computation of SBM.
#'
#' @return a list containing parameters of the estimation of SBM, and the index for the best model.
#' @importFrom blockmodels BM_bernoulli
#' @noRd
estim_SBM<-function(Amat ,cores=2){
  sbm <- BM_bernoulli("SBM_sym",Amat, plotting="", verbosity=0,ncores=cores,
                      explore_min=2, explore_max=20, exploration_factor=5)
  sbm$estimate()
  return(list(estim.sbm=sbm, best=which.max(sbm$ICL)))

}

#' Perform clustering of square matrices according to Stochastic Bloc Modelling
#'
#' @param freq.thresh Threshold of edge detection.
#' @param selected_taxo Taxonomic information about the MGS of the data set.
#' @param force.nb.groups Boolean to force the number of groups.
#' @param edge_freq Edge final frequency, typically the result of an inference aggregation.
#' @param cores Cores for parallel computing.
#' @param as.tibble Boolean to render the value \code{species_groups} as a tibble.
#'
#' @return A list containing
#'\itemize{
#'  \item{groups_features: }{number of groups of the best partition found by SBM.}
#'  \item{groups: }{vector of group membership.}
#'  \item{species_groups: }{list (or tibble if as.tibble=TRUE) of each group's taxonomic information (name of msp, species,
#'   genus, and phylum).}
#'  \item{graph: }{graphical output, the network with nodes coloured according to their group.}
#'}
#' @export
#' @importFrom EMtree draw_network
#' @importFrom tibble as_tibble
#'
#' @examples
#' data("liver")
#' abund<-liver$abundances
#' meta<-liver$metadata
#' taxo=liver$taxonomy
#' counts_from_table<-get_count_table(abund.table=abund, mgs=taxo$msp_name,
#' prev.min=0.9,verbatim=TRUE)$data
#' species<-liver$msp_set_90
#'
#' inf_90<-liver$infer_prev90
#' adapt<-adapt_mean_stability(inf_90, mean.stability=0.7,plot=TRUE)
#'
#' taxo <- liver$taxonomy
#' selected_taxo=dplyr::left_join(data.frame(msp_name=species),taxo, by="msp_name")
#'  norm2<-compute_aggreg_measures(adapt$freqs)$norm2
#'  clustering_output<-clustering_SBM(norm2, selected_taxo=selected_taxo)
clustering_SBM<-function(edge_freq, freq.thresh=0.7,selected_taxo,force.nb.groups=NULL,
                         cores=2, as.tibble=FALSE){
  G=ToSym(edge_freq)
  G[G<freq.thresh]=0
  param.sbm=estim_SBM(G, cores=cores)
  cluster_groups<-get_cluster_groups(param.sbm, force.groups=force.nb.groups, G)
  if(param.sbm$best==1 && is.null(force.nb.groups)){
    g_cluster<-draw_network(G,layout="stress",btw_rank =1, remove_isolated = FALSE,
                           pal_edges = "#31374f")$G
  }else{
    g_cluster<-coloured_graph(G, cluster_groups)
  }
  species_groups<-groups_features(selected_taxo,cluster_groups, plot=FALSE)
  if(as.tibble){
    species_groups<-do.call(rbind, lapply(seq_along(species_groups), function(group){
      res<-species_groups[[group]]
      res$group<-group
      return(res)
    })) %>% as_tibble()
  }
  return(list(best_nb_groups=param.sbm$best,  graph=g_cluster,groups=cluster_groups,
              species_groups=species_groups))
}



#' Extract count information in each group
#'
#' @param data Count data set.
#' @param selected_taxo  Taxonomic information about the MGS of the data set.
#' @param groups Vector of group membership.
#'
#' @return list of count information for each group.
#' @importFrom dplyr select group_by summarise mutate filter
#' @importFrom tibble as_tibble
#' @noRd
table_groups<-function(groups, data){

  count_groups<-data %>% t() %>%as_tibble() %>% mutate(group=groups)%>% mutate(msp_name=colnames(data))

  groups_count<-lapply(1:max(groups), function(num){
    count_groups %>% filter(group==num)
  })
  return(groups_count)
}


#' Graph with CORES colouring
#'
#' @param groups vector of group membership.
#' @param correlation partial correlation matrix.
#' @param species_groups list of each group's taxonomic information
#'
#' @return A list containing
#'\itemize{
#'  \item{graph_final: }{A graph with nodes coloured according to their group.}
#'  \item{graph_sub: }{A list containing separates graphs per group.}
#'}
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom ggraph ggraph create_layout geom_edge_link theme_graph
#' @importFrom tidygraph activate as_tbl_graph
#' @importFrom ggplot2 scale_size_manual scale_color_manual
#' @noRd
draw_graph_CORE <- function(groups, G, species_groups){

  Clustered <- groups
  Non <- which(Clustered==0)
  Clustered[Non] <- 0
  Clustered[-Non] <- 0.1
  Clustered <- as.character(Clustered)

  Cores <- groups+1
  Cores[Cores==1]="Non_classified"

  new_graph <- graph_from_adjacency_matrix(G, weighted = TRUE) %>%
    as_tbl_graph() %>%
    activate(edges) %>%
    filter(weight > 1e-5)%>%
    activate(nodes) %>%
    mutate(Cores=Cores, Clustered=Clustered)

  set.seed(1)
  layout <- create_layout(new_graph, layout = 'igraph', algorithm = 'fr')

  sign_palette <- c(low="#F0FFFF",high = "#0000FF")
  Core_palette=c(brewer.pal(length(unique(Cores))-1,"Paired"),"#CCCCCC")

  graph_final <- ggraph(layout) +
    geom_edge_link(aes(alpha=weight))+
    ggiraph::geom_point_interactive(
      aes(x = x,
          y = y,
          size=Clustered,
          color   = Cores)
    )+
    scale_size_manual(values = c(2,3) )+
    scale_color_manual(values = Core_palette, )+ theme_graph(background = "white")+
    theme(legend.position = "none")

  graph_by_CORE <- lapply(1:max(groups), function(num){
    Gsub<-G[groups==num, groups==num]
    Gsub[Gsub<1e-5] <- 0

    # Read the mean layout

    layout_sub <- layout[which(layout$name %in% species_groups[[num]]$msp_name),]

    ##- Draw the guild

    new_graph <- graph_from_adjacency_matrix(Gsub, weighted = TRUE) %>%
      as_tbl_graph() %>%
      activate(edges) %>%
      filter(weight > 1e-5)%>%
      activate(nodes) %>%
      mutate(Species=species_groups[[num]]$species)

    graph_sub <- ggraph(new_graph, x=layout_sub$x, y=layout_sub$y) +
      geom_edge_link(aes(alpha=weight*15))+
      ggiraph::geom_point_interactive(aes(x = x,
                                          y = y,
                                          color   = Species))+
      theme_graph(background = "white")+
      theme(legend.text = element_text(size = 14, face = "bold"),
            legend.title = element_text(size=16, face="bold"))
  })

  return(list(graph_final=graph_final, graph_sub =graph_by_CORE))
}



#' Perform clustering of square matrices according to CORE-clustering
#'
#' @param data Count data set.
#' @param freq.thresh Threshold of edge detection.
#' @param selected_taxo Taxonomic information about the MGS of the data set.
#' @param edge_freq Edge final frequency, typically the result of an inference aggregation.
#' @param min.thresh Minimal number of variables per cluster
#'
#' @return A list containing
#'\itemize{
#'  \item{groups: }{vector of group membership.}
#'  \item{species_groups: }{list of each group's taxonomic information (name of msp, species,
#'   genus, and phylum).}
#'  \item{graph: }{graphical output, the network with nodes coloured according to their group.}
#'  \item{correlation: }{partial correlation matrix}
#'  \item{table_groups: }{list of msp name per group.}
#'}
#' @export
#' @importFrom EMtree draw_network
#'
#' @examples
#' data("liver")
#' abund<-liver$abundances
#' meta<-liver$metadata
#' taxo=liver$taxonomy
#' counts_from_table<-get_count_table(abund.table=abund, mgs=taxo$msp_name,
#' prev.min=0.9,verbatim=TRUE)$data
#' species<-liver$msp_set_90
#'
#' inf_90<-liver$infer_prev90
#' adapt<-adapt_mean_stability(inf_90, mean.stability=0.7,plot=TRUE)
#'
#' taxo <- liver$taxonomy
#' selected_taxo= dplyr::left_join(data.frame(msp_name=species),taxo, by="msp_name")
#'  norm2<-compute_aggreg_measures(adapt$freqs)$norm2
#'  clustering_output<-clustering_CORE(counts_from_table,norm2, freq.thresh=0.9,
#'  selected_taxo=selected_taxo,min.thresh=4)
clustering_CORE<-function(data, edge_freq, freq.thresh=0.7,selected_taxo,
                          min.thresh){
  qdat <- SPRING::mclr(data)

Sigma <- mixedCCA::estimateR(qdat, type = "trunc", method = "approx",
                             tol = 1e-6, verbose = FALSE)$R

n<-nrow(data)

Parr_Corr<-get_par_corr(edge_freq,Sigma,n,freq.thresh)
correlation <-Parr_Corr$corp_ggm
G <- abs(ToSym(correlation))
dimnames(G)  <- list(rownames(Sigma), colnames(Sigma))

clusters <- principalFunction(G,min.thresh)

Cores <- clusters[[2]][,2]

species_groups <- groups_features(selected_taxo, Cores)
table_groups <- table_groups(Cores,data)

g_cluster=draw_graph_CORE(Cores, G, species_groups)

return(list(groups = Cores, species_groups=species_groups, graph = g_cluster, table_groups=table_groups))
}


