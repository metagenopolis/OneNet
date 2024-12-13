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


#' Graph with group colouring
#'
#' Wrapper of EMtree::draw_network, which uses ggraph and tidygraph.
#'
#' @param G Adjacency matrix, possibly weighted.
#' @param groups Groups on nodes, possibly the result from a clustering.
#' @param palette vector of colours defining the cluster palette.
#'
#' @return A graph with nodes coloured according to their group.
#' @importFrom EMtree draw_network
#' @importFrom grDevices rainbow
#' @noRd
coloured_graph<-function(G, groups, palette){
 
  new_graph <- graph_from_adjacency_matrix(G, weighted = TRUE) %>%
    as_tbl_graph() %>%
    activate(edges) %>%
    filter(weight > 1e-5)%>%
    activate(nodes) %>%
    mutate(Clusters=as.character(groups))
  
  set.seed(1)
  layout <- create_layout(new_graph, layout = 'igraph', algorithm = 'fr')
  
  sign_palette <- c(low="#F0FFFF",high = "#0000FF")
  Clusters_palette=palette[1:length(unique(groups))]
  
  g <- ggraph(layout) +
    geom_edge_link(aes(alpha=weight),  show.legend = F)+
    ggiraph::geom_point_interactive(
      aes(x = x,
          y = y,
          color=Clusters)
    )+
    scale_color_manual(values = Clusters_palette, )+ theme_graph(background = "white")
  
  return(g)
}

#' Extract spcies information in each group
#'
#' @param selected_taxo  Taxonomic information about the MGS of the data set.
#' @param groups Vector of group membership.
#'
#' @return list of taxonomic information for each group.
#' @importFrom dplyr select group_by summarise mutate filter
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal labs
#' @importFrom tibble as_tibble
#' @noRd
groups_features<-function(selected_taxo,  groups){
  taxo_groups<-selected_taxo %>%as_tibble() %>% mutate(group=groups)
  #groups=vecteur des groupes des espèces
  groups_species<-lapply(1:max(groups), function(num){
    taxo_groups %>% filter(group==num) 
  })
  return(groups_species)
}


#' Perform clustering of square matrices according to Stochastic Bloc Modelling
#'
#' @param freq.thresh Threshold of edge detection.
#' @param selected_taxo Taxonomic information about the MGS of the data set.
#' @param edge_freq Edge final frequency, typically the result of an inference aggregation.
#' @param palette vector of colours defining the cluster palette.
#'
#' @return A list containing
#'\itemize{
#'  \item{groups: }{vector of group membership.}
#'  \item{species_groups: }{list (or tibble if as.tibble=TRUE) of each group's taxonomic information (name of msp, species,
#'   genus, and phylum).}
#'  \item{graph: }{graphical output, the network with nodes coloured according to their group.}
#'}
#' @export
#' @importFrom EMtree draw_network
#' @importFrom sbm estimateSimpleSBM
#'
#' @examples
#' data("liver")
#' abund<-liver$abundances
#' meta<-liver$metadata
#' taxo=liver$taxonomy %>% dplyr::select(msp_name, species)
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9,verbatim=TRUE)$data
#' species<-liver$msp_set_90
#' inf_90<-liver$infer_prev90
#' adapt<-adapt_mean_stability(inf_90, mean.stability=0.8,plot=TRUE)
#' selected_taxo=dplyr::left_join(data.frame(msp_name=species), taxo, by="msp_name")
#' mean<-compute_aggreg_measures(adapt$freqs)$mean
#' clustering_output<-clustering_SBM(mean, selected_taxo=selected_taxo)
clustering_SBM<-function(edge_freq, freq.thresh=0.9, selected_taxo, palette = palette_default()){
  G=ToSym(edge_freq)
  G[G<freq.thresh]=0
  param.sbm= G %>% sbm::estimateSimpleSBM("bernoulli", dimLabels = 'tree', estimOptions = list(verbosity = 0, plot = FALSE))
  g_cluster<-coloured_graph(G, param.sbm$memberships, palette)

  species_groups<-groups_features(selected_taxo,param.sbm$memberships)

  return(list(graph=g_cluster,groups=param.sbm$memberships,
              species_groups=species_groups))
}

#' Extract vector of colors
#'
#' @return vector of colors.
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @noRd
palette_default <- function(){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  return(col_vector)
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
#' @param palette vector of colours defining the cluster palette
#'
#' @return A list containing
#'\itemize{
#'  \item{graph_final: }{A graph with nodes coloured according to their group.}
#'  \item{graph_sub: }{A list containing separates graphs per group.}
#'}
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom ggraph ggraph create_layout geom_edge_link theme_graph geom_node_point
#' @importFrom tidygraph activate as_tbl_graph
#' @importFrom ggplot2 scale_size_manual scale_color_manual
#' @importFrom randomcoloR distinctColorPalette
#' @noRd
draw_graph_CORE <- function(groups, G, species_groups, palette){

  Clustered <- groups
  Non <- which(Clustered==0)
  Clustered[Non] <- 0
  Clustered[-Non] <- 0.1
  Clustered <- as.character(Clustered)

  Cores <- groups
  Cores[Cores==0]="Non_classified"

  new_graph <- graph_from_adjacency_matrix(G, weighted = TRUE) %>%
    as_tbl_graph() %>%
    activate(edges) %>%
    filter(weight > 1e-5)%>%
    activate(nodes) %>%
    mutate(Cores=Cores, Clustered=Clustered)

  set.seed(1)
  layout <- create_layout(new_graph, layout = 'igraph', algorithm = 'fr')

  sign_palette <- c(low="#F0FFFF",high = "#0000FF")
  Core_palette= c(palette[1:(length(unique(Cores))-1)],"#CCCCCC")
  
  graph_final <- ggraph(layout) +
    geom_edge_link(aes(alpha=weight), show.legend = F)+
    ggiraph::geom_point_interactive(
      aes(x = x,
          y = y,
          size=Clustered,
          color   = Cores)
    )+
    scale_size_manual(values = c(2,3) )+
    scale_color_manual(values = Core_palette, )+ 
    theme_graph(background = "white")+ 
    guides(size="none")

  
  graph_by_CORE <- lapply(1:max(groups), function(num){
    Gsub<-G[groups==num, groups==num]
    Gsub[Gsub<1e-5] <- 0

    # Read the mean layout

    layout_sub <- layout[which(layout$name %in% species_groups[[num]]$msp_name),]

    ##- Draw the guild
    
    taxa_palette <- distinctColorPalette(50)
    Taxa <-substr(species_groups[[num]][,2,drop=TRUE],1,30)
    
    new_graph <- graph_from_adjacency_matrix(Gsub, weighted = TRUE) %>%
      as_tbl_graph() %>%
      activate(edges) %>%
      filter(weight > 1e-5)%>%
      activate(nodes) %>%
      mutate(Taxa=Taxa)

    graph_sub <- ggraph(new_graph, x=layout_sub$x, y=layout_sub$y) +
      geom_edge_link(aes(alpha=weight*15), show.legend = F)+
      ggiraph::geom_point_interactive(aes(x = x,
                                          y = y,
                                          color   = Taxa,
                                          size=rep("0.1",length(Taxa))))+
      theme_graph(background = "white")+
      scale_color_manual(values = taxa_palette,)+
      scale_size_manual(values=3)+
      theme(legend.text = element_text(size = 14, face = "bold"),
            legend.title = element_text(size=14, face="bold"))+ 
      guides(size="none",color = guide_legend(override.aes = list(size=3)))
  })

  return(list(graph_final=graph_final, graph_sub =graph_by_CORE))
}



#' Perform clustering of square matrices according to CORE-clustering
#'
#' @param data Count data set.
#' @param freq.thresh Threshold of edge detection.
#' @param selected_taxo Taxonomic information about the MGS of the data set.
#' @param edge_freq Edge final frequency, typically the result of an inference aggregation.
#' @param min.thresh Minimal number of variables per cluster.
#' @param palette vector of colours defining the cluster palette.
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
#' taxo=liver$taxonomy %>% dplyr::select(msp_name, species)
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9,verbatim=TRUE)$data
#' species<-liver$msp_set_90
#' inf_90<-liver$infer_prev90
#' adapt<-adapt_mean_stability(inf_90, mean.stability=0.8,plot=TRUE)
#' selected_taxo= dplyr::left_join(data.frame(msp_name=species),taxo, by="msp_name")
#'  mean<-compute_aggreg_measures(adapt$freqs)$mean
#'  clustering_output<-clustering_CORE(counts_from_table, mean, freq.thresh=0.9,
#'  selected_taxo=selected_taxo, min.thresh=3)
clustering_CORE<-function(data, edge_freq, freq.thresh=0.9,selected_taxo,
                          min.thresh, palette= palette_default()){
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

g_cluster=draw_graph_CORE(Cores, G, species_groups, palette)

return(list(groups = Cores, species_groups=species_groups, graph = g_cluster, table_groups=table_groups))
}


