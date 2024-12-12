#' Simulate data from some empirical count dataset.
#'
#' @param real_data Empirical real dataset.
#' @param graph_type Structure type for the conditional dependency structure (among "cluster","erdos", "tree", and "scale-free").
#' @param graph Optional graph to be used. Must have rownames and colnames and reference all features from real_data
#' @param n Number of samples to simulate.
#' @param plot Boolean for displaying the simulated network.
#' @param r For cluster structure: controls the within/between ratio connection probability.
#' @param dens Graph density (for cluster graphs) or edges probability (for erdös-renyi graphs).
#' @param seed Seed number for data generation (rmvnorm).
#' @param must_connect Boolean to force the output graph to be connected.
#' @param k For cluster structure: number of groups
#' @param verbose Boolean controlling the verbosity.
#' @param signed Boolean for simulating both positive and negative partial correlations. Default is to FALSE, which implies only negative partial correlations.
#'
#' @return A list containing the simulated discrete counts, the corresponding true partial
#'  correlation matrix from the latent Gaussian layer of the model and the original graph
#'   structure that was used.
#' @export
#' @importFrom EMtree generator_graph generator_param draw_network
#' @importFrom SPRING synthData_from_ecdf
#' @importFrom igraph is_connected graph_from_adjacency_matrix
#' @importFrom stats cov2cor
#' @examples
#' data("liver")
#' real_data<-liver$abundances
#' meta=liver$metadata
#' taxo=liver$taxonomy
#' count_data<-get_count_table(abund.table=real_data, mgs=taxo$msp_name,  prev.min=0.9)$data
#' simulation<-new_synth_data(count_data, graph_type="cluster",n=50, plot=TRUE)
#' str(simulation, max.level=1)
#' ## Use prespecified matrix and use only part of the data
#' G <- diag(nrow = nrow(real_data)) ## 1990  x 1990 matrix
#' dimnames(G) <- list(taxo$msp_name , taxo$msp_name)
#' simulation<-new_synth_data(count_data, graph = G, n=50, plot=FALSE)
#' str(simulation, max.level=1)
new_synth_data<-function(real_data, graph_type="cluster", must_connect=TRUE, graph = NULL, n=300,
                         plot=TRUE,seed=10010, r=50, dens=4, k=3, verbose=TRUE, signed=FALSE){

  p=ncol(real_data)
  species <- colnames(real_data)
  if (is.null(graph)) {
    if(!must_connect){
      set.seed(seed)
      G<-as.matrix(EMtree::generator_graph(p=p, graph=graph_type, dens=dens/p,r,k))
    }else{
      i<-0; connect<-FALSE
      while(!connect){
        i<-i+1
        set.seed(i)
        G<-as.matrix(EMtree::generator_graph(p=p, graph=graph_type, dens=dens/p,r,k))
        graph<-igraph::graph_from_adjacency_matrix (G)
        connect<-igraph::is_connected(graph)
      }
    }
  } else {
    if (!all(species %in% colnames(graph)) ||  !all(species %in% rownames(graph))) stop("Some species in the abundance dataset do not appear in the provided graph.")
    G <- graph[species, species]
  }
  dimnames(G) <- list(species, species)
  faithful_param<-EMtree::generator_param(G=G, signed=signed)
  parcor<- -cov2cor(as.matrix(faithful_param$omega))

  if(plot){
    Graph_witout_name <- G
    colnames(Graph_witout_name)  = rownames(Graph_witout_name) = NULL
    g<-draw_network(Graph_witout_name, layout="nicely", btw_rank=1)$G
    print(g)
  }
  if(verbose) cat("Simulation from real data ecdf...")
  simu_counts <- simulate_from_ecdf(real_data, Sigma=faithful_param$sigma,
                                    n=n, seed = seed, verbose = verbose)

  simulated_data<-list(counts=simu_counts,
                       par.cor=parcor,
                       G=G)
}

#' Adapted form SPRING::synthData_from_ecdf but faster and more robust
#'
#' Generates synthetic count data based on empirical cumulative distribution (ecdf) of real count data
#'
#' @param real_data Matrix of real count data of size n by p.
#' @param Sigma Covariance structure of size p by p.
#' @param n Number of samples
#' @param seed Seed number for data generation.
#' @param verbose Logical value: if TRUE, iteration and index calculation for each step printed out.
#'
#' @return The vector from the upper triangular part of A.mat.
#' @importFrom stats cov2cor pnorm quantile
#' @noRd
simulate_from_ecdf <- function (real_data, Sigma, n, seed = 10010, verbose = FALSE)
{
  p <- ncol(real_data)
  # zratio <- colMeans(real_data == 0)
  # maxabund <- apply(comm, 2, max)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  mv_norm <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = cov2cor(Sigma))
  mv_unif <- pnorm(mv_norm)
  ## Matrix of simulated counts
  sim_counts <- matrix(0, n, p)
  ## For each species, apply the inverse of the empirical distribution
  ## function to the uniform r.v.
  for (j in 1:p) {
    ## quantile(x, ., type = 1) is the inverse of the empirical distribution function of x
    sim_counts[, j] <- quantile(x = real_data[, j], probs = mv_unif[, j], names = FALSE, type = 1)
  }
  return(sim_counts)
}

#' Manage methods' specificities to get frequency vectors
#'
#' @param lambda Penalty level.
#' @param indiv_output Method's output from the \code{all_inferences()} function.
#' @param method Name of the method.
#'
#' @return The frequency vectors
#' @importFrom dplyr select
#' @noRd
get_vec<-function(lambda,indiv_output,method){
  lambda.seq<-indiv_output$lambda_stab$lambda
  indiv_output<-indiv_output$output
  if(method%in%c("SpiecEasi","Magma")){
    index<-which(lambda.seq==lambda)

    vec_freqs=ToVec(indiv_output$merge[[index]])
  }else if(method=="SPRING"){
    index<-which(lambda.seq==lambda)
    vec_freqs=ToVec(indiv_output$stars$merge[[index]])
  }else if(method%in%c("EMtree")){
    vec=(indiv_output > lambda)
    vec_freqs=colMeans(1 * vec)
  }else if(method%in%c("PLNnetwork","gCoda","ZiLN")){
    vec_freqs=indiv_output %>% filter(abs(Penalty-lambda)<1e-4) %>%
      select(Prob) %>% pull()
  }
  return(vec_freqs)
}

#' Computes precision and recall values
#'
#' @param inf_output Complete ouput from the \code{all_inferences()} function.
#' @param G Original graph structure to retrieve.
#' @param edge_thresh Threshold on selection frequency for edge detection.
#' @param aim_stab Stability threshold.
#'
#' @return A tibble summarising performance (precision, recall, stability, number of edges) for all
#' penalty levels of all individual inferences in the provided \code{inf_output} object.
#'
#' @export
#' @import ggplot2
#' @importFrom dplyr n filter summarise pull select left_join
#' @importFrom tibble as_tibble
#' @examples
#' data("liver")
#' real_data<-liver$abundances
#' meta=liver$metadata
#' taxo=liver$taxonomy
#' count_data<-get_count_table(abund.table=real_data, mgs=taxo$msp_name, prev.min=0.9)$data
#' simulation<-new_synth_data(count_data, graph_type="cluster",n=50, plot=TRUE)
#' some_inferences<-all_inferences(simulation$counts , rep.num=4,cores=2)
#'
#' Eval_inference<-evalQuali(some_inferences,simulation$G)
#' str(Eval_inference)
#' Eval_inference%>% ggplot(aes(TPR, PPV, color=method))+geom_point()+geom_line()+theme_light()
evalQuali<-function(inf_output, G, edge_thresh=0.9, aim_stab=0.9){
  methods = names(inf_output)

  list_quali <- lapply(seq_along(inf_output), function(num_inf) {
    method <- methods[num_inf]
    lambda.seq <- inf_output[[num_inf]]$lambda_stab$lambda
    method_quali <- lapply(lambda.seq, function(lambda) {
      vec_freqs <- get_vec(lambda, indiv_output = inf_output[[num_inf]],
                           method = method)
      freq_data = data.frame(freqs = vec_freqs, here = ToVec(G))
      TPR = freq_data %>% filter(here == 1) %>% summarise(TPR = sum(freqs >
                                                                      edge_thresh)/sum(G)) %>% pull()
      PPV = freq_data %>% filter(freqs > edge_thresh) %>%
        summarise(PPV = sum(here)/n()) %>% pull()
      quali <- data.frame(lambda = lambda, TPR = TPR, PPV = PPV,
                          method = method)
      return(quali)
    })
    method_quali <- do.call(rbind, method_quali) %>% as_tibble()
    method_quali <- left_join(method_quali, inf_output[[num_inf]]$lambda_stab,
                              by = "lambda")
    return(method_quali)
  })
  list_quali <- do.call(rbind, list_quali) %>% as_tibble() %>%
    group_by(method) %>% mutate(stab_point=1*(abs(stability-aim_stab)==min(abs(stability-aim_stab))))


  return(list_quali)
}















