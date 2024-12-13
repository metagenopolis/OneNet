# Robust network inference

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#' From vector to symmetric matrix
#'
#' Makes a symmetric matrix from the vector made of its upper triangular part.
#'
#' @param A.vec Vector obtain with ToVec or the \code{upper.tri()} function.
#'
#' @return The symmetric matrix with A.vec as off-diagonal terms.
#' @noRd
ToSym <- function(A.vec){
  n <- (1+sqrt(1+8*length(A.vec)))/2
  A.mat <- matrix(0, n, n)
  A.mat[upper.tri(A.mat)] <- A.vec
  A.mat <- A.mat + t(A.mat)
  return(A.mat)
}

#' From symmetrix matrix to vector
#'
#' Makes a vector from the upper triangular part of a symmetric matrix.
#'
#' @param A.mat A symmetric matrix.
#'
#' @return The vector from the upper triangular part of A.mat.
#' @noRd
ToVec <- function(A.mat){
  return(suppressMessages(A.mat[upper.tri(A.mat)]))
}


#' Function b(n) from the StARS original paper
#'
#' @param n Sample size.
#'
#' @return A subsample ratio.
#' @noRd
bn<-function(n){
  if(n>144){
    res<-10*sqrt(n)/n
  }else{res<-0.8}
  return(res)
}

#' gCoda resampler
#'
#' Performs StARS with gCoda.
#'
#' @param data Count data set.
#' @param lambda.seq Optional vector of penalties, otherwise computed internally.
#' @param lambda.min.ratio Minimal ratio between largest and smallest values within the grid of lambda penalties. Required for internal computation of lambda.seq.
#' @param covar Optional covariate data set.
#' @param nlambda Length of lambda grid.
#' @param B Number of subsamples.
#' @param cores Number of cores for parallel computation
#' @param subsample.ratio  The subsampling ratio. The default value is 10*sqrt(n)/n when n>144 and 0.8 when n<=144, where n is the sample size.
#'
#' @return A data.frame containing the selection frequencies for each lambda penalty.
#' @importFrom parallel mclapply
#' @importFrom stats as.formula var
#' @noRd
resamp_gcoda<-function(data,lambda.seq=NULL,lambda.min.ratio =1e-2,covar=NULL,
                       nlambda=100, B=100, cores=1,subsample.ratio=NULL){
  p=ncol(data) ; n = nrow(data) ;
  if(is.null(subsample.ratio)) subsample.ratio<-bn(n)
  V = round(subsample.ratio* n)
  x=data
  #recompute lambda.seq
  if(is.null(lambda.seq)){
    x <- x + 0.5;
    x <- x / rowSums(x);
    clr_data<-log(x) - rowMeans(log(x))
    if(is.null(covar)){
      S <- var(clr_data)
    }else{
      string<-paste("x", paste(colnames(covar), collapse=" + "), sep=" ~ ")
      formula<-as.formula(string)
      X = as.matrix(stats::lm(formula, x=T, data=covar)$x)
      model<-stats::lm(clr_data~X)
      U<-model$residuals
      S <- var(U)
    }
    lambda.max <- max(max(S - diag(p)), -min(S - diag(p)));
    lambda.min <- lambda.min.ratio * lambda.max;
    lambda.seq <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  }
  #resampling
  obj=mclapply(1:B,function(b){
    res=tryCatch({
      set.seed(b)
      sample = sample(1:n, V, replace = F)
      data.sample = data[sample, ]
      if(!is.null(covar)){
        covar.samp=data.frame(covar[sample, ])
      }else{covar.samp=NULL}
      res=gcoda(data.sample, counts=T, covar=covar.samp, lambda.seq=lambda.seq, nlambda = nlambda)$path
    }, error=function(e){e}, finally={})
    return(res)
  }, mc.cores=cores)
  #keep only successful optimizations
  good.cores=which(do.call(rbind,lapply(obj,length))==nlambda)
  obj=obj[good.cores]
  #compute edge selection frequencies
  nb.lambda=min(do.call(rbind,lapply(obj, function(x) length(x))))
  if(!is.finite(nb.lambda)) message("Numerical instability in gCoda")
  list_freqs=do.call(rbind,lapply(1:nb.lambda, function(i){
    matfreq=Reduce("+",lapply(obj, function(x) x[[i]]))/B
    res=data.frame(freqs=ToVec(matfreq),
                   pen=lambda.seq[i])
    return(res)
  }))

  return(list_freqs)
}

#' Stability for EMtree
#'
#' Compute stability for each probability threshold on EMtree output.
#'
#' @param Pmat EMtree's Pmat output, from the \code{ResampleEMtree()} function.
#' @param n.levels Number of thresholds.
#' @param edge_thresh Threshold on edge selection frequency from which the edge is considered present. Required to compute the number of selected edges per probability threshold.
#'
#' @return A tibble containing the stability and number of selected edges per probability threshold.
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter
#' @noRd
stab_EMtree<-function(Pmat, n.levels=100, edge_thresh=0.9){
  vec_pr= exp(-seq(8, 1, length.out=n.levels) * log(10))
  list_lambda_stab<-lapply(vec_pr, function(lambda){
    freqs=colMeans(1 * (Pmat > lambda))
    res<-data.frame(lambda=lambda,
                    stability=1-4*mean(freqs*(1-freqs)),
                    nedges=sum(freqs>edge_thresh))
    return(res)
  })
  lambda_stab=do.call(rbind,list_lambda_stab) %>%as_tibble()
  lambda.min = lambda_stab$lambda[which.min(lambda_stab$stability)]
  lambda_stab=lambda_stab %>% filter(lambda>=lambda.min)
  lambda_stab=lambda_stab[,c("stability", "nedges","lambda")]
  return(lambda_stab)
}

#' Stability for penalized inference methods
#'
#' Computes stability for each penalty level.
#'
#' @param stars_output Output from a penalized inference method, containing the \code{merge} and \code{lambda} values.
#' @param edge_thresh Threshold on edge selection frequency from which the edge is considered present. Required to compute the number of selected edges per probability threshold.
#'
#' @return A tibble containing the stability and number of selected edges per penalty level.
#' @noRd
stab_stars<-function(stars_output,edge_thresh=0.9){
  lambda_seq<-stars_output$lambda
  list_lambda_stab<-lapply(1:length(lambda_seq), function(index){
    lambda=lambda_seq[index]#output$est$lambda
    vec_freqs=suppressMessages(ToVec(stars_output$merge[[index]]))
    lambda_stab=data.frame(stability=1-4*mean(vec_freqs*(1-vec_freqs)),
                           nedges=sum(vec_freqs>edge_thresh),lambda=lambda)
    return(lambda_stab)
  })
  lambda_stab=do.call(rbind,list_lambda_stab) %>% as_tibble()
  return(lambda_stab)
}

#' Wrapper for EMtree
#'
#' @param data Count data set.
#' @param rep.num Number of subsamples.
#' @param n.levels Size of probability grid.
#' @param edge_thresh Threshold on selection frequencies to set the number of selected edges.
#' @param Offset Vector of offsets.
#' @param covar Covariate data set.
#' @param spring Boolean: should the latent covariance matrix be estimated using the rank-based estimator from the SPRING method?
#' @param subsample.ratio  The subsampling ratio. The default value is 10*sqrt(n)/n when n>144 and 0.8 when n<=144, where n is the sample size.
#' @param ... Other options used to in the wrappers
#'
#' @return A list containing
#'  \itemize{
#'  \item{output: }{the matrix of all edges probabilities per subsample.}
#'  \item{lambda_stab: }{the tibble gathering the stability and number of detected edges per probability threshold.}
#'}
#' @export
#' @importFrom EMtree ResampleEMtree
#' @importFrom stats as.formula model.matrix
#' @examples
#' \dontrun{
#' data("liver")
#' abund<-liver$abundances
#' meta=liver$metadata
#' taxo=liver$taxonomy
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9, verbatim=FALSE)$data
#' EMtree_output<-just_EMtree(counts_from_table, rep.num=2)
#' str(EMtree_output, max.level=2)
#' }
just_EMtree<-function(data, rep.num=100, n.levels=100, edge_thresh=0.9, Offset=TRUE,
                      covar=NULL, spring=TRUE, subsample.ratio=NULL, ...){
  cat("\nEMtree...")
  t1<-Sys.time()
  p=ncol(data) ; ne = p*(p-1)/2
  if(Offset){ Offset= tryCatch(PLNmodels::compute_offset(data, offset = "GMPR"),
                               error=function(e){
                                 cat("GMPR failed")
                                 PLNmodels::compute_offset(data, offset = "RLE", pseudocounts=1)
                               },finally={})
  }else{Offset=rep(1, nrow(data))}
  if(is.null(subsample.ratio)) subsample.ratio<-bn(nrow(data))
  if (is.null(covar)) { X = matrix(1, nrow = nrow(data), ncol = 1)
  }else{
    string<-paste("~", paste(colnames(covar), collapse=" + "))
    formula<-as.formula(string)
    X = model.matrix(formula, data=covar) }#clean way to deal with categorical variables
  if(spring){
    userfunction<-function(counts, covar_matrix, sample){
      Kcor <- mixedCCA::estimateR(counts[sample,], type = "trunc", method = "approx", tol = 1e-6, verbose = FALSE)$R
    }
  }else{userfunction<-NULL}

  resamp_EMtree<-quiet(ResampleEMtree(data,S = rep.num, cores=1,maxIter = 100,O = Offset,covar_matrix =X,
                                      user_covariance_estimation=userfunction,v = subsample.ratio))
  lambda_stab=stab_EMtree(resamp_EMtree$Pmat,n.levels,edge_thresh)
  t2<-Sys.time() ; diff=difftime(t2, t1)
  cat(paste0(round(diff,3),attr(diff,"units")))
  return(list(output=resamp_EMtree$Pmat, lambda_stab=lambda_stab, nedges=ne))
}

#' Wrapper for SpiecEasi
#'
#' @inheritParams just_EMtree 
#' @param cores Number of cores for parallel computation
#'
#' @return A list containing
#'  \itemize{
#'  \item{output: }{output for SpiecEasi, with the merge, the sequence of lambda and the refit.}
#'  \item{lambda_stab: }{the tibble gathering the stability and number of detected edges per level of penalty.}
#'}
#' @export
#' @importFrom SpiecEasi spiec.easi clr sparseiCov
#' @importFrom pulsar getMaxCov
#' @importFrom huge huge
#' @importFrom stats as.formula lm cor
#' @examples
#' \dontrun{
#' data("liver")
#' abund<-liver$abundances
#' meta=liver$metadata
#' taxo=liver$taxonomy
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9, verbatim=FALSE)$data
#' SpiecEasi_output<-just_spiec(counts_from_table, rep.num=2, cores=1)
#' str(SpiecEasi_output, max.level=2)
#' }
just_spiec<-function(data, rep.num=100, n.levels=100, cores=1, edge_thresh=0.9,covar=NULL,subsample.ratio=NULL, ...){
  cat("\nSpiecEasi...")
  t1<-Sys.time()
  p=ncol(data) ; ne = p*(p-1)/2
  stability=0.9
  if(is.null(subsample.ratio)) subsample.ratio<-bn(nrow(data))
  if(!is.null(covar)){
    clr.matrix <- function(x.f, mar=2, ...) {apply(x.f, mar, clr, ...)}
    #clr transformation with pseudo-count
    U<-t(clr(data+1,mar=1)) #centered log-ratio data transformation
    string<-paste("U", paste(colnames(covar), collapse=" + "), sep=" ~ ")
    formula<-as.formula(string)
    X = as.matrix(lm(formula, x=T, data=covar)$x)
    data<-lm(U~X)$residuals   # spiecEasi is run on residuals of transform data,
    # to take covariates into account
    args<-list()
    args$lambda.max <- quiet(getMaxCov(cor((data))))
    args$lambda.min.ratio<-1e-2
    args$nlambda<-n.levels
    args$lambda <- quiet(getLamPath(args$lambda.max, args$lambda.max*args$lambda.min.ratio,
                                    args$nlambda, log=FALSE))
    args$lambda.min.ratio <- args$nlambda <- args$lambda.max <- NULL
    spiec.out=quiet(pulsar(data, criterion = "stars", fun= huge,fargs = args,
                           thresh = (1-stability),subsample.ratio=subsample.ratio,
                           rep.num = rep.num, ncores = cores))
    output=list(merge=spiec.out$stars$merge, lambda=args$lambda,
                opt.index=spiec.out$stars$opt.index,
                refit=spiec.out$refit)

  }else{

    spiec.out=spiec.easi(data, method="mb",pulsar.select = TRUE,nlambda=n.levels,
                                          sel.criterion="stars",
                                          lambda.min.ratio=1e-2,verbose = FALSE,
                                          pulsar.params=list(ncores=cores,rep.num=rep.num,
                                                             subsample.ratio=subsample.ratio,
                                                             thresh=(1-stability)))

    output=list(merge=spiec.out$select$stars$merge, lambda=spiec.out$lambda,
                opt.index=spiec.out$select$stars$opt.index,
                refit=spiec.out$refit)
  }

  lambda_stab=stab_stars(output,edge_thresh)
  t2<-Sys.time() ; diff=difftime(t2, t1)
  cat(paste0(round(diff,3),attr(diff,"units")))
  return(list(output=output, lambda_stab=lambda_stab, nedges=ne))
}

#' Wrapper for Magma
#'
#' @inheritParams just_EMtree 
#'
#' @return A list containing
#'  \itemize{
#'  \item{output: }{complete output from \code{rMAGMA::magma()}.}
#'  \item{lambda_stab: }{the tibble gathering the stability and number of detected edges per level of penalty.}
#'}
#' @export
#' @importFrom rMAGMA magma
#' @importFrom stats as.formula  model.matrix
#' @examples
#' \dontrun{
#' data("liver")
#' abund<-liver$abundances
#' meta=liver$metadata
#' taxo=liver$taxonomy
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9, verbatim=FALSE)$data
#' Magma_output<-just_magma(counts_from_table, rep.num=2)
#' str(Magma_output, max.level=2)
#' }
just_magma<-function(data, rep.num=100, n.levels=100, edge_thresh=0.9, Offset=TRUE,
                     covar=NULL, subsample.ratio=NULL,...){
  cat("\nMagma...")
  t1<-Sys.time()
  stability=0.9
  p=ncol(data) ; ne = p*(p-1)/2
  if(is.null(subsample.ratio)) subsample.ratio<- bn(nrow(data))
  if(Offset){
    seq_depth = "GMPR"
  }else{seq_depth = "unif"}
  if (is.null(covar)) { X = NULL#matrix(1, nrow = nrow(data), ncol = 1)
  }else{
    string<-paste("~", paste(colnames(covar), collapse=" + "))
    formula<-as.formula(string)
    X = model.matrix(formula, data=covar) }
  output=quiet(magma(data=data, distrib="ZINB",criterion.select="stars",seq_depth=seq_depth,
                     lambda.min.ratio=1e-2,method="mb",X = X,
                     stars.subsample.ratio=subsample.ratio,
                     nlambda=n.levels, rep.num=rep.num,stars.thresh=1-stability) )
  lambda_stab=stab_stars(output,edge_thresh)
  t2<-Sys.time() ; diff=difftime(t2, t1)
  cat( paste0(round(diff,3),attr(diff,"units")))
  output$path<-output$beta<-NULL
  return(list(output=output, lambda_stab=lambda_stab, nedges=ne))
}

#' Wrapper for SPRING
#'
#' @inheritParams just_EMtree 
#' @param seed A numeric seed for subsampling (native from the original SPRING package).
#' @param cores Number of cores for parallel computation
#'
#' @return A list containing
#'  \itemize{
#'  \item{output: }{output value from \code{SPRING::SPRING()}, result from  \code{pulsar::pulsar()} based on StARS criterion.}
#'  \item{lambda_stab: }{the tibble gathering the stability and number of detected edges per level of penalty.}
#'}
#' @export
#' @importFrom SPRING hugeKmb SPRING
#' @importFrom pulsar pulsar
#' @importFrom Rfit rfit
#' @importFrom stats as.formula
#' @examples
#' \dontrun{
#' data("liver")
#' abund<-liver$abundances
#' taxo=liver$taxonomy
#' meta=liver$metadata
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9, verbatim=FALSE)$data
#' spring_output<-just_spring(counts_from_table, rep.num=2, cores=1)
#' str(spring_output, max.level=2)
#' }
just_spring<-function(data, rep.num=100, n.levels=100, cores=1, edge_thresh=0.9,covar=NULL,seed=10010,
                      subsample.ratio=NULL,...){
  cat("\nSPRING...")
  t1<-Sys.time()
  stability=0.9 ; lambda.min.ratio=1e-2
  p=ncol(data) ; ne = p*(p-1)/2
  if(is.null(subsample.ratio)) subsample.ratio<-bn(nrow(data))
  if(!is.null(covar)){
    #modified clr transform
    qdat <- SPRING::mclr(data)
    #------add-in to account for covariates
    qdat_modif=qdat
    sapply(1:ncol(qdat),function(col_ind){
      non_nul=which(qdat[,col_ind]!=0)
      string<-paste("qdat_modif[,col_ind]", paste(colnames(covar), collapse=" + "), sep=" ~ ")
      formula<-as.formula(string)
      qdat_modif[non_nul,col_ind]<<-rfit(formula,data = covar)$residuals[non_nul]
    })
    qdat_modif[qdat_modif!=0]<-qdat_modif[qdat_modif!=0]+abs(min(qdat_modif))
    #------
    #fit sigma tilde
    Kcor <- mixedCCA::estimateR(qdat_modif, type = "trunc", method = "approx", tol = 1e-6, verbose = FALSE)$R
    # generate lambda sequence
    lambda.max <- max(max(Kcor-diag(p)), -min(Kcor-diag(p)))
    lambda.min <- lambda.min.ratio * lambda.max
    lambdaseq <- exp(seq(log(lambda.max), log(lambda.min), length = n.levels))
    #fit stars
    fun <- hugeKmb
    output <- quiet(pulsar(qdat, fun = fun, fargs = list(lambda = lambdaseq, Rmethod = "approx",
                                                                 tol = 1e-6, verbose = FALSE, verboseR = FALSE),
                                   rep.num = rep.num, criterion = 'stars',
                                   subsample.ratio=subsample.ratio, ncores = cores, thresh = 1-stability,
                                   seed=seed))
  }else{
    output <-quiet(SPRING(data, Rmethod = "approx", quantitative = FALSE,
                          lambdaseq = "data-specific", nlambda = n.levels,seed=seed, rep.num = rep.num,
                          ncores=cores, thresh=1-stability,subsample.ratio=subsample.ratio)$output)
  }
  stars_output<-list(merge=output$stars$merge,lambda=output$est$lambda)
  lambda_stab=stab_stars(stars_output,edge_thresh)
  t2<-Sys.time() ; diff=difftime(t2, t1)
  cat(paste0(round(diff,3),attr(diff,"units")))
  output$est<-NULL
  return(list(output=output, lambda_stab=lambda_stab, nedges=ne))
}


#' Resampler for zero-inflation log_Normal
#'
#' @param data count data set
#' @param rep.num Number of resamples
#' @param n.levels Size of the penalty grid.
#' @param lambda.min.ratio Ratio between the maximal and the minimal penalties.
#' @param cores Number of cores.
#' @param subsample.ratio  The subsampling ratio. The default value is 10*sqrt(n)/n when n>144 and 0.8 when n<=144, where n is the sample size.
#'
#' @return A table gathering all edges selection frequency for each penalty level.
#' @importFrom pulsar getLamPath
#' @importFrom huge huge.mb
#' @importFrom stats cor
##' @noRd
resamp_ZiLN<-function(data, rep.num, n.levels,lambda.min.ratio =1e-2, cores,subsample.ratio=NULL){
  p<-ncol(data) ; n<-nrow(data)
  if(is.null(subsample.ratio)) subsample.ratio<-bn(n)
  V = round(subsample.ratio * n)
  #set lambda sequence for all resamples
  Z = infer_Z(data, seq_depth = "TS")
  S = cor(Z)
  lambda_max = max(S[upper.tri(S)])
  lambda_min<-lambda.min.ratio*lambda_max
  pen = getLamPath (min=lambda_min,max = lambda_max, len=n.levels)
  #compute all resamples
  list_paths<- mclapply(1:rep.num, function(b){
    set.seed(b)
    sample = sample(1:n, V, replace = F)
    data.sample = data[sample, ]
    Z = infer_Z(data.sample, seq_depth = "TS")
    mb_Z2 = huge.mb(Z, lambda = pen, verbose = F)
    path = mb_Z2$path

  }, mc.cores=cores)

  list_freqs=do.call(rbind,lapply(1:n.levels, function(i){
    matfreq=Reduce("+",lapply(list_paths, function(x) x[[i]]))/rep.num
    res=data.frame(freqs=ToVec(matfreq),
                   pen=pen[i])
    return(res)
  }))
  return(list_freqs)
}

#' Wrapper for ZiLN
#'
#' @inheritParams just_EMtree 
#' @param cores Number of cores for parallel computation
#'
#' @return A list containing
#'  \itemize{
#'  \item{output: }{merge and the lambda grid.}
#'  \item{lambda_stab: }{the tibble gathering the stability and number of detected edges per level of penalty.}
#'}
#' @export
#' @importFrom dplyr rename group_by summarise pull select
#' @importFrom tibble as_tibble
#' @examples
#' \dontrun{
#' data("liver")
#' abund<-liver$abundances
#' taxo=liver$taxonomy
#' meta=liver$metadata
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9, verbatim=FALSE)$data
#' ZiLN_output<-just_ZiLN(counts_from_table, rep.num=2, cores=1)
#' str(ZiLN_output, max.level=2)
#' }
just_ZiLN<-function(data, rep.num=100,n.levels=100, cores=1,edge_thresh=0.9, subsample.ratio=NULL,... ){
  cat("\nZiLN...")
  t1<-Sys.time()
  p=ncol(data) ; ne = p*(p-1)/2
  list_freqs<-resamp_ZiLN(data,lambda.min.ratio =1e-2,n.levels=n.levels, rep.num=rep.num, cores=cores,
                          subsample.ratio=subsample.ratio)
  lambda.seq= unique(list_freqs%>% select(pen) %>% pull())

  lambda_stab=list_freqs %>% rename(lambda=pen) %>%  group_by(lambda) %>%
    summarise(stability=1-4*mean(freqs*(1-freqs)), nedges=sum(freqs>edge_thresh))
  lambda_stab=lambda_stab[,c("stability", "nedges","lambda")]
  t2<-Sys.time() ; diff=difftime(t2, t1)
  cat( paste0(round(diff,3),attr(diff,"units")))
  list_freqs<-list_freqs %>% as_tibble() %>%
    rename(Prob=freqs, Penalty=pen)
  return(list(output=list_freqs, lambda_stab=lambda_stab, nedges=ne))
}

#' Wrapper for PLNnetwork
#'
#' @inheritParams just_EMtree 
#'
#' @return A list containing
#'  \itemize{
#'  \item{lambda_stab: }{the tibble gathering the stability and number of detected edges per level of penalty.}
#'  \item{lambda.seq: }{the computed lambda sequence.}
#'}
#' @export
#' @importFrom PLNmodels PLNnetwork getBestModel PLNnetwork_param  stability_selection
#' @importFrom dplyr rename group_by summarise select
#' @importFrom tibble as_tibble
#' @importFrom stats model.matrix as.formula
#' @importFrom future plan
#' @examples
#' \dontrun{
#' data("liver")
#' abund<-liver$abundances
#' taxo=liver$taxonomy
#' meta=liver$metadata
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9, verbatim=FALSE)$data
#' PLN_output<-just_PLN(counts_from_table, rep.num=2)
#' str(PLN_output, max.level=2)
#' }
just_PLN<-function(data, rep.num=100, n.levels=100,edge_thresh=0.9, Offset=TRUE, covar=NULL,subsample.ratio=NULL, ...){
  cat("\nPLNnetwork...")
  stability=0.9
  t1<-Sys.time()
  p=ncol(data) ; ne = p*(p-1)/2 
  n = nrow(data)
  if(Offset){ Offset= tryCatch(PLNmodels::compute_offset(data, offset = "GMPR"),
                               error=function(e){
                                 cat("GMPR failed")
                                 PLNmodels::compute_offset(data, offset = "RLE", pseudocounts=1)
                               },finally={})
  }else{Offset=rep(1, n)}
  if(is.null(subsample.ratio)) subsample.ratio<-bn(n)
  if (is.null(covar)) { X = matrix(1, nrow = n, ncol = 1)
  }else{
    string<-paste("~", paste(colnames(covar), collapse=" + "))
    formula<-as.formula(string)
    X = model.matrix(formula, data=covar) }
  
  future::plan(future::multisession, workers = 4)
  network_models<-quiet(PLNnetwork(data~-1+ offset(log(Offset))+.,data=data.frame(X),
                                   control=PLNmodels::PLNnetwork_param(min_ratio=1e-2,n_penalties=n.levels)))

  subs <- replicate(rep.num, sample.int(n, size = subsample.ratio* n), simplify = FALSE)
  PLNmodels::stability_selection(network_models, subsamples = subs)
  model_StARS <- quiet(getBestModel(network_models, "StARS"))
  future::plan("sequential")
  lambda_stab=network_models$stability_path %>% as_tibble() %>% rename(lambda=Penalty) %>%
    group_by(lambda) %>%
    summarise(nedges=sum(Prob>edge_thresh), stability=1-4*mean(Prob*(1-Prob)))
  lambda_stab=lambda_stab[,c("stability", "nedges","lambda")]
  # extract lambda sequence for gCoda
  lambda.seq=network_models$stability_path %>% as_tibble() %>% select(Penalty) %>% pull()
  lambda.seq=sort(unique(lambda.seq), decreasing = TRUE)
  t2<-Sys.time() ; diff=difftime(t2, t1)
  cat(paste0(round(diff,3),attr(diff,"units")))
  output<-network_models$stability_path %>% select(Penalty,Prob)
  return(list(output=output,lambda_stab=lambda_stab, nedges=ne))
}

#' Wrapper for gCoda
#'
#' @inheritParams just_EMtree 
#' @param lambda.seq Grid of lambda penalties.
#' @param cores Number of cores for parallel computation
#'
#' @return A list containing
#'  \itemize{
#'  \item{output: }{merge and the lambda grid.}
#'  \item{lambda_stab: }{the tibble gathering the stability and number of detected edges per level of penalty.}
#'}
#' @export
#' @importFrom dplyr rename group_by summarise pull select
#' @importFrom tibble as_tibble
#' @examples
#' \dontrun{
#' data("liver")
#' abund<-liver$abundances
#' taxo=liver$taxonomy
#' meta=liver$metadata
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9, verbatim=FALSE)$data
#' gcoda_output<-just_gcoda(counts_from_table, rep.num=2, cores=1)
#' str(gcoda_output, max.level=2)
#' }
just_gcoda<-function(data, rep.num=100,n.levels=100, cores=1,edge_thresh=0.9,lambda.seq=NULL, covar=NULL,
                     subsample.ratio=NULL, ...){
  cat("\ngCoda...")
  t1<-Sys.time()
  p=ncol(data) ; ne = p*(p-1)/2
  list_freqs<-quiet(resamp_gcoda(data,lambda.seq=lambda.seq,covar=covar,lambda.min.ratio =1e-2,
                                 nlambda=n.levels, B=rep.num, cores=cores,subsample.ratio=subsample.ratio))
  lambda.seq= unique(list_freqs %>% select(pen) %>% pull())

  lambda_stab=list_freqs %>% rename(lambda=pen) %>%  group_by(lambda) %>%
    summarise(stability=1-4*mean(freqs*(1-freqs)), nedges=sum(freqs>edge_thresh))
  lambda_stab=lambda_stab[,c("stability", "nedges","lambda")]
  t2<-Sys.time() ; diff=difftime(t2, t1)
  cat( paste0(round(diff,3),attr(diff,"units")))
  list_freqs<-list_freqs %>% as_tibble() %>%
    rename(Prob=freqs, Penalty=pen)
  return(list(output=list_freqs, lambda_stab=lambda_stab, nedges=ne))
}



#' Aggregation
#'
#' Compute several aggregations: mean, norme 2, Inverse Variance Weighting (IVW), and number of methods detecting the edge (nb.high.freq).
#'
#' @param data A table of size (number of edges)x(number of methods), filled with edges selection frequencies.
#'
#' @return A data.frame with the different aggregations in column.
#' @export
#'
#' @examples
#' data("liver")
#' inf_90<-liver$infer_prev90
#' adapt<-adapt_mean_stability(inf_90, mean.stability=0.8,plot=TRUE)
#' aggreg<-compute_aggreg_measures(adapt$freqs)
#' head(aggreg)

compute_aggreg_measures<-function(data){
  mean<-rowMeans(data)
  norm2<-apply(data, 1, function(x){
    sqrt(sum(x^2))/(sqrt(sum(rep(1,ncol(data))^2)))})
  IVW<-apply(data, 1, function(x){
    v = x * (1 - x)
    vratio = 1/(v + (v == 0))
    if(min(x)<0.01 || max(x)>0.99){
      indices<-which(x<0.01 | x>0.99)
      vratio[indices] <- 105
    }
    return(sum(x * vratio)/sum(vratio))
  })
  aggreg_measures=data.frame(mean=mean,norm2=norm2,IVW=IVW,
                             nb.high.freq= rowSums(data>0.9))
  return(aggreg_measures)
}



#' Wrapper for partial correlations computation
#'
#' @param edge_freq Edge final frequency, typically the result of an inference aggregation.
#' @param Sigma Estimation of the empirical sample covariance matrix.
#' @param n Number of samples.
#' @param freq.thresh Threshold of edge detection.
#'
#' @return Partial correlations between all variables, in vector format.
#' @export
#' @importFrom ggm fitConGraph parcor
#' @examples
#' data("liver")
#' abund<-liver$abundances
#' meta<-liver$metadata
#' taxo=liver$taxonomy
#' meta=liver$metadata
#' counts_from_table<-get_count_table(abund.table=abund, prev.min=0.9, verbatim=FALSE)$data
#' n<-nrow(counts_from_table)
#'
#' #first estimate Sigma thanks to SPRING model
#' qdat <- SPRING::mclr(counts_from_table)
#' Sigma <- mixedCCA::estimateR(qdat, type = "trunc", method = "approx",tol = 1e-6,
#'  verbose = FALSE)$R
#'
#' # get the norme 2 aggregation from the inference using all methods, with mean stbaility 0.7
#' data("liver")
#' inf_90<-liver$infer_prev90
#' adapt<-adapt_mean_stability(inf_90, mean.stability=0.8,plot=TRUE)
#' mean<-compute_aggreg_measures(adapt$freqs)$mean
#'
#' # compute partial correlations
#'  ParCor<-get_par_corr(mean,Sigma,n,freq.thresh=0.9)
get_par_corr<-function(edge_freq,Sigma,n,freq.thresh=0.9){
  Ghat_pstab=ToSym(1*edge_freq>freq.thresh)
  colnames(Ghat_pstab)=rownames(Ghat_pstab)=colnames(Sigma)=rownames(Sigma)=1:ncol(Sigma)

  #compute partial corr
  fitggm=fitConGraph(amat = Ghat_pstab,S =Sigma,n =n)
  corp_ggm=ToVec(parcor(S = fitggm$Shat))
  res= data.frame(corp_ggm)
  return(res)
}

#' Call methods
#'
#' @param M character indicating with which method the inference should be
#'  performed 
#' @return An output from the desired individual method.
#' @noRd
apply_method=function(M,...){
  call_function=switch(M,
                       PLNnetwork=just_PLN,
                       gCoda = just_gcoda,
                       SpiecEasi= just_spiec,
                       Magma=just_magma,
                       SPRING = just_spring,
                       EMtree = just_EMtree,
                       ZiLN = just_ZiLN)
  do.call(call_function, list( ...))
}
