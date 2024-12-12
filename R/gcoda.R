#-----------
# original gCoda code: https://github.com/huayingfang/gCoda

#' Adapted gCoda function
#'
#' @param x Dataset, either counts or proportions
#' @param counts Boolean for discrete nature of data
#' @param pseudo Value of the pseudo-count
#' @param lambda.min.ratio Ratio between the maximal and minimal value of the lambda sequence
#' @param nlambda Size of the lambda sequence
#' @param ebic.gamma Parameter for model selection with EBIC
#' @param covar Optional covariate data.frame
#' @param lambda.seq Optional pre-computed lambda sequence
#'
#' @return The gCoda fit object
#' @noRd
gcoda <- function(x, counts = T, pseudo = 0.5, lambda.min.ratio = 1e-2,
                  nlambda = 30, ebic.gamma = 0.5, covar=NULL,lambda.seq=NULL) {
  if(counts) {
    x <- x + pseudo;
    x <- x / rowSums(x);
  }
  n <- nrow(x);
  p <- ncol(x);
  # Log transformation for compositional data
  # modification for an easy fix to take covariates into account
  if(is.null(covar)){
    S <- var(log(x) - rowMeans(log(x)))
  }else{
    x<-log(x) - rowMeans(log(x))
    string<-paste("x", paste(colnames(covar), collapse=" + "), sep=" ~ ")
    formula<-as.formula(string)
    X = as.matrix(lm(formula, x=T, data=covar)$x)
    model<-lm(x~X)
    U<-model$residuals
    S <- var(U)
  }
  # Generate lambda via lambda.min.ratio and nlambda
  if(is.null(lambda.seq)){
    lambda.max <- max(max(S - diag(p)), -min(S - diag(p)));
    lambda.min <- lambda.min.ratio * lambda.max;
    lambda <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  }else{lambda=lambda.seq
  nlambda=length(lambda.seq)}
  # Store fit result for gcoda via a series of lambda
  fit <- list();
  fit$lambda <- lambda;
  fit$nloglik <- rep(0, nlambda);
  fit$df <- rep(0, nlambda);
  fit$path <- list();
  fit$icov <- list();
  # Compute solution paths for gcoda
  icov <- diag(p);
  for(i in 1:nlambda) {
    out.gcoda <- gcoda_sub(A = S, iSig = icov, lambda = lambda[i]);
    icov <- out.gcoda$iSig;
    fit$nloglik[i] <- out.gcoda$nloglik;
    fit$icov[[i]] <- icov;
    fit$path[[i]] <- 0 + (abs(icov) > 1e-15);
    diag(fit$path[[i]]) <- 0;
    fit$df[i] <- sum(fit$path[[i]]) / 2;
  }
  # Here I deleted the "Continue if edge density is too small"
  # because of numeric instabilities and inifinite running times
  # Compute EBIC score for lambda selection
  fit$ebic.score <- n * fit$nloglik + log(n) * fit$df +
    4 * ebic.gamma * log(p) * fit$df;
  fit$opt.index <- which.min(fit$ebic.score);
  fit$refit <- fit$path[[fit$opt.index]];
  fit$opt.icov <- fit$icov[[fit$opt.index]];
  fit$opt.lambda <- fit$lambda[fit$opt.index];
  return(fit);
}

#
#'  Optimization for gcoda with given lambda
#'
#' @param A Covariance matrix.
#' @param iSig Inverse covariance matrix.
#' @param lambda Penalty.
#' @param tol_err Tolerance for error.
#' @param k_max Maximal number of iterations.
#' @param verbatim Boolean controlling verbosity.
#' @noRd
gcoda_sub <- function(A, iSig = NULL, lambda = 0.1, tol_err = 1e-4,
                      k_max = 500, verbatim=FALSE) {
  p <- ncol(A);
  if(is.null(iSig)) iSig <- diag(p)
  err <- 1;
  k <- 0;
  fval_cur <- Inf;
  while(err > tol_err && k < k_max) {
    iSig_O <- rowSums(iSig);
    iS_iSig <- 1 / sum(iSig_O);
    iSig_O2 <- iSig_O * iS_iSig;
    A_iSig_O2 <- rowSums(A * rep(iSig_O2, each = p));
    A2 <- A - A_iSig_O2 - rep(A_iSig_O2, each = p) +
      sum(iSig_O2 * A_iSig_O2) + iS_iSig;
    iSig2 <- huge_glasso_mod(S = A2, lambda = lambda);
    fval_new <- obj_gcoda(iSig = iSig2, A = A, lambda = lambda);
    xerr <- max(abs(iSig2 - iSig) / (abs(iSig2) + 1));
    err <- min(xerr, abs(fval_cur - fval_new)/(abs(fval_new) + 1));
    if(is.na(err)){
      if(verbatim)  cat("BREAK")
      break
    }
    k <- k + 1;
    iSig <- iSig2;
    fval_cur <- fval_new;
  }
  nloglik <- fval_cur - lambda * sum(abs(iSig));
  if(k >= k_max) {
    cat("WARNING of gcoda_sub:\n", "\tMaximum Iteration:", k,
        "&& Relative error:", err, "!\n");
  }
  return(list(iSig = iSig, nloglik = nloglik));
}
#----------------------------------------
# Objective function value of gcoda (negative log likelihood + penalty)
#' Title
#'
#' @param iSig Inverse Covariance.
#' @param A Covariance.
#' @param lambda Penalty.
#' @noRd
obj_gcoda <- function(iSig, A, lambda) {
  p <- ncol(A);
  iSig_O <- rowSums(iSig);
  S_iSig <- sum(iSig_O);
  nloglik <- - log(det(iSig)) + sum(iSig * A) + log(S_iSig) -
    sum(iSig_O * rowSums(A * rep(iSig_O, each = p))) / S_iSig;
  pen <- lambda * sum(abs(iSig));
  return(nloglik + pen);
}
#----------------------------------------
#
# Input S must be covariance matrix
# library(huge)
# library(Rcpp)
#sourceCpp("~/Projet_guildes/source/hugeglasso.cpp")
#' Modified huge::huge.glasso for quick preparation
#'
#' @param S Covariance matrix
#' @param lambda Penalty
#' @noRd
huge_glasso_mod <- function(S, lambda) {
  icov <- diag(1/(diag(S) + lambda));
  z <- which(rowSums(abs(S) > lambda) > 1);
  q <- length(z);
  if (q > 0) {
    out.glasso= hugeglasso_sub(S = as.matrix((S[z, z])), W = as.matrix((S[z, z])), T = as.matrix(diag(as.double(q))),
                               d= as.integer(q), ilambda = as.double(lambda), df = as.integer(0),
                               scr=TRUE)
    icov[z, z] = matrix(out.glasso, ncol = q);
  }
  return(icov);
}

