#' Bayesian alignment
#'
#' \code{balign} implements the Bayesian alignment optimization
#'
#' @param data The dataset
#' @param n_items Number of items
#' @param n_ninv Number of noninvariant items
#' @param n_gps Number of groups
#' @param svp Small variance priors
#'
#' @return A vector of length 3.
#'          \item{lambda}{aligned factor loadings estimates}
#'          \item{n_gps}{number of groups}
#'          \item{nu}{aligned intercepts estimates}
#'

balign <- function(data, n_items, n_ninv, n_gps = 2, svp = FALSE){
  bcfa_res <- fit_true(data, n_items, n_ninv, n_gps = 2, svp)
  ## fit the bayesian configural model
  bcfa_res <- fit_config(data, n_items, group = "group", svp = svp)
  blava_draws <- blavInspect(bcfa_res, "draws")
  blava_draw_sam <- map_df(blava_draws, as.data.frame)
  ## reformat the parameter estimates to be used for alignment
  ly_groups <- blava_draw_sam[ , str_which(names(blava_draw_sam), "ly")[1:(n_items*2)]]
  lambda_mat <- lapply(as.list(1:dim(ly_groups)[1]), function(x) ly_groups[x[1],])
  lambda_mat_new <- lapply(lambda_mat, function(x){matrix(sapply(x, as.numeric), ncol = n_items, nrow = n_gps, byrow = TRUE)})
  nu_groups <- blava_draw_sam[ , str_which(names(blava_draw_sam), "Nu")[1:(n_items*2)]]
  nu_mat <- lapply(as.list(1:dim(nu_groups)[1]), function(x) nu_groups[x[1],])
  nu_mat_new <- lapply(nu_mat, function(x){matrix(sapply(x, as.numeric), ncol = n_items, nrow = n_gps, byrow = TRUE)})
  n_group <- blavInspect(bcfa_res, "nobs")
  lambda_aligned <- lapply(1:3000, matrix, data= NA, nrow = n_gps, ncol = n_items)
  nu_aligned <- lapply(1:3000, matrix, data= NA, nrow = n_gps, ncol = n_items)
  # apply alignment on each draw
  for(i in 1:dim(blava_draw_sam)[1]){
    aligned_pars <- invariance.alignment(
      lambda = lambda_mat_new[[i]],
      nu = nu_mat_new[[i]],
      wgt = matrix(sqrt(n_group), ncol = n_items, nrow = n_gps)
    )
    lambda_aligned[[i]] <- aligned_pars$lambda.aligned
    nu_aligned[[i]] <- aligned_pars$nu.aligned
    # print(i)
  }
  list("lambda" = lambda_aligned, n_gps = 2, "nu" = nu_aligned)
}
