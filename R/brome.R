#' Bayesian Region of Measurement Equivalence (ROME) Framework
#'
#' \code{brome} implements the Bayesian ROME framework for Establishing Measurement Equivalence
#'
#' @param data The dataset
#' @param n_items Number of items
#' @param n_ninv Number of noninvariant items
#' @param pninv Proportion of noninvariance
#' @param ROME_lim Preset Region of Measurement Equivalence (ROME)
#' @param eta_lim Range of latent ability
#' @param eta_num Number of latent ability level
#'
#' @return A vector of length 6.
#'          \item{ROME}{preset ROME}
#'          \item{mean}{posterior mean of the expected group difference on total test scores}
#'          \item{sd}{posterior SD of the expected group difference on total test scores}
#'          \item{PDImax}{maximum length of HPDI for the expected group difference on total test scores}
#'          \item{PDImean}{mean length of HPDI for the expected group difference on total test scores}
#'          \item{decision}{decision on whether the scale is practically invariant across the group of interest}
#' @export

brome <- function(n_items, n_ninv, pninv, data, ROME_lim, eta_lim, eta_num){
  ## step1 compute ROME
  pvar <- pooled_sd_test(items = data %>% select(-"group"),  group = data$group)
  rome <- ROME_lim * pvar
  ## step 2 fit the correct partial invariance model
  bcfa_res <- fit_true(data = data, n_items = n_items, n_ninv = n_ninv,
                       pninv = pninv, n_gps = 2, svp = FALSE)
  blava_draws <- blavInspect(bcfa_res, "draws")
  ### Combine results from three chains to one chain
  bpi_draw_sam <- map_df(blava_draws, as.data.frame)
  ## step 3 compute the expected difference in total scale scores
  ly_group <- bpi_draw_sam[ , str_which(names(bpi_draw_sam), "ly")[1:(n_items*2)]]
  ly_groupg1 <- matrix(as.numeric(unlist(ly_group[, 1:n_items])), nrow = 3000, ncol = n_items)
  ly_groupg2 <- matrix(as.numeric(unlist(ly_group[, (n_items+1):(n_items*2)])), nrow = 3000, ncol = n_items)
  sum_ld <- list(rowSums(ly_groupg1), rowSums(ly_groupg2))
  nu_group <- bpi_draw_sam[ , str_which(names(bpi_draw_sam), "Nu")[1:(n_items*2)]]
  nu_groupg1 <- matrix(as.numeric(unlist(nu_group[, 1:n_items])), nrow = 3000, ncol = n_items)
  nu_groupg2 <- matrix(as.numeric(unlist(nu_group[, (n_items+1):(n_items*2)])), nrow = 3000, ncol = n_items)
  sum_nu <- list(rowSums(nu_groupg1), rowSums(nu_groupg2))
  etas <- seq(-eta_lim, eta_lim, length.out = eta_num)
  # ld_list <- list(ly_groupg1, ly_groupg2)
  # nu_list <- list(nu_groupg1, nu_groupg2)
  df_post_diff <- map_df(etas, exp_tot, load = sum_ld, tau = NULL, nu = sum_nu, ordinal = FALSE)
  decision <- vector()
  ## step 4 compare total score and HPDI
  if(all(df_post_diff$ll > -rome) && all(df_post_diff$ul < rome)){
    decision <- TRUE
  }else{
    decision <- FALSE
  }
  c("ROME" = rome, "mean" = mean(df_post_diff$mean), "sd" = mean(df_post_diff$sd),
    "PDImax" = max(abs(df_post_diff$ul - df_post_diff$ll)),
    "PDImean" = mean(abs(df_post_diff$ul - df_post_diff$ll)),
    "decision" = decision)
}
