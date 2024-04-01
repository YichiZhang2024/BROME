#' Bayesian Region of Measurement Equivalence (ROME) with alignment
#'
#' \code{abrome} implements the Bayesian ROME framework with alignment (ABROME)
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

abrome <- function(pninv, n_items, n_ninv, data, ROME_lim, eta_lim, eta_num){
  ## step1 compute ROME
  pvar <- pooled_sd_test(items = data[,1:n_items], group = data$group)
  rome <- ROME_lim * pvar
  ## step 2 bayesian MG-CFA and alignment
  ba_res <- run_ba(data = data, n_items, n_ninv)
  ly_group <- do.call(rbind, ba_res$lambda)
  ly_groupg1 <- ly_group[rownames(ly_group) == "G1",]
  ly_groupg2 <- ly_group[rownames(ly_group) == "G2",]
  nu_group <- do.call(rbind, ba_res$nu)
  nu_groupg1 <- nu_group[rownames(nu_group) == "G1",]
  nu_groupg2 <- nu_group[rownames(nu_group) == "G2",]
  ## step 3 compute expected total scores
  sum_ld <- list(rowSums(ly_groupg1), rowSums(ly_groupg2))
  sum_nu <- list(rowSums(nu_groupg1), rowSums(nu_groupg2))
  etas <- seq(-eta_lim, eta_lim, length.out = eta_num)
  # ld_list <- list(ly_groupg1, ly_groupg2)
  # nu_list <- list(nu_groupg1, nu_groupg2)
  df_post_diff <- map_df(etas, exp_tot, load = sum_ld, tau = NULL, nu = sum_nu, ordinal = FALSE)
  decision <- rep(NA, 1)
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
