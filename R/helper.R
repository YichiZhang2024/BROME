#' Computing the pooled standard deviation of total test scores
#'
#' \code{pooled_sd_test} computes the pooled SD of total test scores
#'
#' @param items Item names for the test.
#' @param group The column for the grouping characteristics
#'
#' @return A single number
#'

pooled_sd_test <- function(items, group){

  ave_gs <- apply(items, 2, function(item) ave(item, group))
  sqrt(sum((rowSums(items) - rowSums(ave_gs))^2)/(nrow(items)-length(unique(group))))
}


#' Computing the expected total test score
#'
#' \code{exp_tot} computes the expected total test score for each group
#'
#' @param eta Latent ability
#' @param si Intercepts for the group of interest
#' @param sl Loadings for the group of interest
#' @param prob the percentage of HPDI
#'
#' @return A vector of length 5.
#'          \item{eta}{latent ability}
#'          \item{mean}{posterior mean total test scores}
#'          \item{sd}{posterior standard deviation of the total test scores}
#'          \item{ll}{lower bound of HPDI}
#'          \item{ul}{upper bound of HPDI}
#'

exp_tot <- function(eta, si, sl, prob = 0.95) {
  tot_sam <- sl * eta + si
  # Compute HPDI:
  hpdi <- coda::HPDinterval(coda::as.mcmc(tot_sam), prob = prob)
  # Return posterior means and HPDI
  tibble(
    eta = eta,
    mean = mean(tot_sam),
    sd = sd(tot_sam),
    ll = hpdi[1],
    ul = hpdi[2]
  )
}


#' Computing the cumulative expected group difference in total test score
#'
#' \code{exp_diff_cum} computes the cumulative expected group difference in
#' total test score for all groups
#'
#' @param eta Latent ability
#' @param si Difference in intercepts
#' @param sl Difference in loadings
#' @param sil The product of the difference in intercepts and difference in loadings
#' @param prob the percentage of HPDI
#'
#' @return A vector of length 5.
#'          \item{eta}{latent ability}
#'          \item{mean}{posterior mean total test scores}
#'          \item{sd}{posterior standard deviation of the total test scores}
#'          \item{ll}{lower bound of HPDI}
#'          \item{ul}{upper bound of HPDI}
#'

exp_diff_cum <- function(eta, si, sl, sil, prob = 0.90) {
  tot_sam <- sqrt(sl * (eta^2) + 2 * eta * sil + si)
  # Compute HPDI:
  hpdi <- coda::HPDinterval(coda::as.mcmc(tot_sam), prob = prob)
  # Return posterior means and HPDI
  tibble(
    eta = eta,
    mean = mean(tot_sam),
    sd = sd(tot_sam),
    ll = hpdi[1],
    ul = hpdi[2]
  )
}
