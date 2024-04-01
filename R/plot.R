#' Plot the expected total test score for each group
#'
#' \code{plot_expt} plots the expected total test score for each group
#' @param df_post The expected total test scores for the group of interest
#'
#' @return A plot that shows the expected total test scores for the group of interest across latent ability level
#'

plot_expt <- function(df_post){
  ggplot(df_post, aes(x = eta, y = mean)) +
    ylab(labs(y="Expected total test scores"))+
    geom_line(data = df_post, aes(linetype = "Mean"))+
    geom_line(aes(x = eta, y = ll, linetype = "HPDI")) +
    geom_line(aes(x = eta, y = ul, linetype = "HPDI")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_linetype_manual(values = c("Mean"="solid", "HPDI" = "dashed"), labels = c("Posterior Means", "95% HPDI")) +
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = "bottom",legend.key.size = unit(2, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=13),
          axis.title=element_text(size=14,face="bold"))+
    theme_bw()
}


#' Plot the expected difference in total test scores across groups
#'
#' \code{plot_exptDiff} plots the expected difference in total test scores across groups
#'
#' @param diff_post The expected difference in total test scores across groups
#' @param ROME The preset region of measurement equivalence (ROME)
#'
#' @return A plot that shows the expected total test scores for the group of interest across latent ability level
#'

plot_exptDiff <- function(diff_post, ROME){
  ggplot(diff_post, aes(x = eta, y = mean)) +
    geom_ribbon(aes(ymin = ll, ymax = ul), alpha = .15, linetype = 2) +
    geom_line(data = diff_post, color = "orange") +
    geom_hline(yintercept = 0, linetype = "dotdash", color = "red") +
    geom_hline(yintercept = ROME[2], linetype = "dotdash", color = "red") +
    labs(y="Expected Group Difference in Total Test Scores",
         x = "Latent ability") +
    ylim(c(0, 0.5)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_brewer(palette = "Set1") +
    papaja::theme_apa(box = TRUE, base_size = 12) +
    theme(legend.position = "bottom", plot.caption=element_text(hjust = 0))+
    annotate(geom="text", x = 2, y = ROME[1] + 0.04, label="ROME",
             color="red") +
    annotate(geom="text", x = 2, y = ROME[2] - 0.04, label="ROME",
             color="red")

}

