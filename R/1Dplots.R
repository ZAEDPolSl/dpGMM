#' Plot  of GMM decomposition for 1D data
#'
#' Function plot the decomposed distribution together with histogram of data. Moreover the cut-off are marked.
#' This plot is also return as regular output of \code{\link{runGMM}}.
#'
#' @param X Vector of 1D data for GMM decomposition.
#' @param dist Output of \code{\link{generate_dist}} function.
#' @param Y Vector of counts, with the same length as "X".
#' Applies only to binned data (Y = NULL, by default).
#' @param threshold Vector with GMM cutoffs.
#' @param pal Name of the RColorBrewer palette used in the figure. By default \code{"Blues"}.
#'
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom grDevices colorRampPalette
#' @importFrom reshape2 melt
#'
#' @seealso \code{\link{runGMM}}
#'
#' @examples
#' \dontrun{
#' data(example)
#'
#' alpha <- c(0.45, 0.5, 0.05)
#' mu <- c(-14, -2, 5)
#' sigma <- c(2, 4, 1.5)
#'
#' dist.plot <- generate_dist(example$Dist, alpha, mu, sigma, 1e4)
#' thr <- find_thr_by_params(alpha, mu, sigma, dist.plot)
#' plot_gmm_1D(example$Dist, dist.plot, Y = NULL, threshold = thr, pal="Dark2")
#' }
#'
#' @export
plot_gmm_1D <- function(X, dist, Y = NULL, threshold = NA, pal = "Blues"){

  # estimating of bin width of histogram
  binwidth <- (max(X) - min(X))/floor(sqrt(length(X)))
  #binwidth = 2.64*IQR(data$V1)*nrow(data)^(-1/3) # alternative
  # extract data from dist
  x_temp <- dist$x
  dist <- dist$dist

  # organize data for line plot
  colnames(dist) <- paste('Comp:', 1:ncol(dist), sep='')
  colnames(dist)[ncol(dist)] <- "Main"


  tmp <- reshape2::melt(dist, id.vars = NULL)
  tmp$xx <- rep(x_temp, ncol(dist))
  tmp$lin <- 1
  tmp$lin[which(tmp$variable == "Main")] <- 0

  if (ncol(dist) == 2){
    col <- c("darkgreen", "grey25")
  } else{
    col <- grDevices::colorRampPalette(brewer.pal(8, pal))(ncol(dist))
    col <- c(col[2:ncol(dist)], "grey25")}


  p <- ggplot() + theme_bw()

  if (!is.null(Y)){
    binwidth <- (max(X) - min(X))/length(X)
    Y <- Y/(sum(Y) * binwidth)
    p <- p + geom_bar(aes(x = X, y = Y), stat = "identity", color = "grey65", fill = "grey65", alpha = 0.2)

  } else{
    p <- p + geom_histogram(aes(x = X, y = after_stat(density)), binwidth = binwidth, color = "black", fill = "grey", alpha = 0.2)
  }

  p <- p + geom_line(aes(x = tmp$xx, y = tmp$value, group = tmp$variable, color = as.factor(tmp$variable), linetype = as.factor(tmp$lin)), linewidth = 1) +
    scale_color_manual(values = col, name="") +
    scale_linetype_manual(values = c("dashed","solid")) + xlab("x") + ylab("Density") + guides(linetype="none")


  threshold <- threshold[is.na(threshold) == F]
  if (sum(!is.na(threshold))){
    p <- p + geom_vline(xintercept = threshold, lty = "dashed", col = "red")
  }

  return(p)
}

#' QQplot of GMM decomposition for 1D data
#'
#' Function return ggplot object with fit diagnostic Quantile-Quantile plot for one normal distribution and fitted GMM.
#' This plot is also return as regular output of \code{\link{runGMM}}.
#'
#' @param X Vector of 1D data for GMM decomposition.
#' @param GModel \code{data.frame} of GMM parameters i.e GModel$alpha, GModel$mu, GModel$sigma (correct \code{colnames} are obligatory).
#'
#' @import ggplot2
#' @importFrom  ggpubr ggarrange
#' @importFrom stats qqplot
#'
#' @examples
#' \dontrun{
#' data(example)
#'
#' alpha <- c(0.45, 0.5, 0.05)
#' mu <- c(-14, -2, 5)
#' sigma <- c(2, 4, 1.5)
#'
#' plot_QQplot(example$Dist, alpha, mu, sigma)
#' }
#'
#' @seealso \code{\link{runGMM}}
#'
#' @export
plot_QQplot <- function(X, alpha, mu, sigma){

  GModel <- data.frame(alpha = alpha,
                       mu = mu,
                       sigma = sigma)

  tor <- generate_norm1D(length(X), GModel$alpha, GModel$mu, GModel$sigma)
  quants <- stats::qqplot(X, tor$Dist, plot.it = F)
  tmp <- data.frame(data = quants$x, theor = quants$y)

  p2 <- ggplot(tmp, aes(theor, data)) + theme_bw() +
    geom_smooth(method = 'lm', formula = y ~ x, color = "red", se = F) +
    geom_point(size = 2.5, shape = 18, color = "#324376") +
    ylab("Data") + xlab("GMM fit") +
    ggtitle("QQ plot: GMM") + theme(plot.title = element_text(hjust = 0.5))

  p1 <- ggplot(tmp, aes(sample = data)) + stat_qq(size=2.5, shape=18, color="#324376") + stat_qq_line(color="red")+
    ylab("Data") + xlab("Normal distribution") + theme_bw()+
    ggtitle("QQ plot: one dist.") + theme(plot.title = element_text(hjust = 0.5))

  return(ggpubr::ggarrange(p1, p2, align="hv"))
}
