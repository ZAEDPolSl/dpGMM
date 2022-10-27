\name{plot_QQplot}
\alias{plot_QQplot}
\title{QQplot for GMM decomposition}

\description{
Function return ggplot object with fit diagnostic Quantile-Quantile plot for one normal distribution and fited GMM.
The QQplot is also return as regular output of runGMM.
}

\usage{
plot_QQplot(
      data,
      GModel
)
}

\arguments{
\item{data}{Vector of data}

\item{GModel}{List of GModel parameters i.e GModel$model$alpha, GModel$model$mu, GModel$model$sigma. Given as output of function gaussian_mixture_vector}
}


\examples{
data<-generate_norm1D(1000, alpha=c(0.2,0.4,0.4), mu=c(-15,0,15), sigma=c(1,2,3))
GModel<-list(model=list(alpha=c(0.2,0.4,0.4), mu=c(-15,0,15), sigma=c(1,2,3)))
plot_QQplot(data,GModel)
}


\seealso{
gaussian_mixture_vector, plot_gmm_1D, generate_norm1D
}
