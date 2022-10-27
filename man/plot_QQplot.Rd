\name{plot_QQplot}
\alias{plot_QQplot}
\title{Diagnostic qqplot for GMM decoposition..}
\usage{
plot_QQplot(
      data,
      GModel
)
}
\arguments{
\item{data}{Vector of data to decompose by GMM.}

\item{GModel}{List of GModel parameters i.e GModel$model$alpha, GModel$model$mu, GModel$model$sigma.Given as output of function gaussian_mixture_vector}

}

\description{
Function return Quantile-Quantile plot for one component normal distribution and fited GMM.
}

\seealso{
gaussian_mixture_vector, plot_gmm_1D
}
