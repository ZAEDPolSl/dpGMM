\name{runGMM}
\alias{runGMM}
\title{Runing of gaussian mixutre modeling (GMM) on data vector.}

\description{
Run of GMM for 1D vector.
}

\usage{
runGMM(
  X,
  KS,
  Y = NULL,
  change = Inf,
  max_iter = 5000,
  SW=NULL,
  IC = "BIC",
  merge = TRUE,
  sigmas.dev = 2.5,
  precision=1e4,
  plot = TRUE,
  col.pal="Blues",
  quick_stop = TRUE,
  signi = 0.05
)
}
\arguments{
\item{X}{Vector of data to decompose by GMM}

\item{KS}{Number of components to check}

\item{Y}{Optional parameter for binned data. Gives value of data counts}

\item{change}{Nie wiem}

\item{max_iter}{Maximum no. of EM iterations}

\item{SW}{Minimum standard devation of component. Default set to range(x)/(5*no.of.components))^2}

\item{IC}{Information Criterion to select best no. of components. Possible "AIC","AICc", "BIC" (default), "ICL-BIC" or "LR"}

\item{merge}{Merging components distance by set standard deviations(sigmas.dev). Default = TRUE}

\item{sigmas.dev}{Maxiumum standard devations of merging}

\item{precision}{Precision of fit plot}

\item{plot}{Plot decomposition with histogram. Default = TRUE}

\item{col.pal}{Palett for ploting components from RColorBrewer}

\item{quick_stop}{Stop the algorithm if LR test is not sigificant for further steps. Default = TRUE}

\item{signi}{Significance level for LR test. Default 0.05}
}



\examples{
data("GMMexample")
runGMM(data,KS=10)
}
