\name{syntheticRcpp-package}
\alias{syntheticRcpp-package}
\alias{syntheticRcpp}
\docType{package}
\title{
syntheticRcpp
}
\author{Hayoung Cheon}

\description{
This package employs synthetic control method to estimate time-varying treatment effect on treated.
}
\details{

The syntheticRcpp is R::synth's implementation in Rcpp. 
The optimal weight W is found with ADMM, and the optimal V is found with BFGS.
ATT trajectory can be drawn with resampling method of inferenece. 
}


\references{
1.	Abadie, A. and Gardeazabal, J. (2003) Economic Costs of Conflict: A Case Study of the Basque Country. American Economic Review 93 (1) 113-132.
2.	Abadie A, Diamond A, Hainmueller J. (2010) Synthetic Control Methods for Comparative Case Studies: Estimating the Effect of California's Tobacco Control Program. Journal of the American Statistical Association 105 (490) 493-505.
}


\examples{
In Github repository, I uploaded sample code as example. For easier comparison, I uploaded basque data which are in R::synth package. 
}
