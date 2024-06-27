\name{generated_data}
\alias{generated_data}
\docType{data}
\title{
This data set demonstrates the \code{po_fit} function and consists of arbitrarily censored data with left truncation.
}
\description{
This data set was generated using the method described in the paper "Regression analysis of arbitrarily censored and left-truncated data under the proportional odds model."
}
\usage{data("generated_data")}
\format{
  This data set contains 500 observations and 2 covariates.
  \describe{
    \item{\code{L}}{a numeric vector of the lower bounds of the observed intervals}
    \item{\code{R}}{a numeric vector of the upper bounds of the observed intervals}
    \item{\code{C}}{a numeric vector of the left truncation times}
    \item{\code{x_1}}{The first covariate, following a standard normal distribution}
    \item{\code{x_2}}{The second covariate, following a Bernoulli(0.5) distribution}
  }
}
\details{
Data is generated following the PO model for failure times, defined as:
\deqn{F(t | \mathbf{x}) = \frac{\Lambda_0(t) \exp(x_1 \beta_1 + x_2 \beta_2)}{1 + \Lambda_0(t) \exp(x_1 \beta_1 + x_2 \beta_2)},}
where \eqn{x_1 \sim N(0,1)}{x_1 ~ N(0,1)} and \eqn{x_2 \sim \text{Bernoulli}(0.5)}{x_2 ~ Bernoulli(0.5)}. The true values of \eqn{\beta_1} and \eqn{\beta_2} are taken from the set \{-1, 0, 1\}, and the true baseline odds function is \eqn{\Lambda_0(t) = \log(1+t) + t^3 + \sin(t)}.

For each subject i, the left truncation time A_i is generated from \eqn{\text{Gamma}(0.5, 1)}{Gamma(0.5, 1)}. The failure time T_i is obtained by numerically solving the equation:
\deqn{F(t_i | a_i, x_i) = P(T_i \leq t_i | T_i \geq a_i, x_i) = u_i,}
where u_i is drawn from a uniform distribution U(0,1). A binary indicator w_i is generated from a Bernoulli distribution with a success probability of 0.2. If w_i = 1, T_i is observed exactly; otherwise, T_i is treated as interval-censored. The observed interval [L_i, R_i] for T_i is determined by generating examination times from cumulative sums of gap times, independently drawn from an exponential distribution with a mean of 0.1, leading to left-truncated and arbitrarily censored data.
}


%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Wang, L., and Wang, L. (2024+)
\emph{Regression analysis of arbitrarily censored and left-truncated data under the proportional odds model}
Submitted for publication.
}

\keyword{datasets}
