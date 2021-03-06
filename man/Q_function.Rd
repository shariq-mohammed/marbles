% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Q_function.R
\name{Q_function}
\alias{Q_function}
\title{Q-Function For EM Algorithm}
\usage{
Q_function(
  Y_t,
  X,
  g_id_t,
  b,
  Delta_inv,
  lambda,
  nusq,
  a1,
  a2,
  mu,
  Lambda_inv,
  delta,
  Psi,
  w,
  d
)
}
\arguments{
\item{Y_t}{response matrix of dimension \eqn{n}x\eqn{p}}

\item{X}{covariate matrix of dimension \eqn{n}x\eqn{g}}

\item{g_id_t}{group IDs vector of size \eqn{p} with entries \eqn{1,\ldots,m}}

\item{b}{coefficients matrix of dimension \eqn{g}x\eqn{p}}

\item{Delta_inv}{inverse covariance matrix of the residual; matrix of
dimension \eqn{p}x\eqn{p}}

\item{lambda}{structure parameter matrix of dimension \eqn{g}x\eqn{m}}

\item{nusq}{prior variance parameter of \eqn{\beta}; matrix of dimension
\eqn{g}x\eqn{p}}

\item{a1}{shape parameter for Gamma prior on \eqn{nusq}}

\item{a2}{rate parameter for Gamma prior on \eqn{nusq}}

\item{mu}{prior mean of \eqn{\lambda}; matrix of dimension \eqn{g}x\eqn{m}}

\item{Lambda_inv}{inverse prior covariance of \eqn{\lambda}; matrix of
dimension \eqn{g}x\eqn{g}}

\item{delta}{degrees of freedom for the inverse-Wishart prior on \eqn{\Delta}}

\item{Psi}{prior covariance for the inverse-Wishart prior on \eqn{\Delta};
matrix of dimension \eqn{p}x\eqn{p}}

\item{w}{probability of selecting associations for the covariates with each
group of responses; matrix of dimension \eqn{g}x\eqn{m}}

\item{d}{a function of w (refer to the manuscript); matrix of dimension
\eqn{g}x\eqn{m}}
}
\value{
computed value for the Q-Function for the EM algorithm
}
\description{
Computes the value for the Q-function from the EM algorithm given the values
of the parameters/arguments.
}
\examples{
n = 100
g = 10
m = 2
p = 5

a1 = 4
a2 = 5
v0 = 0.005
v1 = 10

Y_t = matrix(rnorm(n*p), nrow = n)
X = matrix(rnorm(n*g), nrow = n)
g_id_t = c(1,1,1,2,2)
b = matrix(rnorm(g*p), nrow = g)
Delta_inv = diag(p)
lambda = matrix(rnorm(g*m), nrow = g)
nusq = matrix(rgamma(g*p, shape = 1), nrow = g)

mu = matrix(0, nrow = g, ncol = m)
Lambda_inv = diag(g)
delta = p
Psi = diag(p)
w = matrix(runif(g*m), nrow = g)
d = ((1-w)/v0)+(w/v1)

Q_function(Y_t, X, g_id_t, b, Delta_inv, lambda, nusq,
           a1, a2, mu, Lambda_inv, delta, Psi, w, d)
}
