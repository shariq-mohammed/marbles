#' Q-Function For EM Algorithm
#'
#' Computes the value for the Q-function from the EM algorithm given the values
#' of the parameters/arguments.
#'
#' @param Y_t response matrix of dimension \eqn{n}x\eqn{p}
#' @param X covariate matrix of dimension \eqn{n}x\eqn{g}
#' @param g_id_t group IDs vector of size \eqn{p} with entries \eqn{1,\ldots,m}
#' @param b coefficients matrix of dimension \eqn{g}x\eqn{p}
#' @param Delta_inv inverse covariance matrix of the residual; matrix of
#'                  dimension \eqn{p}x\eqn{p}
#' @param lambda structure parameter matrix of dimension \eqn{g}x\eqn{m}
#' @param nusq prior variance parameter of \eqn{\beta}; matrix of dimension
#'             \eqn{g}x\eqn{p}
#' @param a1 shape parameter for Gamma prior on \eqn{nusq}
#' @param a2 rate parameter for Gamma prior on \eqn{nusq}
#' @param mu prior mean of \eqn{\lambda}; matrix of dimension \eqn{g}x\eqn{m}
#' @param Lambda_inv inverse prior covariance of \eqn{\lambda}; matrix of
#'                   dimension \eqn{g}x\eqn{g}
#' @param delta degrees of freedom for the inverse-Wishart prior on \eqn{\Delta}
#' @param Psi prior covariance for the inverse-Wishart prior on \eqn{\Delta};
#'            matrix of dimension \eqn{p}x\eqn{p}
#' @param w probability of selecting associations for the covariates with each
#'          group of responses; matrix of dimension \eqn{g}x\eqn{m}
#' @param d a function of w (refer to the manuscript); matrix of dimension
#'          \eqn{g}x\eqn{m}
#'
#' @return computed value for the Q-Function for the EM algorithm
#' @export
#'
#' @importFrom stats pnorm
#'
#' @examples
#' n = 100
#' g = 10
#' m = 2
#' p = 5
#'
#' a1 = 4
#' a2 = 5
#' v0 = 0.005
#' v1 = 10
#'
#' Y_t = matrix(rnorm(n*p), nrow = n)
#' X = matrix(rnorm(n*g), nrow = n)
#' g_id_t = c(1,1,1,2,2)
#' b = matrix(rnorm(g*p), nrow = g)
#' Delta_inv = diag(p)
#' lambda = matrix(rnorm(g*m), nrow = g)
#' nusq = matrix(rgamma(g*p, shape = 1), nrow = g)
#'
#' mu = matrix(0, nrow = g, ncol = m)
#' Lambda_inv = diag(g)
#' delta = p
#' Psi = diag(p)
#' w = matrix(runif(g*m), nrow = g)
#' d = ((1-w)/v0)+(w/v1)
#'
#' Q_function(Y_t, X, g_id_t, b, Delta_inv, lambda, nusq,
#'            a1, a2, mu, Lambda_inv, delta, Psi, w, d)

Q_function = function(Y_t, X, g_id_t, b, Delta_inv, lambda, nusq,
                      a1, a2, mu, Lambda_inv, delta, Psi, w, d){
  n = nrow(X)
  g = ncol(X)
  p = ncol(Y_t)
  m = ncol(lambda)

  Q1_1 = -sum(diag(crossprod(Y_t - X%*%b)%*%Delta_inv))/2
  Q1_2 = -(a1-0.5)*sum(log(nusq))
  Q1_3 = -sum(((b^2)*d[,g_id_t])/(2*nusq))
  Q1_4 = -a2*sum(1/nusq)
  Q1_5 = (n+delta+p+1)*log(det(Delta_inv))/2
  Q1_6 = -sum(diag(Psi%*%Delta_inv))/2

  Q2_1 = sum(((1-w)*log(1-pnorm(lambda)))+(w*log(pnorm(lambda))))
  Q2_3 = -sum(sapply(1:m,
                     function(l){
                       z = lambda[,l]-mu[,l]
                       crossprod(z,Lambda_inv/2)%*%z
                     }))
  Q1_1+Q1_2+Q1_3+Q1_4+Q1_5+Q1_6+Q2_1+Q2_3
}
