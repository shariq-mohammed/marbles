#' Bayesian Information Criterion
#'
#' Computes the Bayesian information criterion from the sequential estimation
#' for Bayesian variable selection in multivariate multiple regression setting.
#'
#' @param Y response matrices; list of length \eqn{\tau} with each list component
#'          as a \eqn{n}x\eqn{p} matrix (\eqn{p} could vary for each list
#'          component)
#' @param X covariate matrix of dimension \eqn{n}x\eqn{g}
#' @param g_id group IDs; list of length \eqn{\tau} with each list component as
#'             a vector of size \eqn{p} with entries \eqn{1,\ldots,m} (\eqn{p}
#'             could vary for each list component)
#' @param b estimated coefficients; list of length \eqn{\tau} with each list
#'          component as a \eqn{g}x\eqn{p} matrix (\eqn{p} could vary for each
#'          list component)
#' @param Delta_inv inverse covariance matrices of the response; response
#'                  matrices; list of length \eqn{\tau} with each list component
#'                  as a \eqn{p}x\eqn{p} matrix (\eqn{p} could vary for each list
#'                  component)
#' @param zeta indicator variables identifying associations between response
#'             groups and predictors; response matrices; list of length
#'             \eqn{\tau} with each list component as a \eqn{g}x\eqn{m} matrix
#' @param p_thres threshold on the indicator; defaults to 0.5
#'
#' @return the computed Bayesian information criterion value
#' @export
#'
#' @examples
#' n = 100
#' g = 10
#' m = 2
#' tau = 3
#' p = 5 # p_1=...=p_tau=p
#'
#' Y = lapply(1:tau, function(t) matrix(rnorm(n*p), nrow = n))
#' X = matrix(rnorm(n*g), nrow = n)
#' g_id = lapply(1:tau, function(t) c(1,1,1,2,2))
#' b = lapply(1:tau, function(t) matrix(rnorm(g*p), nrow = g))
#' Delta_inv = lapply(1:tau, function(t) diag(p))
#' zeta = lapply(1:tau, function(t) matrix(rbinom(g*m, size=1, prob=.5), nrow=g))
#'
#' bic(Y, X, g_id, b, Delta_inv, zeta)

bic = function(Y, X, g_id, b, Delta_inv, zeta, p_thres = 0.5){
  n = nrow(X)
  a = numeric()

  for(t in 1:length(Y)){
    n0b = b[[t]]*(zeta[[t]][,g_id[[t]]]>p_thres)
    t1a = crossprod(Y[[t]] - X%*%n0b)
    t1 = -sum(diag(crossprod(Delta_inv[[t]],t1a)))/2
    t2 = (n/2)*log(det(Delta_inv[[t]]))
    t3 = sum(n0b!=0)
    a[t] = t1+t2-(t3*log(n)/2)
  }
  -2*sum(a)
}
