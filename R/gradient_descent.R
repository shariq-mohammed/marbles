#' Gradient Descent Algorithm For Structure Parameter
#'
#' Performs the M-step optimization using gradient descent approach in the
#' EM algorithm for the structure parameter \eqn{\lambda}. Note: The
#' maximization problem was transformed to solve a minimization problem,
#' however, the objective function is returned for the maximization.
#'
#' @param lambda_init initial value of \eqn{\lambda}; vector of size \eqn{g}
#' @param mu prior mean of \eqn{\lambda}; vector of size \eqn{g}
#' @param L prior covariance of \eqn{\lambda}; matrix of dimension
#'          \eqn{g}x\eqn{g}
#' @param L_inv inverse prior covariance of \eqn{\lambda}; matrix of dimension
#'              \eqn{g}x\eqn{g}
#' @param w probability of selecting associations for the covariates; vector of
#'          size \eqn{g}
#' @param alpha step-size for gradient descent; defaults to 0.01 but dynamically
#'              changes through the iterations
#' @param epsilon error threshold; defaults to 1e-10
#' @param maxIter maximum number of iterations; defaults to 1000
#'
#' @return \itemize{
##'  \item{\code{lambda}} {Solution of the gradient descent for the structure
##'                        parameter \eqn{\lambda}}
##'  \item{\code{objFn}} {Values of the objective function (from the gradient
##'                       descent for the structure parameter \eqn{\lambda})
##'                       attained through iterations}
##' }
#' @export
#'
#' @importFrom stats dnorm pnorm
#'
#' @examples
#' g = 10
#'
#' lambda_init = rnorm(g)
#' mu = rep(c(0,1), c(g/2,g/2))
#' L = diag(g)
#' L_inv = solve(L)
#' w = sort(runif(g))
#'
#' res = gradient_descent(lambda_init, mu, L, L_inv, w)

gradient_descent = function(lambda_init, mu, L, L_inv, w, alpha = 0.01,
                            epsilon = 1e-10, maxIter = 1000){
  lambda = lambda_init
  alpha_init = alpha
  g = length(mu)

  diff_norm = epsilon+1
  Iter = 1

  t1 = crossprod(c(lambda-mu), L_inv)%*%c(lambda-mu)/2
  t2 = -sum(((1-w)*log(pnorm(-lambda)))+(w*log(pnorm(lambda))))
  objFn = -t1-t2
  count = 0
  while(count<maxIter){
    count = count+1
    del_pdf = dnorm(lambda)
    del_cdf = pnorm(lambda)
    del2 = del_pdf*(w-del_cdf)/(del_cdf*(1-del_cdf))
    del1 = c(crossprod(c(lambda-mu), L_inv))
    del = alpha*(del1-del2)
    lambda_new = lambda - del

    t1 = crossprod(c(lambda_new-mu), L_inv)%*%c(lambda_new-mu)/2
    t2 = -sum(((1-w)*log(pnorm(-lambda_new)))+(w*log(pnorm(lambda_new))))
    oF = -t1-t2
    if(oF < objFn[Iter]){
      alpha = alpha*0.9
      if(alpha==0) break
      next
    }
    if(oF >= objFn[Iter]){
      Iter = Iter+1
      lambda = lambda_new
      objFn[Iter] = oF
      alpha = alpha*1.1
      if(objFn[Iter]-objFn[Iter-1] < epsilon) break
    }
    diff_norm = sum(del^2)
    if(diff_norm < epsilon) break
  }
  list(lambda = lambda,
       objFn = objFn)
}
