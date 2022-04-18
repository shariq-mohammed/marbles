#' Inclusion Probabilities
#'
#' Computes the numerators of the inclusion probability, that is, the posterior
#' probability that the indicator parameter \eqn{\zeta} is 0 or 1.
#'
#' @param lambda structure parameter matrix of dimension \eqn{g}x\eqn{m}
#' @param b coefficients matrix of dimension \eqn{g}x\eqn{p}
#' @param g_id_t group IDs vector of size \eqn{p} with entries \eqn{1,\ldots,m}
#' @param nusq prior variance parameter of \eqn{\beta}; matrix of dimension
#'             \eqn{g}x\eqn{p}
#' @param v1 scalar factor to multiply to the variance in the prior of
#'           \eqn{\beta}s when corresponding \eqn{\zeta = 1}
#' @param v0 scalar factor to multiply to the variance in the prior of
#'           \eqn{\beta}s when corresponding \eqn{\zeta = 0}
#'
#' @return \itemize{
##'  \item{\code{a_prob}} {numerator of the posterior probability of the
##'                        indicator parameter corresponding \eqn{\zeta=1}}
##'  \item{\code{b_prob}} {numerator of the posterior probability of the
##'                        indicator parameter corresponding \eqn{\zeta=0}}
##' }
#' @export
#'
#' @importFrom stats pnorm dnorm
#'
#' @examples
#' g = 10
#' m = 2
#' p = 5
#' v0 = 0.005
#' v1 = 1
#'
#' # First three columns of Y are first group and next two columns are
#' # second group. Coefficients corresponding to the first three columns
#' # and columns 4,5 are zero and one, respectively.
#' g_id_t = c(1,1,1,2,2)
#' b = matrix(rep(c(0,1), g*c(3,2)), nrow = g)
#' nusq = matrix(rgamma(g*p, shape = 1), nrow = g)
#' lambda = matrix(sort(rnorm(g*m)), nrow = g)
#'
#' res = inclusion_prob(lambda, b, g_id_t, nusq, v1, v0)
#' res$a_prob/(res$a_prob+res$b_prob)

inclusion_prob = function(lambda, b, g_id_t, nusq, v1, v0){
  m = max(g_id_t)
  a_prob = sapply(1:m,
                  function(l){
                    term1 = pnorm(lambda[,l])
                    b_k_m = b[,g_id_t==l]
                    p_k_m = ncol(b_k_m)
                    nusq_k_m = nusq[,g_id_t==l]
                    temp_term2 = sapply(1:p_k_m,
                                        function(l2){
                                          dnorm(b_k_m[,l2],
                                                sd=sqrt(v1*nusq_k_m[,l2]))
                                        })
                    term2 = apply(temp_term2, 1, prod)
                    term1*term2
                  })

  b_prob = sapply(1:m,
                  function(l){
                    term1 = 1-pnorm(lambda[,l])
                    b_k_m = b[,g_id_t==l]
                    p_k_m = ncol(b_k_m)
                    nusq_k_m = nusq[,g_id_t==l]
                    temp_term2 = sapply(1:p_k_m,
                                        function(l2){
                                          dnorm(b_k_m[,l2],
                                                sd=sqrt(v0*nusq_k_m[,l2]))
                                        })
                    term2 = apply(temp_term2, 1, prod)
                    term1*term2
                  })
  list(a_prob = a_prob,
       b_prob = b_prob)
}
