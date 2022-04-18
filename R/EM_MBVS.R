#' EM Estimation For Bayesian Variable Selection
#'
#' Expectation-Maximization based estimation for Bayesian variable selection
#' in multivariate multiple regression setting with grouping structure. Bayesian
#' variable selection can incorporate correlation structures through the priors.
#'
#' @param Y_t response matrix of dimension \eqn{n}x\eqn{p}
#' @param X covariate matrix of dimension \eqn{n}x\eqn{g}
#' @param g_id_t group IDs vector of size \eqn{p} with entries \eqn{1,\ldots,m}
#' @param b_init initial values for \eqn{\beta}; matrix of dimension
#'               \eqn{g}x\eqn{p}. Defaults to (a) OLS estimate OR (b) OLS
#'               estimate with g-inverse (if \eqn{X} is not full rank) OR
#'               (c) random values sampled from normal distribution. The
#'               default choice is made in the priority order (a),(b),(c)
#' @param Delta_init initial value for residual covariance \eqn{\Delta};
#'                   matrix of dimension \eqn{p}x\eqn{p}. Defaults to identity
#'                   matrix
#' @param nusq_init initial value for prior variance parameter of \eqn{\beta}s;
#'                  matrix of dimension \eqn{g}x\eqn{p}. Defaults to 1
#' @param lambda_init initial value for structure parameter; matrix of dimension
#'                    \eqn{g}x\eqn{m}. Defaults to 0
#' @param mu prior mean of \eqn{\lambda}; matrix of dimension \eqn{g}x\eqn{m}.
#'           Defaults to 0
#' @param Lambda prior covariance of the structure parameter \eqn{\lambda};
#'               matrix of dimension \eqn{g}x\eqn{g}. Defaults to covariance
#'               of X
#' @param v0 scalar factor to multiply to the variance in the prior of
#'           \eqn{\beta}s when corresponding \eqn{\zeta = 0}. Defaults to 0.005
#' @param v1 scalar factor to multiply to the variance in the prior of
#'           \eqn{\beta}s when corresponding \eqn{\zeta = 1}. Defaults to
#'           smallest power of 10 greater than the maximum magnitude of
#'           \code{b_init}
#' @param a1 shape parameter for Gamma prior on \eqn{nusq}. Defaults to 4
#' @param a2 rate parameter for Gamma prior on \eqn{nusq}. Defaults to 5
#' @param delta degrees of freedom for the inverse-Wishart prior on
#'              \eqn{\Delta}. Defaults to number of columns of \eqn{X}.
#' @param Psi prior covariance for the inverse-Wishart prior on \eqn{\Delta};
#'            matrix of dimension \eqn{p}x\eqn{p}. Defaults to covariance
#'            of Y_t
#' @param maxIter maximum number of iterations for the EM algorithm; Defaults
#'                to 1000
#' @param epsilon error threshold; defaults to 1e-5
#'
#' @return \itemize{
##'  \item{\code{b}} {estimates of the coefficient \eqn{beta}; matrix of
##'                   dimension \eqn{g}x\eqn{p}}
##'  \item{\code{Delta_inv}} {estimate of the inverse of residual covariance;
##'                           matrix of dimension \eqn{p}x\eqn{p}}
##'  \item{\code{nusq}} {estimates of the prior variance parameter of
##'                      \eqn{\beta}s; matrix of dimension \eqn{g}x\eqn{p}}
##'  \item{\code{lambda}} {estimates of the structure parameter \eqn{\lambda};
##'                        matrix of dimension \eqn{g}x\eqn{m}}
##'  \item{\code{w}} {estimates of probability of selecting associations for
##'                   the covariates with each group of responses; matrix of
##'                   dimension \eqn{g}x\eqn{m}}
##'  \item{\code{iter}} {number of iterations completed by the algorithm}
##'  \item{\code{q_fn}} {values attained by the q-function through the
##'                      iterations}
##'  \item{\code{v1}} {scalar factor to multiply to the variance in the prior
##'                    of \eqn{\beta}s when corresponding \eqn{\zeta = 1}}
##' }
#' @export
#'
#' @importFrom MASS ginv
#' @importFrom stats cov rnorm
#'
#' @examples
#' n = 100
#' g = 10
#' p = 5
#'
#' # First three columns of Y are first group and next two columns are
#' # second group. Coefficients corresponding to the first three columns
#' # and columns 4,5 are zero and one, respectively.
#' g_id_t = c(1,1,1,2,2)
#' beta = matrix(rep(c(0,1), g*c(3,2)), nrow = g)
#' X = matrix(rnorm(n*g), nrow = n)
#' Y_t = matrix(rnorm(n*p, mean = X%*%beta), nrow = n)
#'
#' res = EM_MBVS(Y_t, X, g_id_t)

EM_MBVS = function(Y_t, X, g_id_t, b_init = NULL, Delta_init = NULL,
                   nusq_init = NULL, lambda_init = NULL, mu = NULL,
                   Lambda = NULL, v0 = 0.005, v1 = NULL, a1 = 4, a2 = 5,
                   delta = NULL, Psi = NULL, maxIter = 1000, epsilon = 1e-5){
  # matrix dimensions and group indices
  n = nrow(Y_t)
  p = ncol(Y_t)
  g = ncol(X)
  stopifnot(is.numeric(g_id_t), min(g_id_t)>=1, g_id_t%%1==0)
  m = max(g_id_t)

  # initialize parameter values
  if(is.null(b_init)){
    if(g<=n) b_init = chol2inv(chol(crossprod(X)))%*%t(X)%*%Y_t
    if(g>n){
      b_init = tryCatch({
        ginv(crossprod(X))%*%t(X)%*%Y_t
      }, error=function(cond){
        matrix(rnorm(g*p, sd=v1^2), nrow=g, ncol=p)
      })
    }
  }
  if(is.null(Delta_init)) Delta_init = diag(p)
  if(is.null(lambda_init)) lambda_init = matrix(0, nrow=g, ncol=m)
  if(is.null(nusq_init)) nusq_init = matrix(1, nrow=g, ncol=p)

  # Set hyperparameters
  if(is.null(mu)) mu = matrix(0, nrow=g, ncol=m)
  if(is.null(Lambda)){
    Lambda = cov(X)
    print('Lambda is considered as covariance of X.')
  }
  if(is.null(v1)) v1 = 10^ceiling(log10(max(abs(b_init))))
  if(is.null(delta)) delta = p
  if(is.null(Psi)){
    Psi = diag(p)
    print('Psi is considered as identity matrix.')
  }

  # Compute inverse covaraince of the structure parameter
  if(g<=n) Lambda_inv = chol2inv(chol(Lambda))
  if(g>n) Lambda_inv = ginv(Lambda)

  b = b_init
  Delta_inv = chol2inv(chol(Delta_init))
  nusq = nusq_init
  lambda = lambda_init

  Z = kronecker(X, diag(p))

  inv_temp = 0.01 # temperature parameter
  objfn = numeric()
  diff_norm = numeric()
  iter = 0
  while(iter<maxIter){
    iter = iter+1

    # E-Step
    pip = inclusion_prob(lambda, b, g_id_t, nusq, v1, v0)
    a_prob = pip$a_prob
    b_prob = pip$b_prob

    w = (a_prob^inv_temp)/((a_prob^inv_temp)+(b_prob^inv_temp))
    d = ((1-w)/v0)+(w/v1)

    # M-Step
    ZT_invZ = kronecker(crossprod(X), Delta_inv)
    GM = numeric()
    for(l in 1:m) GM = cbind(GM, nusq[,g_id_t==l]/d[,l])
    GMinv = diag(1/c(t(GM)))
    ZT_invy = kronecker(t(X), Delta_inv)%*%c(t(Y_t))

    bVar.chol = chol(ZT_invZ + GMinv)
    bVar = chol2inv(bVar.chol)
    b_new = matrix(crossprod(bVar, ZT_invy), ncol=p, byrow = T)

    nusq_new = (a2+((b_new^2)*d[,g_id_t]/2))/(a1-0.5)

    Delta_inv_new = (delta+n+p+1)*chol2inv(chol(Psi+crossprod(Y_t-(X%*%b_new))))

    lambda_new = sapply(1:m,
                        function(l){
                          gradient_descent(lambda[,l], mu[,l], Lambda,
                                           Lambda_inv, w[,l])$lambda
                        })

    diff_norm[iter] = sum((b_new-b)^2)+sum((Delta_inv_new-Delta_inv)^2)+
      sum((nusq_new-nusq)^2)+sum((lambda_new-lambda)^2)

    objfn[iter] = Q_function(Y_t, X, g_id_t, b, Delta_inv, lambda, nusq,
                             a1, a2, mu, Lambda_inv, delta, Psi, w, d)

    b = b_new
    Delta_inv = Delta_inv_new
    nusq = nusq_new
    lambda = lambda_new

    inv_temp = inv_temp*1.1
    if(inv_temp>1) inv_temp = 1
    if(iter %% 50 == 0){
      print(paste('Iter = ', iter,'; Obj Fn = ', round(objfn[iter],7), sep= ''))
      print(paste('Temperature', round(inv_temp,4),sep=" = "))
    }
    if(diff_norm[iter] < epsilon){
      print(paste('Estimates converged! Total Iterations = ',iter,'.',sep=''))
      print(paste('Final Objective Function Value = ', round(objfn[iter],7),
                  sep=''))
      break
    }
  }

  pip = inclusion_prob(lambda, b, g_id_t, nusq, v1, v0)
  a_prob = pip$a_prob
  b_prob = pip$b_prob
  w = a_prob/(a_prob+b_prob)

  list(b = b,
       Delta_inv = Delta_inv,
       nusq = nusq,
       lambda = lambda,
       w = w,
       iter = iter,
       q_fn = objfn,
       v1 = v1)
}
