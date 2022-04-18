#' Sequential Estimation And Model Selection
#'
#' Run sequential estimation for Bayesian variable selection in multivariate
#' multiple regression setting with grouping structure for a one or more
#' values of v0. Computes the Bayesian information criterion (BIC) from the
#' estimation for each value of v0. Return the parameter estimates for the
#' model which has the minimum BIC across different choices of v0.
#'
#' @param Y response matrices; list of length \eqn{\tau} with each list
#'          component as a \eqn{n}x\eqn{p} matrix (\eqn{p} could vary for
#'          each list component)
#' @param X covariate matrix of dimension \eqn{n}x\eqn{g}
#' @param g_id group IDs; list of length \eqn{\tau} with each list component as
#'             a vector of size \eqn{p} with entries \eqn{1,\ldots,m} (\eqn{p}
#'             could vary for each list component)
#' @param temporal should the estimation borrow information sequentially ('yes'
#'                 or 'no'); Defaults to 'yes'
#' @param alpha_mix proportion of information to borrow; a scalar in
#'                  \eqn{[0,1]}; Defaults to 0.5
#' @param p_thres threshold on the indicator; defaults to 0.5
#' @param v0seq sequence of values to consider for v0 to perform model
#'              estimation; Defaults to 10 equally spaced numbers between
#'              0.001 and 0.01
#' @param nCores Integer identifying the number of nodes to be forked (for
#'               details see [parallel::makeCluster]). These nodes will be
#'               used to run estimation in parallel for different choices of
#'               v0; Defaults to 1
#' @param ... pass additional arguments to \code{EM_MBVS}
#'
#' @return \itemize{
##'  \item{\code{v0min}} {The value of v0 which has the lowest corresponding BIC}
##'  \item{\code{bic_res}} {A vector of values for BICs corresponding to the
##'                         choice of values of v0}
##'  \item{\code{seq_res}} {A list object for model parameter estimates with the
##'                         following entries:
##'                         \item{\code{b}} {estimates of the coefficient
##'                                          \eqn{beta}; list of length
##'                                          \eqn{\tau} with each list component
##'                                          as a \eqn{g}x\eqn{p} matrix
##'                                          (\eqn{p} could vary for each list
##'                                          component)}
##'                         \item{\code{w}} {estimates of probability of
##'                                          selecting associations for the
##'                                          covariates with each group of
##'                                          responses; list of length \eqn{\tau}
##'                                          with each list component as a
##'                                          \eqn{g}x\eqn{m} matrix}
##'                         \item{\code{lambda}} {estimates of the structure
##'                                               parameter \eqn{\lambda}; list
##'                                               of length \eqn{\tau} with each
##'                                               list component as a
##'                                               \eqn{g}x\eqn{m} matrix}
##'                         \item{\code{nusq}} {estimates of the prior variance
##'                                             parameter of \eqn{\beta}s; list
##'                                             of length \eqn{\tau} with each
##'                                             list component as a
##'                                             \eqn{g}x\eqn{p} matrix (\eqn{p}
##'                                             could vary for each list
##'                                             component)}
##'                         \item{\code{Delta_inv}} {estimate of the inverse of
##'                                                  residual covariance; list
##'                                                  of length \eqn{\tau} with
##'                                                  each list component as a
##'                                                  \eqn{p}x\eqn{p} matrix
##'                                                  (\eqn{p} could vary for
##'                                                  each list component)}
##'                         \item{\code{bic}} {Bayesian information criterion
##'                                            computed across all sequential
##'                                            regression model estimations; a
##'                                            scalar}
##'                         \item{\code{iter}} {number of iterations completed
##'                                             by the algorithm; a list with
##'                                             each component as an integer
##'                                             representing the number of EM
##'                                             iterations for the corresponding
##'                                             model in the sequence}
##'                         \item{\code{q_fn}} {values attained by the q-function
##'                                             through the iterations; a list
##'                                             with each component as a vector
##'                                             of values attained by the
##'                                             q-function through the iterations
##'                                             for the corresponding model in
##'                                             the sequence}
##'                         \item{\code{v1}} {scalar factor to multiply to the
##'                                           variance in the prior of
##'                                           \eqn{\beta}s when corresponding
##'                                           \eqn{\zeta = 1}}
##'                         }
##' }
#' @export
#' @import doParallel foreach iterators parallel
#'
#' @examples
#' n = 100
#' g = 10
#' m = 2
#' tau = 3
#' p = 5 # p_1=...=p_tau=p
#'
#' # First three columns of Y are first group and next two columns are
#' # second group. Coefficients corresponding to the first three columns
#' # and columns 4, 5 are zero and one, respectively. The same structure
#' # is maintained for t=1,...,tau
#' g_id = lapply(1:tau, function(t) c(1,1,1,2,2))
#' beta = lapply(1:tau, function(t) matrix(rep(c(0,1), g*c(3,2)), nrow = g))
#' X = matrix(rnorm(n*g), nrow = n)
#' Y = lapply(1:tau, function(t) matrix(rnorm(n*p, mean = X%*%beta[[t]]),
#'                                      nrow = n))
#'
#' res = seq_est_model_selection(Y, X, g_id)

seq_est_model_selection = function(Y, X, g_id, temporal = 'yes',
                                   alpha_mix = 0.5, p_thres = 0.5,
                                   v0seq = seq(0.001,0.01,length=10),
                                   nCores = 1, ...){
  stopifnot(nCores>0, nCores%%1==0)

  seqEMestimation = function(Y, X, g_id, temporal = 'yes', alpha_mix = 0.5,
                             p_thres = 0.5, ...){
    # dimensions of the response and predictor
    n = nrow(X)
    tau = length(Y)
    g = ncol(X)

    # source the estimation and BIC computation files.
    seq_res = list(b = list(), w = list(), lambda = list(), nusq = list(),
                   Delta_inv = list(), bic = numeric(), iter = list(),
                   q_fn = list(), v1 = numeric())
    # start sequential estimation
    for(t in 1:tau){
      p = ncol(Y[[t]])
      m = max(g_id[[t]])
      r_name = paste('Region ', t, sep='')
      print(r_name)
      # set the prior mean for lambda
      mu = matrix(0, nrow = g, ncol = m)
      if(t>1 & temporal == 'yes') mu = pmax(mu,
                                            alpha_mix*seq_res$lambda[[t-1]])

      res_temp = EM_MBVS(Y[[t]], X, g_id[[t]], mu = mu, ...)
      seq_res$b[[r_name]] = res_temp$b*(res_temp$w[,g_id[[t]]]>p_thres)
      if(!is.null(res_temp$lambda)) seq_res$lambda[[r_name]] = res_temp$lambda
      if(!is.null(res_temp$nusq)) seq_res$nusq[[r_name]] = res_temp$nusq
      if(!is.null(res_temp$Delta_inv)){
        seq_res$Delta_inv[[r_name]] = res_temp$Delta_inv
      }
      if(!is.null(res_temp$w)) seq_res$w[[r_name]] = res_temp$w
      if(!is.null(res_temp$iter)) seq_res$iter[[r_name]] = res_temp$iter
      if(!is.null(res_temp$q_fn)) seq_res$q_fn[[r_name]] = res_temp$q_fn
      if(!is.null(res_temp$v1)) seq_res$v1[r_name] = res_temp$v1
      rm(res_temp)
    }

    # Compute BIC
    seq_res$bic = bic(Y, X, g_id, seq_res$b, seq_res$Delta_inv,
                      seq_res$w, p_thres)
    seq_res
  }

  if(nCores==1){
    seq_res_v0 = list()
    for(k in 1:length(v0seq)){
      # run the model sequential estimation for a specific v0
      # function defined above
      seq_res_v0[[k]] = seqEMestimation(Y, X, g_id, temporal = 'yes',
                                        alpha_mix = 0.5, p_thres = 0.5,
                                        v0 = v0seq[k], ...)
    }
  } else{
    # register a cluster to estimate the model in parallel
    cl = makeCluster(nCores)
    registerDoParallel(cl)
    seq_res_v0 = foreach(k = 1:length(v0seq)) %dopar% {
      # run the model sequential estimation for a specific v0
      # function defined above
      seqEMestimation(Y, X, g_id, temporal = 'yes', alpha_mix = 0.5,
                      p_thres = 0.5, v0 = v0seq[k], ...)
    }
    stopCluster(cl) # close the cluster
  }

  # check the model BICs
  bic_res = numeric()
  for(a in 1:length(v0seq)) bic_res[a] = seq_res_v0[[a]]$bic
  v0_min = which.min(bic_res) # identify the lowest BIC

  list(v0_min = v0seq[v0_min],
       bic_res = bic_res,
       seq_res = seq_res_v0[[v0_min]])
}
