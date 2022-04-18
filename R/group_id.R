#' @title Group Indicators of Principal Component Scores From MRI Data
#'
#' @description The group indicators for the principal component scores computed
#'              from the probability density functions constructed by dividing the
#'              tumor regions in \code{mri_data} into three equal-volume spherical
#'              shells. A list of size 3  corresponding to the three regions. Each
#'              list entry is a vector (size \eqn{p}) of group indicators
#'              (corresponding the list entries of \code{pc_scores}) of PC scores
#'              for each region. \eqn{p} is the number of PCs included across all
#'              imaging sequences in \code{mri_data} for a given spherical shell.
#'              The values for \eqn{p} are 22, 18 and 17 for spherical regions 1,
#'              2 and 3, respectively.
#'
#' @format A list of size 3. Each list entry is a vector. \itemize{
#'         \item{\code{group_id$Region.1}} {a vector of size 22 with group
#'                                          indicators for region 1 (spherical
#'                                          shell 1).}
#'         \item{\code{group_id$Region.2}} {a vector of size 18 with group
#'                                          indicators for region 2 (spherical
#'                                          shell 2).}
#'         \item{\code{group_id$Region.3}} {a vector of size 17 with group
#'                                          indicators for region 3 (spherical
#'                                          shell 3).}
#'         }
"group_id"
