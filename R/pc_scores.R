#' @title Principal Component Scores From MRI Data
#'
#' @description The principal component scores computed from the probability
#'              density functions constructed by dividing the tumor regions in
#'              \code{mri_data} into three equal-volume spherical shells. A list
#'              of size 3 corresponding to the three regions. Each list entry is
#'              a matrix (size 63x\eqn{p}) of PC scores for each region. \eqn{p}
#'              is the number of PCs included across all imaging sequences in
#'              \code{mri_data} for a given spherical shell. The values for
#'              \eqn{p} are 22, 18 and 17 for spherical regions 1, 2 and 3,
#'              respectively. Column names  of the matrices indicate the imaging
#'              sequence name and the PC ID. Row names indicate the subject ID.
#'
#' @format A list of size 3. Each list entry is a matrix with 63 rows. \itemize{
#'         \item{\code{pc_scores$Region.1}} {a 63x22 matrix with PC scores from
#'                                           region 1 (spherical shell 1).}
#'         \item{\code{pc_scores$Region.2}} {a 63x18 matrix with PC scores from
#'                                           region 2 (spherical shell 2).}
#'         \item{\code{pc_scores$Region.3}} {a 63x17 matrix with PC scores from
#'                                           region 3 (spherical shell 3).}
#'         }
"pc_scores"
