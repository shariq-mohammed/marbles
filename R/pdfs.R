#' @title Probability Density Functions From MRI Data
#'
#' @description The probability density functions computed by dividing the tumor
#'              regions in \code{mri_data} into three equal-volume spherical
#'              shells. A list of size 3 corresponding to the three regions.
#'              Each list entry is a list (size 4 corresponding to the four MRI
#'              sequences FLAIR, T1, T1Gd and T2. Each sublist entry is a matrix
#'              of size \eqn{63}x\eqn{1000} corresponding to \eqn{63} subjects.
#'
#' @format A list of size 3. Each list entry is a list of size 4. Each sublist
#'         entry is a matrix of size \eqn{63}x\eqn{1000} corresponding to \eqn{63}
#'         subjects.
"pdfs"
