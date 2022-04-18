#' @title MRI Data For LGG Subjects
#'
#' @description Magnetic Resonance Imaging (MRI) scans from FLAIR, T1, T1Gd and
#'              T2 imaging sequences. A list with entries corresponding to the
#'              four imaging sequences and mask. Each list entry is again a list
#'              (of same length across all imaging sequences) with its entries
#'              corresponding to the 63 LGG subjects. The i-th list entries in
#'              each imaging sequence and the mask correspond to subject i and
#'              are 3D-array of the same dimensions.
#'
#' @format A list of size 5. Each list entry is a list of size 63.
"mri_data"
