#' Circular Disc Plots
#'
#' Circular disc plots of spherical shells indicating the selection/association
#' of the spherical shells with the predictors.
#'
#' @param w estimates of probability of selecting associations for the covariates
#'          with each group of responses; list of length \eqn{\tau} with each
#'          list component as a \eqn{g}x\eqn{m} matrix. This is usually the
#'          output from \code{seq_est_model_selection} (accessed as
#'          \code{seq_est_model_selection$seq_res$w}). List names should be
#'          region IDs.
#' @param im.seq imaging sequence/group names; should be the same size as the number
#'               of columns in \code{w$`Region 1`}, that is, the number of
#'               groups in the response matrix
#' @param X.colnames column names of the predictor matrix (e.g. gene names)
#' @param p_thres threshold on probability of selection to determine inclusion
#'                of an association; defaults to 0.5
#' @param nr number of rows in [ggplot2::facet_wrap]
#' @param layer TRUE or FALSE indicating if the region label is 'Layer' or
#'              'Region'; Defaults to FALSE such that label is 'Region'
#'
#' @return a \code{ggplot2} object of circular disc plot
#' @export
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#'
#' @examples
#' tau = 3
#' m = 4
#' g = 20
#' w = lapply(1:tau, function(t) matrix(runif(m*g), nrow=g))
#' names(w) = paste0('Region ', 1:tau)
#' im.seq = paste0('IM.', 1:m)
#' X.colnames = LETTERS[1:g]
#' print(disc_plot(w, im.seq, X.colnames))
disc_plot = function(w, im.seq, X.colnames, p_thres = 0.5, nr = 2, layer = F){
  mew = melt(w)
  colnames(mew) = c('Gene','Image Sequence','w','Region')
  mew$`Image Sequence` = im.seq[mew$`Image Sequence`]
  mew$w = as.numeric(mew$w > p_thres)
  mew$Gene = as.character(X.colnames[mew$Gene])
  mew$r_num = sapply(1:nrow(mew),
                     function(i) as.numeric(strsplit(mew$Region[i],
                                                     split = ' ')[[1]][2]))
  tau = max(mew$r_num)
  radius = sapply(1:tau, function(t) (sqrt(t)-sqrt(t-1))/sqrt(tau))
  mew$height = radius[mew$r_num]

  if(layer){
    mew$Region = sapply(1:nrow(mew),
                        function(i) paste0('Layer ',
                                           strsplit(mew$Region[i],
                                                    split = ' ')[[1]][2]))
  }
  mew$Region = as.factor(mew$Region)
  mew$aux = rep(cumsum(sort(unique(mew$height),decreasing = T)),
                table(mew$Region))
  mew$Decision = NA
  mew$Decision[mew$w==0] = 'Not Selected'
  mew$Decision[mew$w==1] = 'Selected'

  defaultW = getOption("warn")
  options(warn = -1)
  bp = ggplot(mew, aes(x=Gene, y=aux))+
    geom_tile(aes(fill = Decision, height = height), colour = "gray50")+
    facet_wrap(~`Image Sequence`, nrow = nr)+
    theme_bw()+
    scale_y_discrete(breaks = cumsum(sort(unique(mew$height), decreasing=T)),
                     limits = cumsum(sort(unique(mew$height), decreasing=T)),
                     labels = levels(mew$Region), expand = c(0,0))+
    scale_fill_manual(values = c("white", "cornflowerblue"))+
    theme(axis.title = element_blank(), axis.text.x = element_text(size=5.49),
          axis.text.y = element_text(size=7), legend.position = "bottom",
          legend.title = element_blank(), panel.grid = element_blank(),
          axis.ticks = element_blank())+
    coord_polar(theta = 'x')
  options(warn = defaultW)

  bp
}
