#' Principal Component Analysis on Densities From Images
#'
#' Performs principal component analysis (PCA) on the space of probability
#' density functions (PDFs) generated from the tumor voxels (or other
#' regions of interest) from the image data. Step 1: Divide the tumor region
#' from each imaging sequence into \eqn{\tau} spherical shells of equal volume.
#' Step 2: Estimate the PDFs for each region and imaging sequence. Step 3: For
#' each region and each imaging sequence, compute PCA to obtain principal
#' component scores.
#'
#' @param image_data Image data (MRI scans from various sequences) and the
#'                   segmentation mask. A list with each entry corresponding
#'                   to the imaging sequences and mask (Note: One of the
#'                   entries has to be named 'mask'). Each list entry is
#'                   again a list (of same length across all imaging sequences)
#'                   with its entries corresponding to the subjects. The
#'                   i-th list entries in each imaging sequence and the mask
#'                   correspond to subject i and are 3D-array of the same
#'                   dimensions. For example see \code{mri_data} in the package
#' @param n.regions number of spherical shells to divide the tumor into; defaults
#'                  to 3
#' @param den.nx number of points to evaluate the PDF at; defaults to 1000
#' @param perc.var percentage variance explained by the principal components
#'                 (PCs) to determine the number of PCs to include; defaults to
#'                 99\%
#' @param nCores Integer identifying the number of nodes to be forked (for
#'               details see [parallel::makeCluster]). These nodes will be
#'               used to run the analysis in parallel for different imaging
#'               sequences; Defaults to 1
#'
#' @return \itemize{
##' \item{\code{pdfs}} {A list of size \eqn{\tau = n.regions}. Each list
##'                     entry is a list (size = number of imaging sequences)
##'                     which contains a matrix of size \eqn{n} x \eqn{den.nx}
##'                     of densities corresponding to an imaging sequence,
##'                     where \eqn{n} is the number of subjects.}
##' \item{\code{pc_scores}} {A list of size \eqn{\tau = n.regions}. Each list
##'                          entry is a matrix (size \eqn{n}x\eqn{p}) of PC
##'                          scores for each region. \eqn{n} is the number of
##'                          subjects and \eqn{p} is the number of PCs included
##'                          across all imaging sequences in \code{image_data}.}
##' \item{\code{g_id}} {A list of size \eqn{\tau = n.regions}. Each list entry
##'                     is a vector (size \eqn{p}) of group indicators for each
##'                     region, where \eqn{p} is the number of PCs included
##'                     across all imaging sequences in \code{image_data}.}
##' \item{\code{avg_pdf}} {A list of size of the number of imaging sequences in
##'                        \code{image_data}. Each list entry is a matrix (size
##'                        \eqn{n.regions} x \eqn{den.nx}) of average PDFs.}
##' }
#' @export
#' @import doParallel
#'
#' @examples
#' n = 3
#' # The following will take a couple of minutes for n=3;
#' # Approximately 30 mins for n=63 in mri_data.
#'
#' # Extracting data from first 3 subjects from MRI data
#' image_data = lapply(names(mri_data), function(im) mri_data[[im]][1:n])
#' names(image_data) = names(mri_data)
#'
#' res = pca_on_images(image_data, n.regions = 3, nCores = 4)

pca_on_images = function(image_data, n.regions = 3, den.nx = 1000,
                         perc.var = 99, nCores = 1){

  perc.var = perc.var/100
  stopifnot(perc.var > 0)
  if(perc.var>1) perc.var = 1

  # identify the imaging sequence names in the data
  im.seq = names(image_data)
  stopifnot('mask' %in% im.seq)
  im.seq = im.seq[im.seq!='mask']
  # number of subjects
  n = length(image_data$mask)
  # min and max values for the imaging sequences
  rngs = sapply(im.seq,
                function(im){
                  range(sapply(1:length(image_data[[im]]),
                               function(i) range(image_data[[im]][[i]])))
                })
  rownames(rngs) = c('min','max')

  if(nCores==1) print(cat('Increase the number of cores (nCores) to speed',
                          'up the computation by running PCA in parallel',
                          'for different imaging sequences.',
                          sep='\n'))

  cl = makeCluster(nCores)
  registerDoParallel(cl, cores = nCores)
  dmp = foreach(im = im.seq) %dopar% {
    # initialize variable to save densities for each region across all subjects
    den = lapply(1:n.regions, function(t) matrix(nrow = n, ncol = den.nx))
    names(den) = paste0('Region.', 1:n.regions)
    # identify x values where the densities are computed at
    x = t(sapply(im.seq,
                 function(im) seq(rngs['min',im], rngs['max',im], length = den.nx)))
    # compute the region-wise density for each subject
    for(i in 1:n){
      # identify the tumor and mask arrays
      msk = image_data[['mask']][[i]]
      tum = image_data[[im]][[i]]*msk

      # identify tumor voxel indices, their minimums, maximums and centers
      ind = apply(which(msk!=0, arr.ind = T), 2, range)
      min_r = ind[1,1]; min_c = ind[1,2]; min_s = ind[1,3]
      max_r = ind[2,1]; max_c = ind[2,2]; max_s = ind[2,3]
      cent = round(colMeans(ind))
      cent_r = cent[1]; cent_c = cent[2]; cent_s = cent[3]
      nr = length(min_r:max_r)
      nc = length(min_c:max_c)
      ns = length(min_s:max_s)

      # create three arrays with entries as the row, column and slice indices
      # To match meshgrid function in MATLAB
      rowsInIm = colsInIm = sliceInIm = array(dim=c(nr,nc,ns))
      for(k in 1:nr) rowsInIm[k,,] = k
      for(k in 1:nc) colsInIm[,k,] = k
      for(k in 1:ns) sliceInIm[,,k] = k

      # tumor volume/size
      tot_non_zero_Vox = sum(msk!=0)
      # number of voxels in each region
      vox_Count = round(tot_non_zero_Vox/n.regions)

      # Initialization
      rad = 0 # radius covered
      tot_Vox = 0 # total voxels included
      # list with the voxel identifiers for each region
      circle_Vox = list()
      circle_Vox[[1]] = array(0, dim=c(nr,nc,ns))

      # identify voxels for each spherical region
      s = 1
      while((s<=n.regions) & (tail(tot_Vox,n=1) != tot_non_zero_Vox)){
        s = s+1
        rad[s] = rad[s-1]
        non_zero_Vox = 0

        while(non_zero_Vox < vox_Count){
          rad[s] = rad[s]+0.03
          circle_Vox[[s]] = 1*(((rowsInIm-cent_r)^2) +
                                 ((colsInIm-cent_c)^2) +
                                 ((sliceInIm-cent_s)^2) <= (rad[s]^2))
          tot_Vox[s] = sum((msk*circle_Vox[[s]])>0)
          non_zero_Vox = tot_Vox[s] - tot_Vox[s-1]

          if(tot_Vox[s] == tot_non_zero_Vox) break
        }
      }

      # generate densities for each spherical region
      for(s in 2:length(rad)){
        den_Vox = msk*(circle_Vox[[s]] - circle_Vox[[s-1]])
        nz_vox = tum[den_Vox>0] # non-zero voxel values

        fn1 = density(nz_vox, from = rngs['min',im],
                      to = rngs['max',im], n=den.nx)$y

        # scaling estimated density with the imag-seq max and min
        fn1 = fn1/(rngs['max',im] - rngs['min',im])
        fn1 = fn1/sum(diff(seq(0, 1, length=den.nx))*(fn1[-length(fn1)]+fn1[-1])/2)
        den[[paste0('Region.',s-1)]][i,] = fn1
      }
    }

    # PCA with densities
    # square-root transform of the densities
    sqrt_den = lapply(1:n.regions, function(t) sqrt(den[[paste0('Region.', t)]]))
    names(sqrt_den) = paste0('Region.', 1:n.regions)

    # compute Karcher mean on the sqrt space
    kmean_sqrt = lapply(1:n.regions,
                        function(t){
                          sqgd = sqrt_den[[paste0('Region.', t)]]
                          eps = .5
                          nmv = 100
                          sqgdbar = rep(1, den.nx)

                          iter = 1
                          while(nmv[iter]>1e-10){
                            vbar = rep(0, den.nx)
                            v = matrix(nrow = n, ncol = den.nx)
                            for(i in 1:n){
                              if(abs(sum(sqgdbar-sqgd[i,]))<1e-10){
                                v[i,] = rep(0,m)
                              } else{
                                th = acos(sum(diff(seq(0,1,length=den.nx))*
                                                ((sqgdbar*sqgd[i,])[-den.nx]+
                                                   (sqgdbar*sqgd[i,])[-1])/2))
                                v[i,] = (th/sin(th))*(sqgd[i,]-sqgdbar*cos(th))
                              }
                              vbar = vbar + v[i,]/n
                            }
                            v = eps*vbar

                            nv = sqrt(sum(((v^2)[-length(v^2)]+
                                             (v^2)[-1])/2)/(den.nx-1))
                            sqgdbar = cos(nv)*sqgdbar + sin(nv)*(eps*vbar)/nv

                            iter = iter+1
                            nmv_t = vbar*vbar
                            nmv = append(nmv, sqrt(sum((nmv_t[-length(nmv_t)]+
                                                          nmv_t[-1])/2)/(den.nx-1)))

                          }
                          sqgdbar
                        })
    names(kmean_sqrt) = paste0('Region.', 1:n.regions)

    # compute Karcher mean on the sqrt space
    kmean = lapply(1:n.regions, function(t) kmean_sqrt[[paste0('Region.', t)]]^2)
    names(kmean) = paste0('Region.', 1:n.regions)

    # project sqrt transforms onto the tangent space of Karcher mean
    # compute the sample covariance matrix
    # Compute PC Scores
    pcScores = lapply(1:n.regions,
                      function(t){
                        p = kmean_sqrt[[paste0('Region.', t)]]
                        ieden = matrix(0, nrow = n, ncol = den.nx)
                        for(i in 1:n){
                          q = sqrt_den[[paste0('Region.', t)]][i,]
                          if(abs(sum(p-q))<1e-10){
                            ieden[i,] = rep(0,den.nx)
                          } else{
                            th = acos(sum(diff(seq(0,1,length=den.nx))*
                                            ((p*q)[-den.nx]+ (p*q)[-1])/2))
                            ieden[i,] = (th/sin(th))*(q-p*cos(th))
                          }
                        }

                        K = matrix(0, nrow = den.nx, ncol = den.nx)
                        for(i in 1:n) K = K+tcrossprod(ieden[i,])
                        K = K/(n-1)

                        svdecomp = svd(K)
                        p.n = sum(cumsum((svdecomp$d^2)/sum(svdecomp$d^2))<perc.var)+1
                        ieden%*%(svdecomp$u[,1:p.n])
                      })
    names(pcScores) = paste0('Region.', 1:n.regions)

    list(den = den,
         kmean = kmean,
         pcScores = pcScores)
  }
  stopCluster(cl)
  names(dmp) = im.seq

  # densities
  den = lapply(1:n.regions,
                     function(t){
                       X = list()
                       for(im in im.seq){
                         X[[im]] = dmp[[im]]$den[[paste0('Region.', t)]]
                       }
                       X
                     })
  names(den) = paste0('Region.', 1:n.regions)

  # pc scores
  pc_scores = lapply(1:n.regions,
                     function(t){
                       X = numeric()
                       for(im in im.seq) {
                         X_t = dmp[[im]]$pcScores[[paste0('Region.', t)]]
                         colnames(X_t) = paste(im, 1:ncol(X_t), sep ='.')
                         X = cbind(X, X_t)
                       }
                       rownames(X) = names(image_data[[im]])
                       X
                     })
  names(pc_scores) = paste0('Region.', 1:n.regions)

  # group ids for pc scores
  g_id = lapply(1:n.regions,
                function(t){
                  g = numeric()
                  count=0
                  for(im in im.seq) {
                    count = count+1
                    X_t = dmp[[im]]$pcScores[[paste0('Region.', t)]]
                    g_t = rep(count, ncol(X_t))
                    names(g_t) = paste(im, 1:ncol(X_t), sep ='.')
                    g = append(g, g_t)
                  }
                  g
                })
  names(g_id) = paste0('Region.', 1:n.regions)

  # pdf averages
  avg_pdf = lapply(im.seq,
                   function(im){
                     kmean = numeric()
                     for(t in 1:n.regions){
                       kmean = rbind(kmean,
                                     dmp[[im]]$kmean[[paste0('Region.', t)]])
                     }
                     rownames(kmean) = paste0('Region.', 1:n.regions)
                     kmean
                   })
  names(avg_pdf) = im.seq

  list(pdfs = den, pc_scores = pc_scores, g_id = g_id, avg_pdf = avg_pdf)
}
