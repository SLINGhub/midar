#' Function for Gaussian kernel-based smoothing, for use by `corr_drift_fun`
#' @param tbl Table (`tibble` or `data.frame`) containing the fields `qc_type`, `x` (run order number), and `y` (variable)
#' @param qc_types QC types used for the smoothing (fit) by loess
#' @param span_width Bandwidth of the gaussian kernel function
#' @param ... Additional parameters forwarded to KernSmooth::locpoly
#' @return List with a `data.frame` containing original `x` and the smoothed `y` values, and a `boolean` value indicting whether the fit failed or not not.

fun_gaussiankernel <- function(tbl, qc_types, span_width, ...) {
  arguments <- list(...)
  # browser()
  d_subset <- tbl[tbl$qc_type %in% qc_types, ] |> tidyr::drop_na(.data$y)
  res <- tryCatch(
    {
      fit <- KernSmooth::locpoly(d_subset$x, d_subset$y, bandwidth = span_width, gridsize = nrow(tbl), range.x = c(min(tbl$x), max(tbl$x)))
      fit$y
    },
    error = function(e) {
      print(e$message)
      return(rep(NA_real_, length(tbl$x)))
    }
  )
  list(res = res, fit_error = all(is.na(res)))
}


#
#
# ### Function to clean up serial trend within a batch
# gauss.kernel.smooth = function(tbl,
#                                qc_types,
#                                kernel_size,
#                                outlier_filter,
#                                outlier_ksd = 5,
#                                location_smooth = TRUE,
#                                scale_smooth = FALSE) {
#
#   arguments <- list(...)
#   #browser()
#   d_subset <- tbl[tbl$qc_type %in% qc_types, ] |> tidyr::drop_na(.data$y)
#
#   n = length(d_subset$y) ## number of data points
#
#   ## If outlier filter is turned on, mark an outlier as NA
#   ## Ksd = K times standard deviation of data distribution
#   d_subset$y_train = d_subset$y - median(d_subset$y, na.rm=TRUE)
#
#   if(outlier_filter) {
#     d_subset$y_mean = mean(d_subset$y_train, na.rm=TRUE)
#     d_subset$y_sd = sd(d_subset$y_train, na.rm=TRUE)
#     d_subset$oid = (abs(d_subset$y_train -  d_subset$y_mean) / d_subset$y_sd) > outlier_ksd
#     d_subset$oid[is.na(d_subset$oid)] = FALSE
#     d_subset$y_train.train[d_subset$oid] = NA
#   }
#
#   d_subset$y_est = d_subset$y
#   ## Location parameter smoothing
#   if(location_smooth) {
#     for(i in 1:n) {
#       d_subset$wt = (d_subset$x - d_subset$x[i]) / kernel_size
#       d_subset$wt = dnorm(d_subset$wt, 0, 1)
#       d_subset$wt[is.na(d_subset$yy)] = NA
#       d_subset$y_est[i] = d_subset$y[i] - sum(d_subset$wt * d_subset$y_train, na.rm=TRUE) / sum(d_subset$wt, na.rm=TRUE)
#     }
#   }
#
#   ## If scale parameter smoothing is on, mean-center data,
#   ## estimate point-wise weighted stdev's and scale them
#   ## and add back the overall mean
#   d_subset$y_final =  d_subset$y_est
#
#   if(scale_smooth) {
#     d_subset$v = rep(NA, n) ## point-wise weighted variances
#     d_subset$y_mean = mean(d_subset$y_est, na.rm=TRUE)
#     d_subset$y_est = d_subset$y_est - d_subset$y_mean
#     for(i in 1:n) {
#       if(!is.na(d_subset$y_est[i])) {
#         d_subset$wt = (d_subset$x - d_subset$x[i]) / kernel_size
#         d_subset$wt = dnorm(d_subset$wt, 0, 1)
#         d_subset$wt[is.na(d_subset$y_est)] = NA
#         d_subset$v[i] = sum(d_subset$wt * d_subset$y_est^2, na.rm=TRUE) / sum(d_subset$wt, na.rm=TRUE)
#       }
#     }
#     d_subset$v_mean = mean(d_subset$v, na.rm=TRUE) ## average weighted variances across the data points
#     d_subset$ss = sqrt(d_subset$v)
#     d_subset$ss_mean = mean(d_subset$ss, na.rm=TRUE)
#     d_subset$y_final = d_subset$y_mean + yy.est * ss.mean / ss
#   }
#
#   ## Report the final values
#   d_subset
# }
#
#
# ########### Master function to run in-batch smoothing
# ########### and equalize the mean and variance between batches
#
# batch.correction = function(tab,
#                             kernel.size,
#                             outlier.filter,
#                             outlier.Ksd = 5,
#                             location.smooth,
#                             scale.smooth,
#                             location.crossBatch,
#                             scale.crossBatch) {
#
#   batch = tab$PlateNum
#   batch.order = tab$Order
#   npx = tab$NPX
#
#   ubatch = unique(batch)
#   nbatch = length(ubatch)
#
#   npx.clean = npx ## placeholder
#
#   ### Within batch location and scale normalization
#   for(b in 1:nbatch) {
#     wid = which(batch == ubatch[b])
#     ord = as.numeric(batch.order[wid])
#     tmp = as.numeric(npx[wid])
#     if(!all(is.na(tmp))) {
#       x = ord
#       y = tmp
#       npx.clean[wid] = gauss.kernel.smooth(x, y,
#                                            kernel.size=kernel.size,
#                                            outlier.filter=outlier.filter,
#                                            outlier.Ksd=outlier.Ksd,
#                                            location.smooth=location.smooth,
#                                            scale.smooth=scale.smooth)
#     }
#   }
#
#   ### Cross-batch scale normalization
#   if(location.crossBatch) {
#     tmp = npx.clean
#     loc.batch = rep(NA, nbatch)
#     for(b in 1:nbatch) {
#       id = which(batch == ubatch[b])
#       loc.batch[b] = median(tmp[id], na.rm=TRUE)
#     }
#     loc.batch.mean = mean(loc.batch)
#     for(b in 1:nbatch) {
#       id = which(batch == ubatch[b])
#       xloc = loc.batch[b]
#       npx.clean[id] = (tmp[id] - xloc) + loc.batch.mean
#     }
#   }
#
#   if(scale.crossBatch) {
#     tmp = npx.clean
#     loc.batch = rep(NA, nbatch)
#     sca.batch = rep(NA, nbatch)
#     for(b in 1:nbatch) {
#       id = which(batch == ubatch[b])
#       loc.batch[b] = median(tmp[id], na.rm=TRUE)
#       sca.batch[b] = mad(tmp[id], na.rm=TRUE)
#     }
#     loc.batch.mean = mean(loc.batch)
#     sca.batch.mean = mean(sca.batch)
#     for(b in 1:nbatch) {
#       id = which(batch == ubatch[b])
#       xloc = loc.batch[b]
#       npx.clean[id] = (tmp[id] - xloc) / sca.batch[b] * sca.batch.mean + loc.batch.mean
#     }
#   }
#
#   npx.clean
#
# }
