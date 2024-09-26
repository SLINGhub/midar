# BATCH TREND CORRECTION BASED ON GAUSSIAN KERNEL SMOOTHING
# =========================================================
#
# Original file name: batch-trend-code.R
#
# Version 16. August 2023
# by Hyung Won CHOI, National University of Singapore
#
# ---------------------------------------------------------

### Function to clean up serial trend within a batch
gauss.kernel.smooth = function(xx, yy,
                               kernel.size,
                               outlier.filter=FALSE,
                               outlier.Ksd,
                               location.smooth=TRUE,
                               scale.smooth=FALSE) {

  n = length(yy) ## number of data points

  ## If outlier filter is turned on, mark an outlier as NA
  ## Ksd = K times standard deviation of data distribution
  yy.train = yy - median(yy, na.rm=TRUE)

  if(outlier.filter) {
    mean.y = mean(yy.train, na.rm=TRUE)
    sd.y = sd(yy.train, na.rm=TRUE)
    oid = (abs(yy.train - mean.y) / sd.y) > outlier.Ksd
    oid[is.na(oid)] = FALSE
    yy.train[oid] = NA
  }

  yy.est = yy
  ## Location parameter smoothing
  if(location.smooth) {
    for(i in 1:n) {
      wt = (xx - xx[i]) / kernel.size
      wt = dnorm(wt, 0, 1)
      wt[is.na(yy.train)] = NA
      yy.est[i] = yy[i] - sum(wt * yy.train, na.rm=TRUE) / sum(wt, na.rm=TRUE)
    }
  }

  ## If scale parameter smoothing is on, mean-center data,
  ## estimate point-wise weighted stdev's and scale them
  ## and add back the overall mean
  yy.final = yy.est

  if(scale.smooth) {
    vv = rep(NA, n) ## point-wise weighted variances
    yy.mean = mean(yy.est, na.rm=TRUE)
    yy.est = yy.est - yy.mean
    for(i in 1:n) {
      if(!is.na(yy.est[i])) {
        wt = (xx - xx[i]) / kernel.size
        wt = dnorm(wt, 0, 1)
        wt[is.na(yy.est)] = NA
        vv[i] = sum(wt * yy.est^2, na.rm=TRUE) / sum(wt, na.rm=TRUE)
      }
    }
    vv.mean = mean(vv, na.rm=TRUE) ## average weighted variances across the data points
    ss = sqrt(vv)
    ss.mean = mean(ss, na.rm=TRUE)
    yy.final = yy.mean + yy.est * ss.mean / ss
  }

  ## Report the final values
  yy.final
}


########### Master function to run in-batch smoothing
########### and equalize the mean and variance between batches

batch.correction = function(tab,
                            kernel.size,
                            outlier.filter,
                            outlier.Ksd,
                            location.smooth,
                            scale.smooth,
                            location.crossBatch,
                            scale.crossBatch) {

  batch = tab$Batch
  batch.order = tab$Order
  val = tab$Value

  ubatch = unique(batch)
  nbatch = length(ubatch)

  val.clean = val ## placeholder

  ### Within batch location and scale normalization
  for(b in 1:nbatch) {
    wid = which(batch == ubatch[b])
    ord = as.numeric(batch.order[wid])
    tmp = as.numeric(val[wid])
    if(!all(is.na(tmp))) {
      x = ord
      y = tmp
      val.clean[wid] = gauss.kernel.smooth(x, y,
                                           kernel.size=kernel.size,
                                           outlier.filter=outlier.filter,
                                           outlier.Ksd=outlier.Ksd,
                                           location.smooth=location.smooth,
                                           scale.smooth=scale.smooth)
    }
  }

  ### Cross-batch scale normalization
  if(location.crossBatch) {
    tmp = val.clean
    loc.batch = rep(NA, nbatch)
    for(b in 1:nbatch) {
      id = which(batch == ubatch[b])
      loc.batch[b] = median(tmp[id], na.rm=TRUE)
    }
    loc.batch.mean = mean(loc.batch)
    for(b in 1:nbatch) {
      id = which(batch == ubatch[b])
      xloc = loc.batch[b]
      val.clean[id] = (tmp[id] - xloc) + loc.batch.mean
    }
  }

  if(scale.crossBatch) {
    tmp = val.clean
    loc.batch = rep(NA, nbatch)
    sca.batch = rep(NA, nbatch)
    for(b in 1:nbatch) {
      id = which(batch == ubatch[b])
      loc.batch[b] = median(tmp[id], na.rm=TRUE)
      sca.batch[b] = mad(tmp[id], na.rm=TRUE)
    }
    loc.batch.mean = mean(loc.batch)
    sca.batch.mean = mean(sca.batch)
    for(b in 1:nbatch) {
      id = which(batch == ubatch[b])
      xloc = loc.batch[b]
      val.clean[id] = (tmp[id] - xloc) / sca.batch[b] * sca.batch.mean + loc.batch.mean
    }
  }

  val.clean

}
#
# ###########################
# ### Actual example
# ###########################
# sdata = read.delim("sample_info.txt", header=T, as.is=T, check.names=F)
# qdata = read.delim("sample_data.txt", header=T, as.is=T, check.names=F)
# pdata = read.delim("prot_info.txt", header=T, as.is=T, check.names=F)
#
# sdata$Order = NA
# pp = unique(sdata$PlateNum)
# npp = length(pp)
# for(k in 1:npp) {
#   wid = which(sdata$PlateNum == pp[k])
#   ct = length(wid)
#   sdata$Order[wid] = 1:ct
# }
#
# ######################
# ####### CHEK2 example
# ######################
# wid = which(pdata$Assay == "CHEK2")
# ydata = as.numeric(qdata[wid, ])
#
# tab = data.frame(sdata[,3:4], NPX=ydata, stringsAsFactors = FALSE)
# data.corr = batch.correction(tab,
#                              kernel.size=10,
#                              outlier.filter=TRUE,
#                              outlier.Ksd=5,
#                              location.smooth=TRUE,
#                              scale.smooth=FALSE,
#                              location.crossBatch=TRUE,
#                              scale.crossBatch=TRUE)
#
# par(mfrow=c(2,1))
# borders = c(0,cumsum(table(tab$PlateNum))) + 0.5
# nborders = length(borders)
# ### Before correction
# rg = range(ydata)
# plot(ydata, cex=.2, ylab="Reported Data", main="CHEK2", ylim=rg)
# for(k in 1:nborders) abline(v=borders[k], lty=2, col=2)
# ### After correction
# plot(data.corr, cex=.2, ylab="Normalized Data", main="CHEK2", ylim=rg)
# for(k in 1:nborders) abline(v=borders[k], lty=2, col=2)
#
# ######################
# ####### HK2 example
# ######################
# wid = which(pdata$Assay == "HK2")
# ydata = as.numeric(qdata[wid, ])
#
# tab = data.frame(sdata[,3:4], NPX=ydata, stringsAsFactors = FALSE)
# data.corr = batch.correction(tab,
#                              kernel.size=10,
#                              outlier.filter=TRUE,
#                              outlier.Ksd=5,
#                              location.smooth=TRUE,
#                              scale.smooth=TRUE,
#                              location.crossBatch=TRUE,
#                              scale.crossBatch=TRUE)
#
# par(mfrow=c(2,1))
# borders = c(0,cumsum(table(tab$PlateNum))) + 0.5
# nborders = length(borders)
# ### Before correction
# rg = range(ydata)
# plot(ydata, cex=.2, ylab="Reported Data", main="HK2", ylim=rg)
# for(k in 1:nborders) abline(v=borders[k], lty=2, col=2)
# ### After correction
# plot(data.corr, cex=.2, ylab="Normalized Data", main="HK2", ylim=rg)
# for(k in 1:nborders) abline(v=borders[k], lty=2, col=2)
#
#
