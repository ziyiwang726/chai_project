# This FDP is without cutoff
computeFDP<-function(lFDR,trueIndex,alpha) {
  df<-data.frame(cbind(1:length(lFDR),lFDR))
  df<-df[order(df[,2]),]
  df$cumFDR<-cumsum(df[,2])
  df$cumFDR<-df$cumFDR/(1:nrow(df))
  df<-df[order(df[,1]),]
  select<-which(df$cumFDR<=alpha)
  if(length(select)==0){
    FDP<-0
  }else{
    FDP<-1-length(intersect(select,trueIndex))/length(select)
  }
  FDP
}


computePower<-function(lFDR,trueIndex,alpha) {
  df<-data.frame(cbind(1:length(lFDR),lFDR))
  df<-df[order(df[,2]),]
  df$cumFDR<-cumsum(df[,2])
  df$cumFDR<-df$cumFDR/(1:nrow(df))
  df<-df[order(df[,1]),]
  select<-which(df$cumFDR<=alpha)
  length(intersect(select,trueIndex))/length(trueIndex)
}

computeFDPIndex<-function(index,trueIndex) {
  if(length(index)==0){
    FDP<-0
  }else{
    FDP<-1-length(intersect(index,trueIndex))/length(index)
  }
  FDP
}


computePowerIndex<-function(index,trueIndex) {
  length(intersect(index,trueIndex))/length(trueIndex)
}


ratio_z_given_x_custom<- function(z, post_w, m, v) {
  log_pz <- log(
    sum(post_w * dnorm(z, mean=m, sd=sqrt(v), log=FALSE))
  )
  log_phi_z <- dnorm(z, mean=0, sd=1, log=TRUE)
  # log ratio
  log_ratio <- exp(log_pz - log_phi_z)
  return(log_ratio)  # still can be a big or small #, but less likely Inf
}

minRatioForX_custom <- function(post_w, m, v, lower=-10, upper=10) {
  f <- function(z) ratio_z_given_x_custom(z, post_w, m, v)
  res <- optimize(f, interval=c(lower, upper), maximum=FALSE)
  list(z_star = res$minimum, ratio_star = res$objective)
}

lFDRselect <- function(obj_or_lfdr, q = 0.05, max_lFDR = 1) {
  if (is.list(obj_or_lfdr) && !is.null(obj_or_lfdr$avgFDR)) {
    ord <- obj_or_lfdr$ord
    lFDR_sorted <- if (!is.null(obj_or_lfdr$lFDR_sorted)) {
      obj_or_lfdr$lFDR_sorted
    } else if (!is.null(obj_or_lfdr$clFDR_sorted)) {
      obj_or_lfdr$clFDR_sorted
    } else {
      stop("CHAI object does not contain sorted lFDR/clFDR values.")
    }
    avgFDR <- obj_or_lfdr$avgFDR
  } else {
    lfdr <- as.numeric(obj_or_lfdr)
    ord <- order(lfdr)
    lFDR_sorted <- lfdr[ord]
    avgFDR <- cumsum(lFDR_sorted) / seq_along(lFDR_sorted)
  }

  valid <- which(is.finite(avgFDR) & is.finite(lFDR_sorted) &
                   avgFDR <= q & lFDR_sorted <= max_lFDR)
  if (!length(valid)) return(integer(0))
  ord[seq_len(max(valid))]
}
