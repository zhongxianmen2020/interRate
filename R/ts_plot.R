#' @title interRate
#'
#' @description plot the time series draw from the conditional distributions.
#'
#' @param est A datafram of the time series
#'
#' @examples
#' ts_plot(est)
#' @export

ts_plot = function(est){

  m = dim(est)[2]

  if (m==3) {
  op=par(mfrow = c(m,2),
         #omg=c(2,2,0,0)+0.05,
         mar=c(2.5,2.5,3.0,3.0) +0.05)

  nn=1:dim(est)[1]
  plot(nn, est[,1], type="l", lwd=2, col="blue",xlab="", ylab="",
       main=expression(paste("Time series of ", alpha)), cex.main=1.5,
       cex.lab=1.9, cex.axis = 1.5)
  hist(est[,1], lwd=2, col="blue",xlab="", ylab="",
       main=expression(paste("Histogram of ", alpha)), cex.main=1.5,
       cex.lab=1.9, cex.axis = 1.5)


  plot(nn, est[,2], type="l", lwd=2, col="blue",
       main=expression(paste("Time series of ", beta)), cex.main=1.5,
       cex.lab=1.9, cex.axis = 1.5)
  hist(est[,2],  lwd=2, col="blue",xlab="", ylab="",
       main=expression(paste("Histogram of ", beta)), cex.main=1.5,
       cex.lab=1.9, cex.axis = 1.5)


  plot(nn, est[,3], type="l", lwd=2, col="blue",
       main=expression(paste("Time series of ", sigma)), cex.main=1.5,
       cex.lab=1.9, cex.axis = 1.5)

  hist(est[,3],  lwd=2, col="blue",xlab="", ylab="",
       main=expression(paste("Histogram of ", sigma)), cex.main=1.5,
       cex.lab=1.9, cex.axis = 1.5)

}

  if (m == 4) {
    op=par(mfrow = c(m,2),
           #omg=c(2,2,0,0)+0.05,
           mar=c(2.5,2.5,3.0,3.0) +0.05)

    nn=1:dim(est)[1]
    plot(nn, est[,1], type="l", lwd=2, col="blue",xlab="", ylab="",
         main=expression(paste("Time series of ", alpha)), cex.main=1.5,
         cex.lab=1.9, cex.axis = 1.5)
    hist(est[,1], lwd=2, col="blue",xlab="", ylab="",
         main=expression(paste("Histogram of ", alpha)), cex.main=1.5,
         cex.lab=1.9, cex.axis = 1.5)


    plot(nn, est[,2], type="l", lwd=2, col="blue",
         main=expression(paste("Time series of ", beta)), cex.main=1.5,
         cex.lab=1.9, cex.axis = 1.5)
    hist(est[,2],  lwd=2, col="blue",xlab="", ylab="",
         main=expression(paste("Histogram of ", beta)), cex.main=1.5,
         cex.lab=1.9, cex.axis = 1.5)

    plot(nn, est[,3], type="l", lwd=2, col="blue",
         main=expression(paste("Time series of ", sigma)), cex.main=1.5,
         cex.lab=1.9, cex.axis = 1.5)

    hist(est[,3],  lwd=2, col="blue",xlab="", ylab="",
         main=expression(paste("Histogram of ", sigma)), cex.main=1.5,
         cex.lab=1.9, cex.axis = 1.5)

    plot(nn, est[,4], type="l", lwd=2, col="blue",
         main=expression(paste("Time series of ", gamma)), cex.main=1.5,
         cex.lab=1.9, cex.axis = 1.5)

    hist(est[,4], lwd=2, col="blue",xlab="", ylab="",
         main=expression(paste("Histogram of ", gamma)), cex.main=1.5,
         cex.lab=1.9, cex.axis = 1.5)

  }


}
# stop here

