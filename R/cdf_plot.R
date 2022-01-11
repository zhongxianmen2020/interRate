#' @title interRate
#'
#' @description pit (model_cdf) plots
#'
#' @param model_cdf model_cdf
#'
#' @examples
#' cdf_plot(model_cdf)
#' @export



cdf_plot = function(model_cdf){

  x=model_cdf
  z= (x-min(x))/(max(x)-min(x))
 #hist1=hist(z, breaks =15)

  # dev.off()
  par(mar=c(2.5, 2.5, 2.5, 2.5))
  # par(mar = rep(2,4))
  par(mfrow= c(2,1))

  op=par(mfrow = c(2,1),
         #omg=c(2,2,0,0)+0.05,
         mar=c(2.5,2.5,3.0,3.0) +0.05)

  hist1= hist(z, breaks =15, xlab="", ylab="",main="", lwd=2, cex.axis=1.2)
  #plot(hist1,  main="", xlab="", ylab="", lwd=2, cex.axis=1.2)
  mu1=mean(hist1$counts)
  sd1=sd(hist1$counts)
  xx1=seq(0,1, by=0.01)
  y_up= rep(mu1 +sd1*1.96, length(xx1))
  y_down= rep(mu1 -sd1*1.96, length(xx1))


  title(main="CDF comparison", cex.lab=1.2, col.lab="black")


  lines(xx1, y_up, lwd=2, col="blue", main="CDF comparison", cex.axis=2)
  lines(xx1, y_down, lwd=2, col="blue", main="CDF comparison", cex.axis=2)


  a=ecdf(model_cdf)

  plot(a, col="blue", main ="", cex.axis=1.2, xlab="", ylab="", xlim=c(0,1))
  x=seq(from=0, to=1, by= 0.01)

  yy=x
  lines(x, yy, type='l', col="red", lwd=3)

  title(main="CDF comparison", cex.lab=1.2, col.lab="black")

  legend("topleft", legend=c("Theoretical CDF", "Empirical CDF"),
         col=c("red", "blue"), lty=c(1,1),
         lwd=2,  cex=0.9)

  ks.test(model_cdf, "punif",0,1)
}
