#' @title interRate
#'
#' @description QQ plot of the standardized residuals.
#'
#' @param data list of the standardized residuals. 
#'
#' @examples
#' QQ_plot(df)
#' @export


QQ_plot = function(ress){
  op=par(mfrow = c(1,1),
         #omg=c(2,2,0,0)+0.05,
         mar=c(2.5,2.5,3.0,3.0) +0.05)
  
  min1 = min(ress)
  max1 = max(ress)
  
  qqnorm(ress, lwd=2, cex.axis=1.2, main ="Normal QQ plot",cex.lab=1.6, 
         xlim=c(min1, max1), ylim=c(min1, max1))
  #qqline(ress, col=2, lwd=2)
  
  abline(0,1, col="red", lwd=3)
  
  ks.test(ress, "pnorm",0,1)
}