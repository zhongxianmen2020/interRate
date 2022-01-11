#' @title interRate
#'
#' @description Fits four types of interest rate models.
#'
#' @param data A data set object
#' @param burn_in The sample size of the burn in before estimating the parameters
#' @param N Sample size of the totally sampled time series
#' @param model_type Specify the type of the model that will fitted.(e.g. OU, MOU, CIR, or CEV)
#'
#' @return (1) The estimated parameters and their confidence intervals,
#' (2) the sampled time series after a burn_in period for parameter estimation,
#' (3) the probability integral transform of the fitted model based on the training data,
#' (4) the standardized residuals , and (5) the AIC an BIC.
#'
#' @import nimble
#' @import HDInterval
#' @examples
#' interRate(data, burn_in, N, model_type)
#' @export


interRate = function(data, burn_in, N, model_type){
  y = data
  if (sum(is.na(y)) != 0){
    return("There are some missing values in the data set!")
  }


  if (!((is.numeric(burn_in)) & (burn_in%%1 ==0)& (burn_in>0))){
    return("burn_in has to be a positive integer!")
  }


  if (!((is.numeric(N)) & (N%%1 ==0)& (N>0))){
    return("The sample size of MCMC time sereis has to be a positive integer!")
  }


  if(burn_in >= N){
    return ("Burin_in can not be larger than the simulated sample size")
  }

  if (!(toupper(model_type) %in% c("OU", "MOU", "CIR", "CEV"))){
    return("model_type has to be of of the values c('OU', 'MOU', 'CIR', 'CEV')!")
  }

  if (toupper(model_type) %in% c("OU")){

    NN=length(y)

    alpha1=0.001
    beta1= 0.002
    sigma1=0.001

    alpha_est1=0
    beta_est1 =0
    sigma_est1 =0

    n = burn_in

    m = N

    est = matrix(rep(0, (m-n)*3), nrow=m-n)


    for (i in 1:m){
      #print(i)
      c=0
      d=0
      for (t in 1:(NN-1))
      {
        #c=c+1/y[t]^2
        c=c+1
        d=d+(y[t+1]-(1-beta1)*y[t])
      }

      a=d/c
      b=sigma1/sqrt(c)
      alpha1=rnorm(1, a, b)


      c=0
      d=0
      for (t in 1:(NN-1))
      {
        c=c+y[t]^2
        d=d+(alpha1 +y[t] -y[t+1]  )*y[t]
      }

      a=d/c
      b=sigma1/sqrt(c)
      beta1=rnorm(1, a, b)

      c=0
      for (t in 1:(NN-1))
      {
        c=c+ (y[t+1]-alpha1 -(1-beta1)*y[t])^2
      }

      c=c/2

      sigma1=sqrt(rinvgamma(1,(NN / 2-1), c))

      #d1 <- data.frame(i, alpha1, beta1, sigma1)
      #write.table(d1, "cir-mcmc-R.txt", row.names = FALSE,
      #            col.names = FALSE, append = TRUE)

      #print(c(alpha1, beta1, sigma1))

      if (i >n){

        alpha_est1=alpha_est1 +alpha1
        beta_est1 = beta_est1 +beta1
        sigma_est1 = sigma_est1 +sigma1
        est[i-n,]=c(alpha1, beta1, sigma1)
      }
    }
    alpha_est1=alpha_est1/(m-n)
    beta_est1 = beta_est1/(m-n)
    sigma_est1 = sigma_est1/(m-n)


    # print(paste0("Estimated parameter values monthly y"))
    # print(c(beta_est1, alpha_est1/beta_est1, sigma_est1, mean(y)))
    # hdi(est)
    #
    # print(paste0("Estimated parameter values monthly y"))
    # print(c(alpha_est1, beta_est1, sigma_est1, mean(y)))
    # hdi(est)

    mean(est[,1])
    mean(est[,2])
    mean(est[,3])


    sd(est[,1])
    sd(est[,2])
    sd(est[,3])

    estimate = c(alpha_est1, beta_est1, sigma_est1)
    sd1 = apply(est,2, sd)
    cf = t(hdi(est))
    df= data.frame(estimate, sd1, cf)
    colnames(df) = c("Est.", "Std.", "HPD CI(95%)_lower", "HPD CI(95%)_upper")

    df_simu_time_series = data.frame(est)
    colnames(df_simu_time_series) = c("alpha", "beta", "simga")




    L_y=length(y)
    L1=L_y -1
    ress=rep(0, L1)
    mu=rep(0, L1)
    vol=rep(0, L1)
    rr=rep(0, L1)

    aic=0
    bic=0
    c=0
    for (i in 2:L_y){
      mu[i-1] = y[i-1] + alpha_est1 -beta_est1*y[i-1]
      rr[i-1] = y[i]-mu[i-1]
      vol[i-1] = sigma_est1
      ress[i-1] = rr[i-1]/vol[i-1]
      c = c -2*log(dnorm(rr[i-1], 0, vol[i-1]))
    }

    bic = c + 2 * log(L1)
    aic = c + 2 * 3


    model_cdf=rep(0, L1)
    for (i in 2:L_y){
      model_cdf[i] = pnorm(y[i], y[i-1] +alpha_est1 -beta_est1 *y[i-1], sigma_est1)
    }
    #ks.test(model_cdf, "punif",0,1)

    lt = list("estimate" = df, "standard_residuals"=ress , "mcmc_series" = df_simu_time_series, "Aic" = aic, "Bic"=bic, "empir_cdf"=model_cdf)

    #newlist <- list(lt,mn)

    return(lt)


  }


  # now we fit MOU model

  if (toupper(model_type) %in% c("MOU")){

    NN=length(y)

    alpha1=0.001
    beta1= 0.002
    sigma1=0.001

    alpha_est1=0
    beta_est1 =0
    sigma_est1 =0

    n = burn_in

    m = N

    est = matrix(rep(0, (m-n)*3), nrow=m-n)


    for (i in 1:m){
      #print(i)
      c=0
      d=0
      for (t in 1:(NN-1))
      {
        #c=c+1/y[t]^2
        c=c+1/y[t]^2
        d=d+(y[t+1]-(1-beta1)*y[t])/y[t]^2
      }

      a=d/c
      b=sigma1/sqrt(c)
      alpha1=rnorm(1, a, b)


      c=0
      d=0
      for (t in 1:(NN-1))
      {
        c=c+1
        d=d+(alpha1 +y[t] -y[t+1]  )/y[t]
      }

      a=d/c
      b=sigma1/sqrt(c)
      beta1=rnorm(1, a, b)

      c=0
      for (t in 1:(NN-1))
      {
        c=c+ (y[t+1]-alpha1 -(1-beta1)*y[t])^2 /y[t]^2
      }

      c=c/2

      sigma1=sqrt(rinvgamma(1,(NN / 2-1), c))

      #d1 <- data.frame(i, alpha1, beta1, sigma1)
      #write.table(d1, "cir-mcmc-R.txt", row.names = FALSE,
      #            col.names = FALSE, append = TRUE)

      #print(c(alpha1, beta1, sigma1))

      if (i >n){

        alpha_est1=alpha_est1 +alpha1
        beta_est1 = beta_est1 +beta1
        sigma_est1 = sigma_est1 +sigma1
        est[i-n,]=c(alpha1, beta1, sigma1)
      }
    }
    alpha_est1=alpha_est1/(m-n)
    beta_est1 = beta_est1/(m-n)
    sigma_est1 = sigma_est1/(m-n)


    # print(paste0("Estimated parameter values monthly y"))
    # print(c(beta_est1, alpha_est1/beta_est1, sigma_est1, mean(y)))
    # hdi(est)
    #
    # print(paste0("Estimated parameter values monthly y"))
    # print(c(alpha_est1, beta_est1, sigma_est1, mean(y)))
    # hdi(est)

    mean(est[,1])
    mean(est[,2])
    mean(est[,3])


    sd(est[,1])
    sd(est[,2])
    sd(est[,3])

    estimate = c(alpha_est1, beta_est1, sigma_est1)
    sd1 = apply(est,2, sd)
    cf = t(hdi(est))
    df= data.frame(estimate, sd1, cf)
    colnames(df) = c("Est.", "Std.", "HPD CI(95%)_lower", "HPD CI(95%)_upper")

    df_simu_time_series = data.frame(est)
    colnames(df_simu_time_series) = c("alpha", "beta", "simga")




    L_y=length(y)
    L1=L_y -1
    ress=rep(0, L1)
    mu=rep(0, L1)
    vol=rep(0, L1)
    rr=rep(0, L1)

    aic=0
    bic=0
    c=0
    for (i in 2:L_y){
      mu[i-1] =y[i-1] + alpha_est1 -beta_est1*y[i-1]
      rr[i-1] = y[i]-mu[i-1]
      vol[i-1]= sigma_est1*y[i-1]
      ress[i-1] = rr[i-1]/vol[i-1]
      c=c -2*log(dnorm(rr[i-1], 0, vol[i-1]))
    }

    bic = c + 2 * log(L1)
    aic = c + 2 * 3


    model_cdf=rep(0, L1)
    for (i in 2:L_y){
      model_cdf[i] = pnorm(y[i], y[i-1] +alpha_est1 -beta_est1 *y[i-1], sigma_est1*y[i-1])
    }
    #ks.test(model_cdf, "punif",0,1)

    lt = list("estimate" = df, "standard_residuals"=ress , "mcmc_series" = df_simu_time_series, "Aic" = aic, "Bic"=bic, "empir_cdf"=model_cdf)

    #newlist <- list(lt,mn)

    return(lt)


  }


  # now we fit CIR model


  if (toupper(model_type) %in% c("CIR")){

    NN=length(y)

    alpha1=0.001
    beta1= 0.002
    sigma1=0.001

    alpha_est1=0
    beta_est1 =0
    sigma_est1 =0

    n = burn_in

    m = N

    est = matrix(rep(0, (m-n)*3), nrow=m-n)


    for (i in 1:m){
      #print(i)
      c=0
      d=0
      for (t in 1:(NN-1))
      {
        #c=c+1/y[t]^2
        c=c+1/y[t]
        d=d+(y[t+1]-(1-beta1)*y[t])/y[t]
      }

      a=d/c
      b=sigma1/sqrt(c)
      alpha1=rnorm(1, a, b)


      c=0
      d=0
      for (t in 1:(NN-1))
      {
        c=c+y[t]
        d=d+alpha1 +y[t] -y[t+1]
      }

      a=d/c
      b=sigma1/sqrt(c)
      beta1=rnorm(1, a, b)

      c=0
      for (t in 1:(NN-1))
      {
        c=c+ (y[t+1]-alpha1 -(1-beta1)*y[t])^2/y[t]
      }

      c=c/2

      sigma1=sqrt(rinvgamma(1,(NN / 2-1), c))

      #d1 <- data.frame(i, alpha1, beta1, sigma1)
      #write.table(d1, "cir-mcmc-R.txt", row.names = FALSE,
      #            col.names = FALSE, append = TRUE)

      #print(c(alpha1, beta1, sigma1))

      if (i >n){

        alpha_est1=alpha_est1 +alpha1
        beta_est1 = beta_est1 +beta1
        sigma_est1 = sigma_est1 +sigma1
        est[i-n,]=c(alpha1, beta1, sigma1)
      }
    }
    alpha_est1=alpha_est1/(m-n)
    beta_est1 = beta_est1/(m-n)
    sigma_est1 = sigma_est1/(m-n)


    # print(paste0("Estimated parameter values monthly y"))
    # print(c(beta_est1, alpha_est1/beta_est1, sigma_est1, mean(y)))
    # hdi(est)
    #
    # print(paste0("Estimated parameter values monthly y"))
    # print(c(alpha_est1, beta_est1, sigma_est1, mean(y)))
    # hdi(est)

    mean(est[,1])
    mean(est[,2])
    mean(est[,3])


    sd(est[,1])
    sd(est[,2])
    sd(est[,3])

    estimate = c(alpha_est1, beta_est1, sigma_est1)
    sd1 = apply(est,2, sd)
    cf = t(hdi(est))
    df= data.frame(estimate, sd1, cf)
    colnames(df) = c("Est.", "Std.", "HPD CI(95%)_lower", "HPD CI(95%)_upper")

    df_simu_time_series = data.frame(est)
    colnames(df_simu_time_series) = c("alpha", "beta", "simga")




    L_y=length(y)
    L1=L_y -1
    ress=rep(0, L1)
    mu=rep(0, L1)
    vol=rep(0, L1)
    rr=rep(0, L1)

    aic=0
    bic=0
    c=0
    for (i in 2:L_y){
      mu[i-1] =y[i-1] + alpha_est1 -beta_est1*y[i-1]
      rr[i-1] = y[i]-mu[i-1]
      vol[i-1]= sigma_est1 *sqrt(y[i-1])
      ress[i-1] = rr[i-1]/vol[i-1]
      c=c -2*log(dnorm(rr[i-1], 0, vol[i-1]))
    }

    bic = c + 2 * log(L1)
    aic = c + 2 * 3


    model_cdf=rep(0, L1)
    for (i in 2:L_y){
      model_cdf[i] = pnorm(y[i], y[i-1] +alpha_est1 -beta_est1 *y[i-1], sigma_est1*sqrt(y[i-1]))
      #ks.test(model_cdf, "punif",0,1)
    }
    lt = list("estimate" = df, "standard_residuals"=ress , "mcmc_series" = df_simu_time_series, "Aic" = aic, "Bic"=bic, "empir_cdf"=model_cdf)

    #newlist <- list(lt,mn)

    return(lt)


  }


  # now we fit CEV model


  if (toupper(model_type) %in% c("CEV")){

    NN=length(y)

    alpha1=0.001
    beta1= 0.002
    sigma1=0.001
    gama1=0.001

    alpha_est1=0
    beta_est1 =0
    sigma_est1 =0
    gama_est1 =0

    n = burn_in

    m = N

    est = matrix(rep(0, (m-n)*4), nrow=m-n)


    for (i in 1:m){
      #print(i)
      c=0
      d=0
      for (t in 1:(NN-1))
      {
        c=c+1/y[t]^(gama1*2)
        d=d+(y[t+1]-(1-beta1)*y[t])/y[t]^(gama1*2)
      }

      a=d/c
      b=sigma1/sqrt(c)
      alpha1=rnorm(1, a, b)


      c=0
      d=0
      for (t in 1:(NN-1))
      {
        c=c+y[t]^((1-gama1)*2)
        d=d+(alpha1 +y[t] -y[t+1]  )/y[t]^(gama1*2-1)
      }

      a=d/c
      b=sigma1/sqrt(c)
      beta1=rnorm(1, a, b)

      c=0
      for (t in 1:(NN-1))
      {
        c=c+ (y[t+1]-alpha1 -(1-beta1)*y[t])^2 /y[t]^(gama1*2)
      }

      c=c/2

      sigma1=sqrt(rinvgamma(1,(NN / 2-1), c))

      bmin=-10
      bmax=10

      for (t in 1:(NN-1))
      {
        u1=runif(1,0,1)*exp(-( y[t+1]-alpha1 -(1-beta1)*y[t])^2/2/(sigma1 *y[t]^gama1)^2 )
        cc1= (y[t+1]-alpha1 -(1-beta1)*y[t])^2/2/sigma1^2

        u2= runif(1,0,1)/y[t]^gama1

        if(log(y[t])>0){
          bmin=max(bmin, log(-cc1/log(u1))/(2*log(y[t])))
          bmax=min(bmax, -log(u2)/log(y[t]) )
        }

        if(log(y[t])<0){
          bmin=max(bmin, -log(u2)/log(y[t]))
          bmax=min(bmax, log(-cc1/log(u1))/(2*log(y[t])) )
        }

      }

      gama1=runif(1, bmin, bmax)



      if (i >n){
        alpha_est1=alpha_est1 +alpha1
        beta_est1 = beta_est1 +beta1
        sigma_est1 = sigma_est1 +sigma1
        gama_est1 = gama_est1 +gama1
        est[i-n,]=c(alpha1, beta1, sigma1, gama1)
      }
    }
    alpha_est1=alpha_est1/(m-n)
    beta_est1 = beta_est1/(m-n)
    sigma_est1 = sigma_est1/(m-n)
    gama_est1 = gama_est1/(m-n)


    mean(est[,1])
    mean(est[,2])
    mean(est[,3])


    sd(est[,1])
    sd(est[,2])
    sd(est[,3])

    estimate = c(alpha_est1, beta_est1, sigma_est1, gama_est1)
    sd1 = apply(est,2, sd)
    cf = t(hdi(est))
    df= data.frame(estimate, sd1, cf)
    colnames(df) = c("Est.", "Std.", "HPD CI(95%)_lower", "HPD CI(95%)_upper")

    df_simu_time_series = data.frame(est)
    colnames(df_simu_time_series) = c("alpha", "beta", "simga", "gamma")




    L_y=length(y)
    L1=L_y -1
    ress=rep(0, L1)
    mu=rep(0, L1)
    vol=rep(0, L1)
    rr=rep(0, L1)

    aic=0
    bic=0
    c=0
    for (i in 2:L_y){
      mu[i-1] =y[i-1] + alpha_est1 -beta_est1*y[i-1]
      rr[i-1] = y[i]-mu[i-1]
      vol[i-1]= sigma_est1*y[i-1]^gama_est1
      ress[i-1] = rr[i-1]/vol[i-1]
      c=c -2*log(dnorm(rr[i-1], 0, vol[i-1]))
    }

    bic = c + 2 * log(L1)
    aic = c + 2 * 4


    model_cdf=rep(0, L1)
    for (i in 2:L_y){
      model_cdf[i] = pnorm(y[i], y[i-1] +alpha_est1 -beta_est1 *y[i-1], sigma_est1*y[i-1]^gama_est1)
      #ks.test(model_cdf, "punif",0,1)
    }
    lt = list("estimate" = df, "standard_residuals"=ress , "mcmc_series" = df_simu_time_series, "Aic" = aic, "Bic"=bic, "empir_cdf"=model_cdf)

    #newlist <- list(lt,mn)

    return(lt)


  }
}
