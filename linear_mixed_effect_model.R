rm(list = ls())
library(dplyr)
library(brms)
library(msm)
library(furrr)

set.seed(123)

datasim = function(
    N, # num of clusters
    m, # the maximum number of observations per group
    p, # num of baseline covariates
    b_z, # coefficient of treatment
    b_ixyz, # coefficient of baseline covariates include intercept
    b_sig # standard_deviation of random intercept
    
){
  T = sample(1:m,N, replace = TRUE) # number of observations per group
  Nn = sum(T) # total number of observation
  x = sapply(1:p, function(r)rbinom(N,1,0.5)) # baseline covariates
  x = cbind(x, sapply(1:p, function(r)rnorm(N,1,0.3)))
  z = c() ## treatment indicator
  for(i in 1: length(T)){
    z = c(z,rbinom(T[i],1,0.5))
  }
  b = rnorm(N,0.5,b_sig)   ### random intercept
  time = c() ## treatment indicator
  for(i in 1: length(T)){
    time = c(time, sort(runif(T[i],3,20)))
  }
  time_z1 <- cos(0.5*pi * time / 10) + 1  ## effect of time on trt 1
  time_z0 <- 2 / (1 + 1.5 * exp(time/4))+0.5 ## effect of time on trt 0
  y0 = rtnorm(N,2,1, lower = 3) # baseline outcome
  y = numeric()
  y_pre = numeric()
  id = rep(1:N, times = T)
  c_id = 1
  for (i in 1:N){
    y_temp=numeric()
    for (j in 1:T[i]){
      if (j == 1) {
        y[c_id] = t(c(1,x[i,],y0[i],z[c_id]))%*% b_ixyz+
          time_z1[c_id]*z[c_id]+time_z0[c_id]*(1-z[c_id])+b[i]+ rnorm(1,0,0.2)
      } else {
        y[c_id] = t(c(1,x[i,],y[c_id-1],z[c_id]))%*% b_ixyz +
          time_z1[c_id]*z[c_id]+time_z0[c_id]*(1-z[c_id])+b[i] + rnorm(1,0,0.2) 
        # Subsequent outcomes based on previous
      }
      y_temp = c(y_temp, y[c_id])
      c_id = c_id+1
    }
    y_pre = c(y_pre,y0[i],y_temp[-length(y_temp)])
    
  }
  data = data.frame(id=rep(1:N, times = T), x1=rep(x[,1], times =T), 
                    x2=rep(x[,2], times =T), x3=rep(x[,3], times =T),
                    x4=rep(x[,4], times =T), time = time, 
                    trt = z, y_pre=y_pre, y = y )
  
  return(data)
}

### estimate treatment effect
est_y = function(fit,timepoint, prev_y,baseline){
  nobs = nrow(baseline)
  niter = nrow(prev_y)
  est = matrix(NA, niter, nobs)
  baseline$time = timepoint
  samples = as.matrix(fit)
  fixdata = model.matrix(~x1+x2+x3+x4+y_pre+trt, data = baseline)
  sm_int = posterior_epred(fit, newdata = baseline)-samples[,1:7]%*%t(fixdata)
  
  for(i in 1: niter){
    fixdata[,6] = c(prev_y[i,])
    est[i,] = fixdata%*%samples[i,1:7]+sm_int[i,]
  } 
  return(est)}

### priors setting
priors <- c(
  set_prior("normal(0, 1)", class = "b"), # For regression coefficients
  set_prior("normal(0, 1)", class = "Intercept"), # For the global intercept
  set_prior("cauchy(0, 1)", class = "sigma")
)


### function to fit data
lmm_model = function(tt){
  
  print(tt)
  data = datasim( N = 100, m = 10, p = 2, 
                  b_ixyz = c(0.2, 0.03, 0.06, 0.09, 0.12, 0.3, 0.5), b_sig = 0.25)
  data$trt=factor(data$trt)
  
  model1 = brm(y~1+x1+x2+x3+x4+y_pre+trt+trt:time+ (1|id), family = gaussian(), 
               data=data, prior = priors, chains = 2, iter = 60000, thin=10,silent = TRUE,
               warmup = 15000,control = list(adapt_delta = 0.99, max_treedepth = 5))
  
  ### treatment = 1
  trt1 = data
  trt1$trt = 1
  trt1$trt = factor(trt1$trt, levels = levels(factor(data$trt)))
  
  baseline1 = trt1%>%
    group_by(id) %>%      
    arrange(time, .by_group = TRUE) %>%  
    slice(1) %>%
    ungroup()%>%
    mutate(time = 6)
  
  y6_1 = posterior_epred(model1, newdata = baseline1)
  y12_1 = est_y(model1, 12, y6_1, baseline1)
  y18_1 = est_y(model1, 18, y12_1, baseline1)
  
  ### treatment = 0
  trt0 = data
  trt0$trt = 0
  trt0$trt = factor(trt0$trt, levels = levels(factor(data$trt)))
  
  baseline0 = trt0%>%
    group_by(id) %>%      
    arrange(time, .by_group = TRUE) %>%  
    slice(1) %>%
    ungroup()%>%
    mutate(time = 6)
  
  y6_0 = posterior_epred(model1, newdata = baseline0)
  y12_0 = est_y(model1, 12, y6_0, baseline0)
  y18_0 = est_y(model1, 18, y12_0, baseline0)
  
  pred_data = data.frame(y6_0 = rowMeans(y6_0), y12_0 = rowMeans(y12_0), y18_0 = rowMeans(y18_0),
                         y6_1 = rowMeans(y6_1), y12_1 = rowMeans(y12_1), y18_1 = rowMeans(y18_1))
  
  return(pred_data)
}



lmm <- list()
for (i in 1:20) {
  plan(multisession, workers = 25)
  lmm[[i]] <-
    future_map(as.list((1 + (i - 1) * 25):(25 + (i - 1) * 25)), lmm_model)
}

