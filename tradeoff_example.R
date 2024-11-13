library(xtable)

make_data <- function(n=10, trt_prop=.5, d_prop_trt=.1, d_prop_ctrl=.2, d_x=0, dtx_y=0, dx_y=0){
  max.time <- 2
  trt.time <- 2

  dat <- expand.grid(id = c(1:n), tp = 1:max.time) %>% arrange(id,tp) %>% group_by(id) %>%
    mutate(int=rnorm(1,0,sd=0.25), # random intercept
           p.trt=trt_prop, # probability of treatment
           trt=rbinom(1, 1, p.trt), # treatment
           p.d.trt=d_prop_trt, # probability in subgroup given treated
           p.d.ctrl=d_prop_ctrl, # probability in subgroup given not treated
           d=rbinom(1,1, p.d.trt*trt + p.d.ctrl*(1-trt)), # subgroup indicator
           x=rnorm(1, mean = 1 - 0.5*trt + .2*d*d_x +.1*trt*d*d_x, sd = 1),
           post=I(tp >= trt.time), # indicator of post-treatment period
           treated=I(post == 1 & trt == 1) # time-varying indicator if treated or not
    ) %>% 
    ungroup()
  
  dat <- dat %>% mutate(
    y0 = 1 + x*(tp+.2*d*dx_y+.1*d*tp*dtx_y) + trt + int + ((tp - 2.5)^2)/10,
    y1 = 1 + x*(tp+.2*d*dx_y+.1*d*tp*dtx_y) + trt + int + (1-.2*d) + ((tp - 2.5)^2)/10,
    y = 1 + x*(tp+.2*d*dx_y+.1*d*tp*dtx_y) + trt + int + treated*(1-.2*d) + ((tp - 2.5)^2)/10) %>%
    group_by(id) %>% mutate(y.diff = y - lag(y)) %>% ungroup()
  
  return(dat)
}

params <- data.frame(d_x=  c(0,1,1,1),
                     dtx_y=c(0,0,0,1),
                     dx_y= c(0,0,1,1))

results <- numeric(0)

for(i in 1:dim(params)[1]){
  newdat <- do.call(make_data,c(list(n=200, trt_prop=.5, d_prop_trt=.1, d_prop_ctrl=.2),
                                as.list(params[i,])))
  
  # 1
  twfe_1 <- lm(y ~ trt*post*d, newdat)
  # summary(twfe_1)$coefficients[c("trt:postTRUE","trt:postTRUE:d"),]
  
  # 2
  twfe_2 <- lm(y ~ trt*post*d + x*post, newdat)
  # summary(twfe_2)$coefficients[c("trt:postTRUE","trt:postTRUE:d"),]
  
  # 3
  twfe_3 <- lm(y ~ trt*post*d + x*post + x*d, newdat)
  # summary(twfe_3)$coefficients[c("trt:postTRUE","trt:postTRUE:d"),]
  # Note: adjusting by x*d, similar to adjusting by x, does nothing for confounding
  # important part of confounding is the time-varying nature of the covariate? 
  
  # 4
  twfe_4 <- lm(y ~ trt*post*d + x*post*d, newdat)
  # summary(twfe_4)$coefficients[c("trt:postTRUE","trt:postTRUE:d"),]
  
  results <- rbind(results,c(summary(twfe_1)$coefficients["trt:postTRUE:d","Estimate"], 
    summary(twfe_2)$coefficients["trt:postTRUE:d","Estimate"],
    summary(twfe_3)$coefficients["trt:postTRUE:d","Estimate"],
    summary(twfe_4)$coefficients["trt:postTRUE:d","Estimate"])
  )
}

print(xtable(results))


# Error/bias from sample size
props <- seq(.1,.7,.05)
stderrs <- rep(NA, length(props))
ests <- rep(NA, length(props))

for(i in 1:length(props)){
  newdat <- do.call(make_data,c(list(n=200, trt_prop=.5, d_prop_trt=props[i], d_prop_ctrl=props[i]),
                                as.list(params[4,])))
  twfe_4 <- lm(y ~ trt*post*d + x*post*d, newdat)
  stderrs[i] <- summary(twfe_4)$coefficients["trt:postTRUE:d","Std. Error"]
  ests[i] <- summary(twfe_4)$coefficients["trt:postTRUE:d","Estimate"]
}

plot(props, stderrs, xlab = "Proportion in subgroup D", ylab = "Standard Error of trt:post:d estimate")
plot(props, ests, xlab = "Proportion in subgroup D", ylab = "trt:post:d estimate")


# Model misspecification
dtx_coefs <- seq(0,10,1)
stderrs <- rep(NA, length(dtx_coefs))
ests <- rep(NA, length(dtx_coefs))

for(i in 1:length(dtx_coefs)){
  newdat <- do.call(make_data, list(n=200, trt_prop=.5, d_prop_trt=.2, d_prop_ctrl=.2, 
                                    d_x=1, dtx_y=dtx_coefs[i], dx_y = 1) )
  twfe_2 <- lm(y ~ trt*post*d + x*post, newdat)
  stderrs[i] <- summary(twfe_2)$coefficients["trt:postTRUE:d","Std. Error"]
  ests[i] <- summary(twfe_2)$coefficients["trt:postTRUE:d","Estimate"]
}

plot(dtx_coefs, stderrs, xlab = "True dtx coefficient", ylab = "Standard Error of trt:post:d estimate")
plot(dtx_coefs, ests, xlab = "True dtx coefficient", ylab = "trt:post:d estimate")




