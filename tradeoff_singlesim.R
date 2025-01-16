library(xtable)
library(ggplot2)
library(broom)
library(modelr)

# Error/bias from sample size - original version (no replications) -----------------------------------------------

props <- seq(.1,.7,.05)
stderrs4 <- rep(NA, length(props))
stderrs2 <- rep(NA, length(props))
ests4 <- rep(NA, length(props))
ests2 <- rep(NA, length(props))
rmses4 <- rep(NA, length(props))
rmses2 <- rep(NA, length(props))

for(i in 1:length(props)){
  newdat <- do.call(make_data,c(list(n=200, trt_prop=.5, d_prop_trt=props[i], d_prop_ctrl=props[i]),
                                as.list(params[4,])))
  twfe_4 <- lm(y ~ trt*post*d + x*post*d, newdat)
  stderrs4[i] <- summary(twfe_4)$coefficients["trt:postTRUE:d","Std. Error"]
  ests4[i] <- summary(twfe_4)$coefficients["trt:postTRUE:d","Estimate"]
  rmses4[i] <- sqrt(mean(twfe_4$residuals^2))
  
  # misspecified model
  twfe_2 <- lm(y ~ trt*post*d + x*post, newdat)
  stderrs2[i] <- summary(twfe_2)$coefficients["trt:postTRUE:d","Std. Error"]
  ests2[i] <- summary(twfe_2)$coefficients["trt:postTRUE:d","Estimate"]
  rmses2[i] <- sqrt(mean(twfe_2$residuals^2))
}

bias4 <- ests4 + .2
bias2 <- ests2 + .2

plot(props, stderrs4, xlab = "Proportion in subgroup D", ylab = "Standard Error of trt:post:d estimate",col=3)
points(props, stderrs2,col=2)
legend('topright', legend = c('Misspecified model', 'Correct model'), col = c(2, 3), lty = 1)

plot(props, bias2, xlab = "Proportion in subgroup D", ylab = "trt:post:d bias",col=2)
plot(props, bias4,col=3)
legend('topright', legend = c('Misspecified model', 'Correct model'), col = c(2, 3), lty = 1)

plot(props, rmses4, xlab = "Proportion in subgroup D", ylab = "trt:post:d RMSE",ylim = c(.22,.32), col=3)
points(props, rmses2,col=2)
legend('topright', legend = c('Misspecified model', 'Correct model'), col = c(2, 3), lty = 1)




# Model misspecification --------------------------------------------------

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
