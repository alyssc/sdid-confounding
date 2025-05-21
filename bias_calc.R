# Can we calculate the bias theoretically, based on parameters, 
# in order to identify what contributes to subgroup-specific confounding? 

source("tradeoff_fn.R", local = knitr::knit_global())



# Decrease noise (cannot remove noise -- error)
misspec_bias <- vary_param_plot(param_name = "alpha_1", 
                                param_vals = seq(.1,1,.1),
                                param_defaults = data.frame(noise=.001,x_noise=.001),
                                NREP = 1)

empirical_bias <- misspec_bias$mean_bias

# theoretical bias
alpha_1 <- seq(.1,1,.1)
alpha_3 <- .1
beta_3 <- .1
theoretical_bias <- beta_3*(alpha_1+alpha_3)

comp_bias <- data.frame(empirical = empirical_bias, theoretical = theoretical_bias,
                        difference = empirical_bias-theoretical_bias)

xtable(comp_bias, digits = 4)


# Exploration: Why are empirical and theoretical results different? 
nreps <- 20
emp.biases <- rep(NA, nreps)
int.biases <- rep(NA, nreps) # intermediate result. uses gammas which are empirical
gamma_3s <- rep(NA, nreps)
gamma_1s <- rep(NA, nreps)
gamma_1misspecs <- rep(NA, nreps)

for(i in c(1:nreps)){
  alpha_1=-.5
  alpha_3=.1
  beta_3=.1
  data <- make_data(alpha_1=alpha_1,alpha_3=alpha_3,beta_3=beta_3)
  correct <- lm(y ~ trt*post*d+x*post*d, data = data)
  misspec <- lm(y ~ trt*post*d+x*post, data = data)
  gamma_1_mis <- misspec$coefficients["postTRUE:x"]
  gamma_1 <- correct$coefficients["postTRUE:x"]
  gamma_3 <- correct$coefficients["postTRUE:d:x"]
  
  gamma_1misspecs[i] <- gamma_1_mis
  gamma_1s[i] <- gamma_1
  gamma_3s[i] <- gamma_3
  
  d_ATT_mis <- misspec$coefficients["trt:postTRUE"]+misspec$coefficients["trt:postTRUE:d"]
  d_ATT_cor <- correct$coefficients["trt:postTRUE"]+correct$coefficients["trt:postTRUE:d"]
  emp.biases[i] <- d_ATT_mis-d_ATT_cor
  
  # theoretical bias 1 is not obtainable from the data
  int.biases[i] <- -(gamma_1_mis-gamma_1-gamma_3) *(alpha_1+alpha_3)
  th.bias <- (beta_3) * (alpha_1+alpha_3)
}


hist(th.bias-int.biases)
hist(emp.biases-th.bias)
hist(emp.biases-int.biases, breaks=10) #centered at 0, which is good

mean(emp.biases)
mean(th.bias-int.biases)
mean(emp.biases-th.bias)
mean(emp.biases-int.biases, breaks=10) 

summary(gamma_3s) #correctly specified model always correct gets beta.3
summary(gamma_1s) #correctly specified model always correct gets beta.1
summary(gamma_1misspecs) 
#misspecified model overshoots, which makes sense bc true answer is between beta.1 and beta.3

# the problem is that gamma_1mispec is not gamma_1
