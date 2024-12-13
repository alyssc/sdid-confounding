library(xtable)
library(ggplot2)
library(broom)
library(modelr)

make_data <- function(n=100, trt_prop=.5, d_prop_trt=.1, d_prop_ctrl=.2, d_x=1, dtx_y=1, dx_y=1, noise=1, trt_effect_noise=.5){
  max.time <- 2
  trt.time <- 2
  
  # number of units in treated, treated&d, and ctrl&d respectively
  trt_n <- round(trt_prop*n)
  d_n_trt <- round( round(trt_prop*n)*d_prop_trt )
  d_n_ctrl <- round( round((1-trt_prop)*n)*d_prop_ctrl )
  # print(c(trt_n,d_n_trt,d_n_ctrl))

  dat <- expand.grid(id = c(1:n), tp = 1:max.time) %>% arrange(id,tp) %>% group_by(id) %>%
    mutate(int=rnorm(1,0,sd=noise), # random intercept
           trt=1*I(id <= trt_n), # treatment
           d=1*I(id <= d_n_trt | (id > trt_n & id <= trt_n+d_n_ctrl) ), # subgroup indicator
           x=rnorm(1, mean = 1 - 0.5*trt + .2*d*d_x +.1*trt*d*d_x, sd = 1),
           post=I(tp >= trt.time), # indicator of post-treatment period
           treated=I(post == 1 & trt == 1) # time-varying indicator if treated or not
    ) %>% 
    ungroup()
  
  dat <- dat %>% mutate(
    d_trt_effect = rnorm(1,-.2,sd=trt_effect_noise),
    y0 = 1 + x*(tp+.2*d*dx_y+.1*d*tp*dtx_y) + trt + int + ((tp - 2.5)^2)/10,
    y1 = y0 + (1+d_trt_effect*d),
    y = treated*y1+(1-treated)*y0  ) %>%
    group_by(id) %>% mutate(y.diff = y - lag(y)) %>% ungroup()
  
  return(dat)
}

# make_data with df of parameters as input
make_data_df <- function(params){
  return( do.call(make_data, as.list(params) ) )
}


# DGMs x regression models table ------------------------------------------

params <- data.frame(d_x=  c(0,1,1,1),
                     dtx_y=c(0,0,0,1),
                     dx_y= c(0,0,1,1))

results <- numeric(0)


for(i in 1:dim(params)[1]){
  newdat <- do.call(make_data,c(list(n=200),
                                as.list(params[i,])))
  twfe_1 <- lm(y ~ trt*post*d, newdat)
  # summary(twfe_1)$coefficients[c("trt:postTRUE","trt:postTRUE:d"),]
  twfe_2 <- lm(y ~ trt*post*d + x*post, newdat)
  # summary(twfe_2)$coefficients[c("trt:postTRUE","trt:postTRUE:d"),]
  twfe_3 <- lm(y ~ trt*post*d + x*post + x*d, newdat)
  # summary(twfe_3)$coefficients[c("trt:postTRUE","trt:postTRUE:d"),]
  # Note: adjusting by x*d, similar to adjusting by x, does nothing for confounding
  # important part of confounding is the time-varying nature of the covariate? 
  twfe_4 <- lm(y ~ trt*post*d + x*post*d, newdat)
  # summary(twfe_4)$coefficients[c("trt:postTRUE","trt:postTRUE:d"),]
  
  results <- rbind(results,c(summary(twfe_1)$coefficients["trt:postTRUE:d","Estimate"], 
    summary(twfe_2)$coefficients["trt:postTRUE:d","Estimate"],
    summary(twfe_3)$coefficients["trt:postTRUE:d","Estimate"],
    summary(twfe_4)$coefficients["trt:postTRUE:d","Estimate"])
  )
}

print(xtable(results))



# Functions to calculate error/bias across params ---------------------------------------------

true_satt <- function(data){
  att <- data %>% 
    filter(trt==1, tp==1) %>%
    summarize(mean = mean(d_trt_effect))
  return(att$mean[[1]])
}

simanalyze <- function(params){
  sim_data <- params %>% 
    group_by(scenario,replicate) %>%
    nest() %>% mutate(simdat=map(data,make_data_df)) %>%
    unnest(cols=c(simdat, data))
  
  ests <- sim_data %>% 
    group_by(scenario,replicate) %>%
    nest() %>% 
    mutate(true_satt=map_dbl(data, true_satt)) %>% 
    crossing(misspec = c(0,1)) %>%
    mutate(model = 
             case_when(misspec == 0 ~ map(data, ~lm(y ~ trt*post*d+x*post*d, data = .) ),
                       misspec == 1 ~ map(data, ~lm(y ~ trt*post*d+x*post, data = .) )
             )) %>% 
    mutate(coefs = map(model, broom::tidy)) %>%
    unnest(coefs) %>% 
    filter(term == "trt:postTRUE:d") %>% 
    select(-c(statistic, p.value, term, data, model)) %>%
    mutate(bias = estimate-true_satt)
  
  return(list('data'=sim_data,'ests'=ests))
}


# Error across sample size ------------------------------------------------

NREP <- 100

vary_d_prop <- data.frame(expand.grid(
  n=100, 
  d_prop_trt=seq(.02,.2,.02), 
  d_x=1, dtx_y=1, dx_y=1, noise=1,trt_effect_noise=0,
  replicate = 1:NREP)) %>%
  mutate(d_prop_ctrl=d_prop_trt,
         scenario=rep(1:(n()/NREP), NREP)) 

simdat <- simanalyze(vary_d_prop)

scenario_to_count <- simdat$data %>% 
  ungroup() %>%
  filter(replicate ==1, post==0) %>%
  select(scenario, d) %>%
  group_by(scenario) %>%
  summarize(d_count = sum(d))

ests <- simdat$ests %>%
  inner_join(scenario_to_count, by="scenario")

ggplot(ests, aes(x=factor(d_count), y=bias, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Bias across subgroup size, nrep=%s, n=100", NREP),
       x="# of units in subgroup D", y="Bias of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/bias_vary_d_prop.png",width=6,height=4)

ggplot(ests, aes(x=factor(d_count), y=std.error, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Std error across subgroup size, nrep=%s, n=100", NREP),
       x="# of units in subgroup D", y="Std error of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/std_error_vary_d_prop.png",width=6,height=4)


# Error across sample size; plus heterogenous treatment effect

vary_d_prop <- data.frame(expand.grid(
  n=100, 
  d_prop_trt=seq(.02,.2,.02), 
  d_x=1, dtx_y=1, dx_y=1, noise=1,trt_effect_noise=0,
  replicate = 1:NREP)) %>%
  mutate(d_prop_ctrl=d_prop_trt,
         scenario=rep(1:(n()/NREP), NREP)) 

simdat <- simanalyze(vary_d_prop)

scenario_to_count <- simdat$data %>% 
  ungroup() %>%
  filter(replicate ==1, post==0) %>%
  select(scenario, d) %>%
  group_by(scenario) %>%
  summarize(d_count = sum(d))

ests <- simdat$ests %>%
  inner_join(scenario_to_count, by="scenario")

ggplot(ests, aes(x=factor(d_count), y=bias, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Bias across subgroup size, nrep=%s, n=100", NREP),
       x="# of units in subgroup D", y="Bias of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/bias_vary_d_prop_hetero_effect.png",width=6,height=4)

ggplot(ests, aes(x=factor(d_count), y=std.error, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Std error across subgroup size, nrep=%s, n=100", NREP),
       x="# of units in subgroup D", y="Bias of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/std_error_vary_d_prop_hetero_effect.png",width=6,height=4)


# Varying D control group size only

NREP=100
vary_d_ctrl <- data.frame(expand.grid(
  n=100, trt_prop=.5, 
  d_prop_trt=.02,
  d_prop_ctrl=seq(.02,.3,.05),
  d_x=1, dtx_y=1, dx_y=1,noise=1,
  replicate = 1:NREP)) %>%
  mutate(scenario=rep(1:(n()/NREP), NREP)) 

simdat <- simanalyze(vary_d_ctrl)
vals = seq(.02,.3,.05)
scenario_vals <- data.frame(
  vals = vals,
  scenario = 1:length(vals)
)

ests <- simdat$ests  %>%
  inner_join(scenario_vals, by="scenario")

ggplot(ests, aes(x=factor(vals), y=bias, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Keeping prop_d_trt=.02 constant, nrep=%s, n=100", NREP),
       x="prop_d_ctrl", y="Bias of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/bias_vary_d_ctrl.png",width=6,height=4)

ggplot(ests, aes(x=factor(vals), y=std.error, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Keeping prop_d_trt=.02 constant, nrep=%s, n=100", NREP),
       x="prop_d_ctrl", y="Std error of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/std_error_vary_d_ctrl.png",width=6,height=4)




NREP=100
vary_d_ctrl <- data.frame(expand.grid(
  n=100, trt_prop=.5, 
  d_prop_trt=.05,
  d_prop_ctrl=seq(.02,.3,.05),
  d_x=1, dtx_y=1, dx_y=1,noise=1,
  replicate = 1:NREP)) %>%
  mutate(scenario=rep(1:(n()/NREP), NREP)) 

simdat <- simanalyze(vary_d_ctrl)
vals = seq(.02,.3,.05)
scenario_vals <- data.frame(
  vals = vals,
  scenario = 1:length(vals)
)

ests <- simdat$ests  %>%
  inner_join(scenario_vals, by="scenario")

ggplot(ests, aes(x=factor(vals), y=bias, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Keeping prop_d_trt=.05 constant, nrep=%s, n=100", NREP),
       x="prop_d_ctrl", y="Bias of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/bias_vary_d_ctrl.png",width=6,height=4)

ggplot(ests, aes(x=factor(vals), y=std.error, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Keeping prop_d_trt=.05 constant, nrep=%s, n=100", NREP),
       x="prop_d_ctrl", y="Std error of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/std_error_vary_d_ctrl.png",width=6,height=4)


# Does increasing noise do anything?
NREP=100
vary_noise <- data.frame(expand.grid(
  n=100, trt_prop=.5, 
  d_prop_trt=.04,
  d_prop_ctrl=.04,
  d_x=1, dtx_y=1, dx_y=1,noise=seq(1,10,1),
  replicate = 1:NREP)) %>%
  mutate(scenario=rep(1:(n()/NREP), NREP)) 

simdat <- simanalyze(vary_noise)
vals <- seq(1,10,1)
scenario_vals <- data.frame(
  vals = vals,
  scenario = 1:length(vals)
)

ests <- simdat$ests  %>%
  inner_join(scenario_vals, by="scenario")

ggplot(ests, aes(x=factor(vals), y=bias, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Adding noise, prop_d_trt/ctrl=.04, nrep=%s, n=100", NREP),
       x="SD of unit-level noise", y="Bias of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/bias_varynoise.png",width=6,height=4)


# What about increasing treatment effect heterogeneity to d-specific trt effect? 

NREP=50
vary_trteffectnoise <- data.frame(expand.grid(
  n=100, trt_prop=.5, 
  d_prop_trt=.1,
  d_prop_ctrl=.1,
  d_x=1, dtx_y=1, dx_y=1,noise=1,
  trt_effect_noise = seq(1,10,1),
  replicate = 1:NREP)) %>%
  mutate(scenario=rep(1:(n()/NREP), NREP)) 

simdat <- simanalyze(vary_trteffectnoise)
vals <- seq(1,10,1)
scenario_vals <- data.frame(
  vals = vals,
  scenario = 1:length(vals)
)

ests <- simdat$ests  %>%
  inner_join(scenario_vals, by="scenario")

ggplot(ests, aes(x=factor(vals), y=bias, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Adding d-trt effect heterogeneity, prop_d_trt/ctrl=.04, nrep=%s, n=100", NREP),
       x="SD of unit-level d-trt effect", y="Bias of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/bias_varytrteffect.png",width=6,height=4)


# Multiple covariates ----------------------------------------------------

make_data_multX <- function(n=100, trt_prop=.5, d_prop_trt=.1, d_prop_ctrl=.2, noise=1, trt_effect_noise=.5){
  max.time <- 2
  trt.time <- 2
  
  # number of units in treated, treated&d, and ctrl&d respectively
  trt_n <- round(trt_prop*n)
  d_n_trt <- round( round(trt_prop*n)*d_prop_trt )
  d_n_ctrl <- round( round((1-trt_prop)*n)*d_prop_ctrl )
  # print(c(trt_n,d_n_trt,d_n_ctrl))
  
  dat <- expand.grid(id = c(1:n), tp = 1:max.time) %>% arrange(id,tp) %>% group_by(id) %>%
    mutate(int=rnorm(1,0,sd=noise), # random intercept
           trt=1*I(id <= trt_n), # treatment
           d=1*I(id <= d_n_trt | (id > trt_n & id <= trt_n+d_n_ctrl) ), # subgroup indicator
           x1=rnorm(1, mean = 1 - 0.5*trt, sd = 1),
           x2=rnorm(1, mean = 1.2 - 0.6*trt + .2*d +.1*trt*d, sd = 1),
           x3=rnorm(1, mean = 1.2 - 0.6*trt + .2*d +.1*trt*d, sd = 1),
           x4=rnorm(1, mean = 1.2 - 0.6*trt + .2*d +.1*trt*d, sd = 1),
           x5=rnorm(1, mean = 1.2 - 0.6*trt + 1*d +1*trt*d, sd = 1),
           post=I(tp >= trt.time), # indicator of post-treatment period
           treated=I(post == 1 & trt == 1) # time-varying indicator if treated or not
    ) %>% 
    ungroup()
  
  dat <- dat %>% mutate(
    d_trt_effect = rnorm(1,-.2,sd=trt_effect_noise),
    y0 = 1 + x1*tp + x2*tp + x3*(tp +.2*d+.1*d*tp)+  3*x4*(tp+.2*d+.1*d*tp) + x5*(tp+.2*d+.1*d*tp) + trt + int + ((tp - 2.5)^2)/10,
    y1 = y0 + (1+d_trt_effect*d),
    y = treated*y1+(1-treated)*y0  ) %>%
    group_by(id) %>% mutate(y.diff = y - lag(y)) %>% ungroup()
  
  return(dat)
}

# make_data with df of parameters as input
make_data_df_multX <- function(params){
  return( do.call(make_data_multX, as.list(params) ) )
}

data <- make_data_multX()

small <- lm(y~trt*post*d+x1+x2+x3+x4+x5, 
            data=data)
correct <- lm(y~trt*post*d+x1*post+x2*post+x3*post+x4*post*d+x5*post*d, 
              data=data)
full <- lm(y~trt*post*d+x1*post*d+x2*post*d+x3*post*d+x4*post*d+x5*post*d, 
                     data=data)

NREP=100
params <- data.frame(expand.grid(
  d_prop_trt=.1,
  d_prop_ctrl=.1,
  replicate = 1:NREP))


simanalyze_multX <- function(params){
  sim_data <- params %>% 
    group_by(replicate) %>%
    nest() %>% mutate(simdat=map(data,make_data_df_multX)) %>%
    unnest(cols=c(simdat, data))
  
  ests <- sim_data %>% 
    group_by(replicate) %>%
    nest() %>% 
    mutate(true_satt=map_dbl(data, true_satt)) %>% 
    crossing(model_type = c("smallest","smaller_x3","smaller_x4","smaller_x5", "correct", "fuller1","fuller2","fuller3", "fullest")) %>%
    mutate(model_type = factor(model_type, levels=c("smallest","smaller_x3","smaller_x4","smaller_x5", "correct", "fuller1","fuller2","fuller3", "fullest"))) %>% 
    mutate(model = 
             case_when(
                       model_type == "smallest" ~ map(data, ~lm(y ~ trt*post*d+x1*post+x2*post+x3*post+x4*post+x5*post, data = .) ),
                       model_type == "smaller_x3" ~ map(data, ~lm(y ~ trt*post*d+x1*post+x2*post+x3*post*d+x4*post+x5*post, data = .) ),
                       model_type == "smaller_x4" ~ map(data, ~lm(y ~ trt*post*d+x1*post+x2*post+x3*post+x4*post*d+x5*post, data = .) ),
                       model_type == "smaller_x5" ~ map(data, ~lm(y ~ trt*post*d+x1*post+x2*post+x3*post+x4*post+x5*post*d, data = .) ),
                       model_type == "correct" ~ map(data, ~lm(y ~ trt*post*d+x1*post+x2*post+x3*post*d+x4*post*d+x5*post*d, data = .) ),
                       model_type == "fuller1" ~ map(data, ~lm(y ~ trt*post*d+x1*post*d+x2*post+x3*post*d+x4*post*d+x5*post*d, data = .) ),
                       model_type == "fuller2" ~ map(data, ~lm(y ~ trt*post*d+x1*post+x2*post*d+x3*post*d+x4*post*d+x5*post*d, data = .) ),
                       model_type == "fullest" ~ map(data, ~lm(y ~ trt*post*d+x1*post*d+x2*post*d+x3*post*d+x4*post*d+x5*post*d, data = .) )
             )) %>% 
    mutate(coefs = map(model, broom::tidy)) %>%
    unnest(coefs) %>% 
    filter(term == "trt:postTRUE:d") %>% 
    select(-c(statistic, p.value, term, data, model)) %>%
    mutate(bias = estimate-true_satt)
  
  return(list('data'=sim_data,'ests'=ests))
}


simdat <- simanalyze_multX(params)

ests <- simdat$ests

ggplot(ests, aes(x=model_type, y=std.error)) + 
  geom_boxplot() +
  labs(title = "Standard error of misspecified and correct models",
       x="Model Type", y="Standard Error") 

ggsave("plots/mult_varyX_se.png",width=6,height=4)


ggplot(ests, aes(x=model_type, y=bias)) + 
  geom_boxplot() +
  labs(title = "Bias of misspecified and correct models",
       x="Model Type", y="Bias") 

ggsave("plots/mult_varyX_bias.png",width=6,height=4)

means <- ests %>% 
  group_by(model_type) %>%
  summarize(mean_bias = round(mean(bias),5), mean_se = round(mean(std.error),5))
print(xtable(means),include.rownames=FALSE)

means

pre_data <- data %>% 
  filter(post==FALSE) 
lm(y~x3,pre_data)
lm(y~x4,pre_data)
lm(y~x5,pre_data)





# Error/bias from sample size - original version -----------------------------------------------

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





