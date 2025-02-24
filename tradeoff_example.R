library(xtable)
library(ggplot2)
library(broom)
library(modelr)
library(tidyverse)
library(margins)


source("tradeoff_fn.R", local = knitr::knit_global())


# DGMs x regression models table ------------------------------------------

params <- data.frame(d_x=  c(0,1,1,1),
                     xdt_y=c(0,0,0,1),
                     xd_y= c(0,0,1,1))

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


# Error across sample size ------------------------------------------------

NREP <- 100

vary_d_prop <- data.frame(expand.grid(
  n=100, 
  d_prop_trt=seq(.02,.2,.02), 
  d_x=1, xdt_y=1, xd_y=1, noise=1,trt_effect_noise=0,
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
  d_x=1, xdt_y=1, xd_y=1, noise=1,trt_effect_noise=0,
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
  d_x=1, xdt_y=1, xd_y=1,noise=1,
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
  d_x=1, xdt_y=1, xd_y=1,noise=1,
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
  d_x=1, xdt_y=1, xd_y=1,noise=seq(1,10,1),
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
  d_x=1, xdt_y=1, xd_y=1,noise=1,
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

# How does bias from model misspecification change with changing X~d?

NREP=50
vary_dx <- data.frame(expand.grid(
  n=100, trt_prop=.5, 
  d_prop_trt=.1,
  d_prop_ctrl=.1,
  d_x=seq(0,5,1), xdt_y=1, xd_y=1,
  replicate = 1:NREP)) %>%
  mutate(scenario=rep(1:(n()/NREP), NREP)) 

simdat <- simanalyze(vary_dx)
vals <- seq(0,5,1)
scenario_vals <- data.frame(
  vals = vals,
  scenario = 1:length(vals)
)

ests <- simdat$ests  %>%
  inner_join(scenario_vals, by="scenario")

ggplot(ests, aes(x=factor(vals), y=bias, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Increasing relationship X~d", NREP),
       x="d_x", y="Bias of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/bias_varydx.png",width=6,height=4)

# How does bias from model misspecification change with changing X~d*trt?

NREP=50
vary_dgx <- data.frame(expand.grid(
  n=100, trt_prop=.5, 
  d_prop_trt=.1,
  d_prop_ctrl=.1,
  dg_x=seq(0,8,1), xdt_y=1, xd_y=1,
  replicate = 1:NREP)) %>%
  mutate(scenario=rep(1:(n()/NREP), NREP)) 

simdat <- simanalyze(vary_dgx)
vals <- seq(0,8,1)
scenario_vals <- data.frame(
  vals = vals,
  scenario = 1:length(vals)
)

ests <- simdat$ests  %>%
  inner_join(scenario_vals, by="scenario")

ggplot(ests, aes(x=factor(vals), y=bias, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Increasing relationship X~dg", NREP),
       x="dg_x", y="Bias of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/bias_vary_dg_x.png",width=6,height=4)

# How does bias from model misspecification change with changing Y~Xt?

NREP=50
vary_xy <- data.frame(expand.grid(
  n=100, trt_prop=.5, 
  d_prop_trt=.1,
  d_prop_ctrl=.1,
  xt_y=seq(1,6), xdt_y=1, xd_y=1,
  replicate = 1:NREP)) %>%
  mutate(scenario=rep(1:(n()/NREP), NREP)) 

simdat <- simanalyze(vary_xy)
vals <- seq(0,6,1)
scenario_vals <- data.frame(
  vals = vals,
  scenario = 1:length(vals)
)

ests <- simdat$ests  %>%
  inner_join(scenario_vals, by="scenario")

ggplot(ests, aes(x=factor(vals), y=bias, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Increasing relationship Y~xt", NREP),
       x="x_y", y="Bias of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/bias_vary_xty.png",width=6,height=4)



NREP=50
vary_xy <- data.frame(expand.grid(
  n=100, trt_prop=.5, 
  d_prop_trt=.1,
  d_prop_ctrl=.1,
  xdt_y=seq(1,6), xt_y=1, xd_y=1,
  replicate = 1:NREP)) %>%
  mutate(scenario=rep(1:(n()/NREP), NREP)) 

simdat <- simanalyze(vary_xy)
vals <- seq(0,6,1)
scenario_vals <- data.frame(
  vals = vals,
  scenario = 1:length(vals)
)

ests <- simdat$ests  %>%
  inner_join(scenario_vals, by="scenario")

ggplot(ests, aes(x=factor(vals), y=bias, fill=factor(misspec))) + 
  geom_boxplot() +
  labs(title = sprintf("Increasing relationship Y~dtx", NREP),
       x="x_dtx", y="Bias of subgroup-specific estimate") +
  scale_fill_discrete(name = "Model type", labels = c("TWFE4 (correct)","TWFE2 (misspecified)"))

ggsave("plots/bias_vary_xdt_y.png",width=6,height=4)

# Decrease noise
# Error if we remove noise (perfect fit)
misspec_bias <- vary_param_plot(param_name = "alpha_1", 
                  param_vals = seq(.1,1,.1),
                  param_defaults = data.frame(noise=.001,x_noise=.001),
                  NREP = 50)

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

data <- make_data()
# defaults: 
# alpha_0=1, alpha_1=-.5, alpha_2=.2, alpha_3=.1, 
# beta_1=1, beta_2=.2, beta_3=.1
correct <- lm(y ~ trt*post*d+x*post*d, data = data)
misspec <- lm(y ~ trt*post*d+x*post, data = data)
gamma_1_mis <- misspec$coefficients["postTRUE:x"]
gamma_1 <- correct$coefficients["postTRUE:x"]
gamma_3 <- correct$coefficients["postTRUE:d:x"]

d_ATT_mis <- misspec$coefficients["trt:postTRUE"]+misspec$coefficients["trt:postTRUE:d"]
d_ATT_cor <- correct$coefficients["trt:postTRUE"]+correct$coefficients["trt:postTRUE:d"]
bias <- d_ATT_mis-d_ATT_cor
# misspec$coefficients["trt:postTRUE:d"]-correct$coefficients["trt:postTRUE:d"]

bias-b
bias-a

a <- -(gamma_1_mis-gamma_1-gamma_3)*(alpha_1+alpha_3)
b <- (beta_3) * (alpha_1+alpha_3)
b-a

# Multiple covariates ----------------------------------------------------

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





