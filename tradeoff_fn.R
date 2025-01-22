library(xtable)
library(ggplot2)
library(broom)
library(modelr)


# Data generation -- single covariate -------------------------------------

make_data <- function(n=100, trt_prop=.5, 
                      d_prop_trt=.1, d_prop_ctrl=.2, 
                      d_x=1,dg_x=0,
                      noise=.2, trt_effect_noise=0, 
                      xt_y=1, dtx_y=0, dx_y=0){
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
           x=rnorm(1, mean = 1 - 0.5*trt + .2*d*d_x+.1*d*trt*dg_x, sd = 1),
           post=I(tp >= trt.time), # indicator of post-treatment period
           treated=I(post == 1 & trt == 1) # time-varying indicator if treated or not
    ) %>% 
    ungroup()
  
  dat <- dat %>% mutate(
    d_trt_effect = rnorm(1,-.2,sd=trt_effect_noise),
    y0 = 1 + x*(xt_y * tp+.2*d*dx_y+.1*d*tp*dtx_y) + trt + int + ((tp - 2.5)^2)/10,
    y1 = y0 + (1+d_trt_effect*d),
    y = treated*y1+(1-treated)*y0  ) %>%
    group_by(id) %>% mutate(y.diff = y - lag(y)) %>% ungroup()
  
  return(dat)
}

# make_data with df of parameters as input
make_data_df <- function(params){
  return( do.call(make_data, as.list(params) ) )}

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
    filter(term == "trt:postTRUE:d")%>% 
    select(-c(statistic, p.value, data, model)) %>%
    mutate(bias = estimate-true_satt)
  
  return(list('data'=sim_data,'ests'=ests))
}
  

# Data generation -- multiple covariates ----------------------------------

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


# Analyzing data -- multiple covariates -----------------------------------


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

