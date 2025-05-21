
source("tradeoff_fn.R", local = knitr::knit_global())

# New graphs

vary_gam <- vary_param_plot(param_name="gamma", param_vals = seq(0,2,.5), NREP=50)
vary_gam[2]
vary_gam[3]

vary_param_plot(param_name="alpha_2", param_vals = seq(0,2,.5), NREP=100)
vary_param_plot(param_name="alpha_3", param_vals = seq(0,2,.5), NREP=100)
vary_param_plot(param_name="beta_3", param_vals = seq(0,1,.3), NREP=100)
vary_param_plot(param_name="beta_2", param_vals = seq(0,1,.3), NREP=50)


# One problem here is that we can't control d, X sample size. 
# Q: How does varying these parameters affect the sample size? 
# Generate two data sets: alpha_2=.1, 2

# Looking at how treated proportion changes with alpha_2
nreps=100
prop_trt1 <- rep(NA,nreps)
prop_trt2 <- rep(NA,nreps)

for(i in c(1:nreps)){
  data1 <- make_data(n=100, alpha_2=0.1)
  data2 <- make_data(n=100, alpha_2=10) 
  
  prop_trt1[i] <- data1 %>% 
    filter(tp==1) %>%
    summarize(prop_trt=mean(trt)) %>% pull(prop_trt)
  prop_trt2[i]=data2 %>% 
    filter(tp==1) %>%
    summarize(prop_trt=mean(trt)) %>% pull(prop_trt)
}

hist(prop_trt1)
hist(prop_trt2)




