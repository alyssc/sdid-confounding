
max.time <- 2
trt.time <- 2

# number of units in d
d_n <- round( d_prop*n )

dat <- expand.grid(id = c(1:n), tp = 1:max.time) %>% arrange(id,tp) %>% group_by(id) %>%
  mutate(int=rnorm(1,0,sd=noise), # random intercept
         d=1*I(id <= d_n ), # subgroup indicator
         v0=rnorm(1, 1, 2), 
         s0=rnorm(1,3,2), 
         w0=rnorm(1,alpha_1*d,1), 
         trt=rbinom(1,1,exp(v0+ alpha_2*w0)/(1+exp(v0+alpha_2*w0))), 
         post=I(tp >= trt.time), # indicator of post-treatment period
         treated=I(post == 1 & trt == 1) # time-varying indicator if treated or not
  ) %>% 
  ungroup()

dat <- dat %>% mutate(
  d_trt_effect = rnorm(1,-.2,sd=trt_effect_noise),
  y0 = 1 + x*( beta_1*tp + beta_2*d + beta_3*d*tp) + s0 + trt + int + ((tp - 2.5)^2)/10, #untreated PO
  y1 = y0 + (1+d_trt_effect*d), # treated PO
  y = treated*y1+(1-treated)*y0  ) %>% # observed outcome at time tp
  group_by(id) %>% mutate(y.diff = y - lag(y)) %>% ungroup()

table(dat$trt)
table(dat$trt,dat$d)


dat <- make_data(alpha_0=-2.5, alpha_2=1.5)
table(dat$trt)
table(dat$trt,dat$d)



