

make_data <- function(n=10, trt_prop=.5, d_prop_trt=.1, d_prop_ctrl=.2, interact_dtx=0){
  max.time <- 2
  trt.time <- 2
  
  dat <- expand.grid(id = 1:n, tp = 1:max.time) %>% arrange(id,tp) %>% group_by(id) %>%
    mutate(int=rnorm(1,0,sd=0.25), # random intercept
           p.trt=trt_prop, # probability of treatment
           trt=rbinom(1, 1, p.trt), # treatment
           p.d.trt=d_prop_trt, # probability in subgroup given treated
           p.d.ctrl=d_prop_ctrl, # probability in subgroup given not treated
           d=rbinom(1,1, p.d.trt*trt + p.d.ctrl*(1-trt)), # subgroup indicator
           x=rnorm(1, mean = 1 - 0.5*trt + .2*d +.1*trt*d, sd = 1),
           post=I(tp >= trt.time), # indicator of post-treatment period
           treated=I(post == 1 & trt == 1) # time-varying indicator if treated or not
    ) %>% 
    ungroup()
  
  dat <- dat %>% mutate(
    y0 = 1 + x*(1+.5*tp+.2*d+.1*d*tp*interact_dtx) + trt + int + ((tp - 2.5)^2)/10,
    y1 = 1 + x*(1+.5*tp+.2*d+.1*d*tp*interact_dtx) + trt + int + (1-.2*d) + ((tp - 2.5)^2)/10,
    y = 1 + x*(1+.5*tp+.2*d+.1*d*tp*interact_dtx) + trt + int + treated*(1-.2*d) + ((tp - 2.5)^2)/10) %>%
    group_by(id) %>% mutate(y.diff = y - lag(y)) %>% ungroup()
  
  return(dat)
}

newdat <- make_data(n=100,interact_dtx = 3)

# 2
twfe_2 <- lm(y ~ trt*post*d + x*post, newdat)
summary(twfe_2)$coefficients[c("trt:postTRUE","trt:postTRUE:d"),]

# 3
twfe_3 <- lm(y ~ trt*post*d + x*post + x*d, newdat)
summary(twfe_3)$coefficients[c("trt:postTRUE","trt:postTRUE:d"),]

# 4
twfe_4 <- lm(y ~ trt*post*d + x*post*d, newdat)
summary(twfe_4)$coefficients[c("trt:postTRUE","trt:postTRUE:d"),]


