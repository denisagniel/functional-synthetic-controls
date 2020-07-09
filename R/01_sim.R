library(tidyverse)
library(here)
library(glue)
library(augsynth)
library(fdapace)
library(longsurr)
library(refund)
library(fscm)
library(clustermq)

simfn <- function(n, m, s, k, run = 0) {
  library(tidyverse)
  library(here)
  library(glue)
  library(augsynth)
  library(fdapace)
  library(longsurr)
  library(refund)
  library(fscm)
  library(clustermq)
  set.seed(run)
  ds <- sim_data(n = n, m = m, sigma2 = s, delta = 1, K = k)
  
  trt_ds <- ds %>%
    select(id, trt) %>%
    unique
  pre_ds <- ds %>%
    filter(tt < 10)
  post_ds <- ds %>%
    filter(tt == 10)
  
  lin_fit <- fsc(id, trt, tn, y, pre_ds, post_ds, linear = TRUE)
  fg_fit <- fsc(id, trt, tn, y, pre_ds, post_ds, linear = FALSE)
  sc_base <- augsynth(form = y ~ trt, unit = id, time = tn, data = ds, t_int = m,
                      progfunc = 'None',
                      scm = TRUE)
  asc_base <- augsynth(form = y ~ trt, unit = id, time = tn, data = ds, t_int = m,
                       progfunc = 'ridge',
                       scm = TRUE)
  
  tibble(fsc = lin_fit$sc_est,
         afscl = lin_fit$asc_est,
         afscg = fg_fit$asc_est,
         scm_est = summary(sc_base)$average_att[1] %>% unlist,
         ascm_est = summary(asc_base)$average_att[1] %>% unlist)
}


sim_params <- expand.grid(n = c(15, 30, 50, 100, 500),
                          m = c(11, 20, 30, 100),
                          s = c(0.1, 1, 2, 5),
                          k = c(2, 4, 8, 16),
                          run = 1:3)

# tst <- sim_params %>% sample_n(2)
# # 
# tst
# Q_rows(tst, simfn)

options(
  clustermq.defaults = list(ptn="short",
                            log_file="Rout/log%a.log",
                            time_amt = "12:00:00"
  )
)
sim_res <- Q_rows(sim_params, simfn, 
                  fail_on_error = FALSE,
                  n_jobs = 20)
saveRDS(sim_res, here('results/01_sim-results.rds'))
