remotes::install_github('denisagniel/fscm')
# remotes::install_github('denisagniel/tPACE')
library(tidyverse)
library(here)
library(glue)
library(augsynth)
library(fdapace)
library(longsurr)
library(refund)
library(fscm)
library(clustermq)
library(synthdid)

simfn <- function(n, m, s, k, g, run = 0) {
  library(tidyverse)
  library(here)
  library(glue)
  library(augsynth)
  library(fdapace)
  library(longsurr)
  library(refund)
  library(fscm)
  library(clustermq)
  library(synthdid)
  
  set.seed(run)
  ds <- sim_data(n = n, m = m, sigma2 = s, delta = 1, K = k, gamma = g)
  
  yy <- ds %>%
    select(id, tn, y) %>%
    spread(tn, y) %>%
    arrange(-id) %>%
    select(-id) %>%
    as.matrix
  trt_ds <- ds %>%
    select(id, trt) %>%
    unique
  pre_ds <- ds %>%
    filter(tt < 10)
  post_ds <- ds %>%
    filter(tt == 10)
  
  sdid <- synthdid_estimate(yy, n-1, m-1)
  sdid_w <- attr(sdid, 'weights')
  sdid_lw <- sdid_w$lambda
  
  lin_fit <- fsc(id, trt, tn, y, pre_ds, post_ds, linear = TRUE)
  lw_fit <- fsc(id, trt, tn, y, pre_ds, post_ds, linear = TRUE, wts = sdid_lw)
  fg_fit <- fsc(id, trt, tn, y, pre_ds, post_ds, linear = FALSE)
  fgw_fit <- fsc(id, trt, tn, y, pre_ds, post_ds, linear = FALSE, wts = sdid_lw)
  liny_fit <- fsc(id, trt, tn, y, pre_ds, post_ds, linear = TRUE, include_y = TRUE)
  sc_base <- augsynth(form = y ~ trt, unit = id, time = tn, data = ds, t_int = m,
                      progfunc = 'None',
                      scm = TRUE)
  asc_base <- augsynth(form = y ~ trt, unit = id, time = tn, data = ds, t_int = m,
                       progfunc = 'ridge',
                       scm = TRUE)
  
  tibble(n, m, s, k, g, delta = 1,
                       fsc = lin_fit$sc_est,
                       fscw = lw_fit$sc_est,
                       afscl = lin_fit$asc_est,
                       afsclw = lw_fit$asc_est,
                       afscg = fg_fit$asc_est,
                       afscgw = fgw_fit$asc_est,
                       afscy = liny_fit$asc_est,
                       scm = summary(sc_base)$average_att[1] %>% unlist,
                       ascm = summary(asc_base)$average_att[1] %>% unlist,
                       scm2 = sc_estimate(yy, n-1, m-1),
                       did = did_estimate(yy, n-1, m-1),
                       sdid = synthdid_estimate(yy, n-1, m-1),
  )
}


sim_params <- expand.grid(n = c(15, 30, 100),
                          m = c(11, 25, 50, 100),
                          s = c(0.01, 0.1, 0.2),
                          k = c(2, 8, 16),
                          g = c(1.1, 2),
                          run = 1:1000)

# tst <- sim_params %>% sample_n(2)
# 
# tst
# Q_rows(tst, simfn, n_jobs = 1)
# tst <- sim_params %>%
#   sample_n(1)
# tst
# with(tst, simfn(n = n, m = m, s = s, k = k, g= g, run = run)) %>% data.frame
options(
  clustermq.defaults = list(ptn="medium",
                            log_file="Rout/log%a.log",
                            time_amt = "48:00:00"
  )
)
sim_res <- Q_rows(sim_params, simfn, 
                  fail_on_error = FALSE,
                  n_jobs = 250)
saveRDS(sim_res, here('results/01_sim-results.rds'))
