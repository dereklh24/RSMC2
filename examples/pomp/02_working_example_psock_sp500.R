# *********************
# Set RSMC2 path
#*******************
#library("RSMC2")
devtools::load_all()
#**************************
#* Generated Data
#**************************

library(qrmdata)
data("SP500")

# Convert SP500 to flat tbl
library(dplyr)
library(xts)
library(zoo)

sp500_data <-
  tibble(date = as.Date(index(SP500)), close = as.double(SP500)) %>%
  mutate(lclose = log(close)) %>%
  mutate(ret = c(NA, diff(lclose))) %>%
  tail(2520) %>%
  mutate(t = 1:2520)

## sp500_data_in <- select(sp500_data)
end_T         <- nrow(sp500_data_in)

#**************************
#* Likelihood functions----
#**************************

library(pomp)

## Covariate table which has previous days' observation (for leverage effect)
ret_prev_covtbl <- covariate_table(ret_prev = c(42.0, head(sp500_data$ret, -1)), times = sp500_data$t)

pomp_base <- pomp(data = sp500_data, times = "t", t0 = sp500_data$t[1],
                  rinit    = Csnippet("h = rnorm(theta_h, pow(psi_h,2) + exp(lomega_h));"),
                  rmeasure = Csnippet("ret = mu + rnorm(0, exp(h * 0.5));"),
                  dmeasure = Csnippet("double l = dnorm(ret, mu, exp(h * 0.5), give_log);
                                        lik = (l > -1e30) ? l : -1e30;"),
                  rprocess = discrete_time(
                    step.fun = Csnippet("double sigma2_h = (ret_prev == 42.0) ? pow(psi_h, 2) + exp(lomega_h) : exp(lomega_h);
                                         double e = (ret_prev == 42.0) ? 0.0 : (ret_prev - mu) / exp(h * 0.5);
                                         h = h + kappa_h * (theta_h - h) + psi_h * e;
                                         h += rnorm(0.0, sqrt(sigma2_h));"),
                   delta.t = 1),
                  paramnames = c("mu", "theta_h", "kappa_h", "lomega_h", "psi_h"),
                  obsnames = c("ret"),
                  statenames = c("h"),
                  covar = ret_prev_covtbl
                  )

est_params <-
  c(mu = 0.000130108072007712, theta_h = -9.30001267751156, kappa_h = 0.0275880630220889, 
    lomega_h = -3.80862087929197, psi_h = -0.1997928622724)

system.time(
  pomp_node_test <- pomp_node$new(est_params, 10000, 2520, pomp_base)
)

system.time(
  pomp_node_test$run_pf(1:2520)
)
#**************************
#* Prior functions----
#**************************
prior_data <- list(
  mu      = list(mean = 0, var = 1e-8),
  theta_h = list(mean = 0, var = 100),
  kappa_h = list(mean = 1, var = 10),
  psi_h   = list(mean = 0, var = 10),
  lomega_h = list(alpha = 2.00005, beta = .0100005)
)

prior_sampler <- function(n_p, prior_data) {
  mu       <- rnorm(n_p, prior_data$mu$mean, sqrt(prior_data$mu$var))
  theta_h  <- rnorm(n_p, prior_data$theta_h$mean, sqrt(prior_data$theta_h$var))
  kappa_h  <- qnorm(
    runif(n_p, pnorm(0, prior_data$kappa_h$mean, sqrt(prior_data$kappa_h$var)),
          pnorm(2, prior_data$kappa_h$mean, sqrt(prior_data$kappa_h$var))),
    prior_data$kappa_h$mean, sqrt(prior_data$kappa_h$var))

  lomega_h <- log(invgamma::rinvgamma(n_p,
                                      prior_data$lomega_h$alpha,
                                      prior_data$lomega_h$beta))
  psi_h    <- rnorm(n_p, prior_data$psi_h$mean, sqrt(prior_data$psi_h$var))

  x        <- cbind(mu, theta_h, kappa_h, lomega_h, psi_h)
  return(x)
}

prior_sampler_dens <- function(params, prior_data) {
  dens <-
    dnorm(params[, 'mu'], prior_data$mu$mean,
          sqrt(prior_data$mu$var), log = TRUE) +
    dnorm(params[, 'theta_h'], prior_data$theta_h$mean,
          sqrt(prior_data$theta_h$var), log = TRUE) +
    dnorm(params[, 'kappa_h'], prior_data$kappa_h$mean,
          sqrt(prior_data$kappa_h$var), log = TRUE) +
    dnorm(params[, 'psi_h'], prior_data$kappa_h$mean,
          sqrt(prior_data$kappa_h$var), log = TRUE) +
    invgamma::dinvgamma(exp(params[, 'lomega_h']),
                        prior_data$lomega_h$alpha,
                        prior_data$lomega_h$beta, log=TRUE) +
    params[, 'lomega_h']

  return(dens)
}

prior_sampler <- make_closure(prior_sampler,
                              prior_data = prior_data)
prior_density <- make_closure(prior_sampler_dens,
                              prior_data = prior_data)

#**************************
#* PMCMC functions----
#**************************


parameter_bounds <-
  list(
    lower = list(
      mu = function(x) -Inf, theta_h = function(x) -Inf, lomega_h = function(x) -Inf,
      kappa_h = function(x) 0, psi_h = function(x) -Inf),
    upper = list(
      mu = function(x) Inf, theta_h = function(x) Inf, lomega_h = function(x) Inf,
      kappa_h =  function(x) 2, psi_h = function(x) Inf)
  )

pmcmc_proposal_f <- function(
  parameters,
  parameter_features,
  parameter_bounds,
  t,
  iter,
  #C      = (2.38^2) / 5
  C       = 1.13288
) {

  library(abind)

  n_parameters   <- nrow(parameters)
  V              <- parameter_features$cov * C

  lower          <- parameter_bounds$lower[colnames(parameters)]
  upper          <- parameter_bounds$upper[colnames(parameters)]

  new_parameters <-
    lapply(
      seq_len(n_parameters),
      function(idx) {
        r_nltmvtnorm_jump(
          parameters[idx,],
          mu =  parameters[idx,],
          sigma = V,
          lower = lower,
          upper = upper
        )
      }

    ) %>%
    abind::abind(along = 0)

  colnames(new_parameters) <- colnames(parameters)

  return(new_parameters)

}

pmcmc_proposal_dens_f <- function(
  old_parameters,
  new_parameters,
  parameter_bounds,
  parameter_features,
  t,
  iter,
  C       = 1.13288
) {
  library(abind)
  n_parameters   <- nrow(old_parameters)
  V              <- parameter_features$cov * C

  lower          <- parameter_bounds$lower[colnames(old_parameters)]
  upper          <- parameter_bounds$upper[colnames(old_parameters)]
  densities <- lapply(
    seq_len(n_parameters),
    function(idx) {
      d_nltmvtnorm_jump(
        old_parameters[idx,  ],
        new_parameters[idx,  ],
        mu = old_parameters[idx,  ],
        sigma = V,
        lower = lower,
        upper = upper
      )

    }

  ) %>% abind::abind(along=1)
  return(densities)
}

pmcmc_proposal_cmp <-
  make_closure(pmcmc_proposal_f, parameter_bounds = parameter_bounds)

pmcmc_proposal_dens_cmp <-
  make_closure(pmcmc_proposal_dens_f, parameter_bounds = parameter_bounds)


n_parameters <- 19*3
n_particles  <- 500
# Run for ten years
batch_t <- 2520


set.seed(1231231)
priors        <- prior_sampler(n_parameters)
set.seed(NULL)


library(parallel)
cl <- makePSOCKcluster(15)
invisible(clusterEvalQ(cl, devtools::load_all()))

batch_results <- density_tempered_pf(
  cluster_object              = cl,
  parameters                  = priors,
  prior_function              = prior_density,
  ## particle_mutation_lik_function = transition_lik_cmp
  pfilter_spec                = pomp_base,
  parameter_proposal_function = pmcmc_proposal_cmp,
  parameter_proposal_density  = pmcmc_proposal_dens_cmp,
  n_particles                 = n_particles,
  end_T                       = batch_t,
  gamma_threshold             = .5,
  pmcmc_mh_steps              = 20,
  grid_steps                  = 1e-13,
  adaptive_grid               = .01,
  check_lik_iter = c(1,2,5)
)

