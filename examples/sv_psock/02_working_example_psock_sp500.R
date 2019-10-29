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
  tail(2520)

sp500_data_in <- select(sp500_data, y1 = ret)
end_T         <- nrow(sp500_data_in)

#**************************
#* Likelihood functions----
#**************************

transition_lik_f <- function(parameters, particles, weights, t, y1,
                             util = list(exp_mean = exp_mean, calc_sampling_weights = calc_sampling_weights)){
  # Calculate our parameter values from the reparameterization ----
  sigma2_h             <- parameters['psi_h']^2 + exp(parameters['lomega_h'])
  rho_h                <- parameters['psi_h'] / sqrt(sigma2_h)
  kappa_h              <- parameters['kappa_h']
  theta_h              <- parameters['theta_h']
  uncond_var           <- sigma2_h/ (2 * kappa_h - kappa_h^2)

  if (t == 1) {
    # Generate the array of particles
    particles <- array(
      data = NA_real_,
      dim  = c(nrow(particles), 1),
      dimnames = list(NULL, c("h"))
    )
    mean_1               <- theta_h
    sigma2_1             <- sigma2_h
  } else {
    # Reweight particles
    assertthat::assert_that(!all(is.infinite(weights)))
    prev_w      <- util$calc_sampling_weights(weights)
    a           <- sample.int(length(prev_w), replace=TRUE, prob=prev_w)
    particles[] <- particles[a, ]

    e_1       <- (y1[t-1] - parameters['mu']) / exp(particles[, 'h'] * .5)
    sigma2_1  <- exp(parameters['lomega_h'])
    mean_1    <- particles[, 'h'] + kappa_h * (theta_h - particles[, 'h']) + parameters['psi_h'] * e_1
  }
  # Sample particles from transition equation
  particles[, 'h'] <- mean_1 + rnorm(nrow(particles)) * sqrt(sigma2_1)
  w                <- dnorm(y1[t] - parameters["mu"], mean = 0, sd = exp(particles[, 'h'] * .5), log = TRUE)
  if(all(is.infinite(w))) {
    lik <- -Inf
  } else {
    lik <- util$exp_mean(w, return_log = TRUE)
}
  return(list(particles = particles, w = w, lik = lik))
}

transition_lik_cmp <- make_closure(transition_lik_f,
                               y1         = sp500_data_in$y1,
                               util       = list(exp_mean = make_closure(exp_mean),
                                                 calc_sampling_weights = make_closure(calc_sampling_weights)))
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


#********************
# * Forecasting
#********************

# Note: t is fixed at the last day of available y (so can use info up to and including t). Horizon represents how many days
# after t we are forecasting
forecast_f <- function(parameters, particles, weights, t, horizon, fcast_window, extract_n = NULL, y1,
                             util = list(exp_mean = exp_mean, calc_sampling_weights = calc_sampling_weights)){
  # Calculate our parameter values from the reparameterization ----
  sigma2_h             <- parameters['psi_h']^2 + exp(parameters['lomega_h'])
  rho_h                <- parameters['psi_h'] / sqrt(sigma2_h)
  kappa_h              <- parameters['kappa_h']
  theta_h              <- parameters['theta_h']
  uncond_var           <- sigma2_h/ (2 * kappa_h - kappa_h^2)

  if (t == 1) {
    # Generate the array of particles
    particles <- array(
      data = NA_real_,
      dim  = c(nrow(particles), 1),
      dimnames = list(NULL, c("h"))
    )
    mean_1               <- theta_h
    sigma2_1             <- sigma2_h
  } else {
    if(!is.null(weights)) {
      assertthat::assert_that(!all(is.infinite(weights)))
      prev_w      <- util$calc_sampling_weights(weights)
      a           <- sample.int(length(prev_w), replace=TRUE, prob=prev_w)
      particles[] <- particles[a, ]
    }

    # If it is the first time we are forecasting, create a column for the leverage effect
    if(horizon == 1) {
      if(!('e_1' %in% colnames(particles))) {
        particles <- cbind(particles, e_1 = (y1[t] - parameters["mu"]) / (exp(particles[, 'h'] * 0.5)))
      } else {
        particles[, 'e_1'] <- (y1[t] - parameters["mu"]) / (exp(particles[, 'h'] * 0.5))
      }
    }

    sigma2_1  <- exp(parameters['lomega_h'])
    mean_1    <- particles[, 'h'] + kappa_h * (theta_h - particles[, 'h']) + parameters['psi_h'] * particles[, 'e_1']
  }

  particles[, 'h']   <- mean_1 + rnorm(nrow(particles), mean = 0, sd = sqrt(sigma2_1))
  particles[, 'e_1'] <- rnorm(nrow(particles))

  # If our current horizon isn't in our window, we can just skip this part to save time
  if(horizon %in% fcast_window) {
    # Forecasts
    current_sd <- exp(particles[, 'h'] * 0.5)

    # Density
    dens <-
      util$exp_mean(dnorm(y1[t+horizon] - parameters['mu'], 0, current_sd, log = TRUE),
                    return_log = TRUE)

    # Samples
    if(!is.null(extract_n)) {
      fcast_sample <- parameters['mu'] + sample(particles[, 'e_1'] * current_sd, extract_n, replace=TRUE)
    } else {
      fcast_sample <- NULL
    }
  } else {
    dens <- NULL
    fcast_sample <- NULL
  }

  return(list(particles = particles, fcast_dens = dens, fcast_sample = fcast_sample))
}

forecast_f_cmp <- make_closure(forecast_f,
                                   y1         = sp500_data_in$y1,
                                   util       = list(exp_mean = make_closure(exp_mean),
                                                     calc_sampling_weights = make_closure(calc_sampling_weights)))




n_parameters <- 19*3
n_particles  <- 500
# Run for one year
batch_t <- 252


set.seed(1231231)
priors        <- prior_sampler(n_parameters)
set.seed(NULL)


library(parallel)
cl <- makePSOCKcluster(15)
#invisible(clusterEvalQ(cl, library("RSMC2")))
invisible(clusterEvalQ(cl, devtools::load_all()))

batch_results <- density_tempered_pf(
  cluster_object              = cl,
  parameters                  = priors,
  prior_function              = prior_density,
  particle_mutation_lik_function = transition_lik_cmp,
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


# ***************************
# Online filtering
# ***************************
online_results <- smc2_particle_filter(
  cl,
  batch_results$parameters,
  n_particles = 1000,
  start_T = batch_t + 1,
  end_T = end_T,
  particle_mutation_lik_function = transition_lik_cmp,
  prior_function                 = prior_density,
  parameter_proposal_function    = pmcmc_proposal_cmp,
  parameter_proposal_density     = pmcmc_proposal_dens_cmp,
  extract_variables              = "h", save_steps = 100,
  forecast_mutation_lik_function = forecast_f_cmp,
  forecast_settings              =  list(
        window     = c(1,2,5,10,22),
        moments    = c(1,2,3,4),
        quantiles  = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
        forecast_n = 10
        ),
  save_prefix = "temp_online_")

