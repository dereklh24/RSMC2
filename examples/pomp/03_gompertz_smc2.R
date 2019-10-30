## ----prelims,echo=FALSE,cache=FALSE--------------------------------------
library(ggplot2)
library(knitr)
theme_set(theme_bw())
library(pomp)
stopifnot(packageVersion("pomp")>="2.1")
options(
  keep.source=TRUE,
  stringsAsFactors=FALSE,
  encoding="UTF-8",
  scipen=5
)
set.seed(1332379783L)


## ----gompertz-init,results="hide"----------------------------------------
library(pomp)
gompertz() -> gomp
theta <- coef(gomp)
theta.true <- theta

## ----gompertz-sim,include=FALSE------------------------------------------
gomp %>%
  window(start=1) %>%
  simulate(seed=340398091L) -> gomp


## ----gompertz-mif2-1,results='hide'--------------------------------------
## library(foreach)
## library(doParallel)
## registerDoParallel()

estpars <- c("r","sigma","tau")

## prior_data <- list(
##   r = list(lmean = 0, sigma = 1),
##   sigma = list(lme)
## )

prior_sampler <- function(n_p) {
  X <- rlnorm(n_p * 3, 0, 1) %>%
    matrix(nrow=n_p, dimnames = list(NULL, c("r", "sigma", "tau")))
  K <- 1.0
  X_0 <- 1.0

  cbind(K, X, X_0)
}

prior_sampler_dens <- function(params) {
  rowSums(dlnorm(params[, c("r", "sigma", "tau")], 0.0, 1.0, log=TRUE))
}

parameter_bounds <-
  list(
    lower = list(
      r      = function(x) 0.0,
      sigma  = function(x) 0.0,
      tau    = function(x) 0.0
    ),
    upper = list(
      r      = function(x) Inf,
      sigma  = function(x) Inf,
      tau    = function(x) Inf
    ))


pmcmc_proposal_f <- function(
  parameters,
  parameter_features,
  parameter_bounds,
  t,
  iter,
  C       = 1.888133
) {

  library(abind)

  n_parameters   <- nrow(parameters)
  V              <- parameter_features$cov[c("r", "sigma", "tau"), c("r", "sigma", "tau")] * C

  lower          <- parameter_bounds$lower[colnames(V)]
  upper          <- parameter_bounds$upper[colnames(V)]

  parameters_update <- parameters[, c("r", "sigma", "tau")]

  new_parameters <-
    lapply(
      seq_len(n_parameters),
      function(idx) {
        r_nltmvtnorm_jump(
          parameters_update[idx,],
          mu =  parameters_update[idx,],
          sigma = V,
          lower = lower,
          upper = upper
        )
      }

    ) %>%
    abind::abind(along = 0)

  colnames(new_parameters) <- c("r", "sigma", "tau")


  return(cbind(K = 1, new_parameters, X_0 = 1))

}

pmcmc_proposal_dens_f <- function(
  old_parameters,
  new_parameters,
  parameter_bounds,
  parameter_features,
  t,
  iter,
  C       = 1.888133
) {
  library(abind)
  n_parameters   <- nrow(old_parameters)
  V              <- parameter_features$cov[c("r", "sigma", "tau"), c("r", "sigma", "tau")] * C

  lower          <- parameter_bounds$lower[colnames(V)]
  upper          <- parameter_bounds$upper[colnames(V)]
  densities <- lapply(
    seq_len(n_parameters),
    function(idx) {
      d_nltmvtnorm_jump(
        old_parameters[idx, c("r", "sigma", "tau")],
        new_parameters[idx, c("r", "sigma", "tau")],
        mu = old_parameters[idx,  c("r", "sigma", "tau")],
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

n_parameters <- 15*3
n_particles  <- 500
batch_t <- 100


priors <- prior_sampler(n_parameters)

library(parallel)
cl <- makePSOCKcluster(15)
invisible(clusterEvalQ(cl, devtools::load_all()))

batch_results <- density_tempered_pf(
  cluster_object              = cl,
  parameters                  = priors,
  prior_function              = prior_sampler_dens,
  ## particle_mutation_lik_function = transition_lik_cmp
  pfilter_spec                = gomp,
  parameter_proposal_function = pmcmc_proposal_cmp,
  parameter_proposal_density  = pmcmc_proposal_dens_cmp,
  n_particles                 = n_particles,
  end_T                       = batch_t,
  gamma_threshold             = 0.5,
  pmcmc_mh_steps              = 20,
  grid_steps                  = 1e-13,
  adaptive_grid               = .01,
  check_lik_iter = c(1,2,5)
)
