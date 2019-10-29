# POMP node R6 class ----
library(R6)
#' Class that runs and maintains state of a particle filter (POMP)
#' @export
#' @import R6
#' @import data.table
#' @importFrom compiler cmpfun
#' @importFrom data.table data.table setkey setDT
#' @import pomp
#' @importFrom pomp pfilter
#' @include utilities.R
#' @details This is a class that facilitates SMC using a "pfilter" object from the POMP package
pomp_node <- R6Class(
  "pomp_node",
  inherit = particle_node,
  public  = list(
    # particle_node object members: ----
    pomp  = NULL,
    pfilterd_pomp = NULL,
    # Utility functions
    # particle_node methods: -----
    # r6 Constructor
    ## The "pomp oj"

    initialize = function(parameters, n_particles, end_T, pomp, save_history = FALSE){
      self$parameters                   <- parameters
      self$n_particles                  <- as.integer(n_particles)
      self$end_T                        <- as.integer(end_T)
      self$lik_t                        <- rep(NA_real_, end_T)
      self$particle_ess                 <- rep(NA_real_, end_T)
      self$current_particles            <- matrix(nrow = n_particles, ncol = 1)
      self$w                            <- rep(0.0, n_particles)
      self$save_history                 <- save_history
      self$pomp                         <- pomp

      self$pomp@params <- parameters
    },
    # Function that advances the particles forward from t-1 to t and returns the likelihood
    pf_one_step = function(t) {
      error("This feature not yet implemented for the pomp backend...")
      return(NULL)
    },
    # Main particle filter function. T_range is a range of t values (e.g. 1:100). 
    run_pf      = function(T_range) {
      assertthat::assert_that(T_range[1] == 1)
      self$pfilterd_pomp <- pomp::pfilter(self$pomp, Np = self$n_particles)
      return(self$pfilterd_pomp@loglik)
    },
    forecast_pf = function(fcast_start, fcast_windows, fcast_extract_n, forecast_mutation_lik_function_c) {
      error("Forecasting not yet implemented for pomp backend...")
  },
    # Useful for SMC^2 to get a filtered distribution of latent variables without returning all of them
    subsample_latent_var_hist = function(T_range, n_x, extract_var) {
      error("Subsampling not yet implemented for pomp backend")
    },
    subsample_latent_var_current = function(n_x, extract_var) {
      error("Subsampling not yet implemented for pomp backend")
    }
  ))
