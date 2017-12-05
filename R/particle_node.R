# particle node R6 class ----
library(R6)
#' Class that runs and maintains state of a particle filter
#' @export
#' @import R6
#' @import data.table
#' @importFrom compiler cmpfun
#' @importFrom data.table data.table setkey setDT
#' @include utilities.R
#' @details This is a class that facilitates SMC. transition_function:(x_{t-1}, w_{t-1}) => (x_t, w_t, lik_t)
particle_node <- R6Class(
  "particle_node",
  public  = list(
    # particle_node object members: ----
    parameters                   = NULL, 
    n_particles                  = NULL,
    lik                          = NULL,
    w                            = NULL,
    current_particles            = NULL,
    current_ancestors            = NULL,
    resampled_particles          = NULL,
    end_T                        = NULL,
    lik_t                        = NULL,
    particle_ess                 = NULL,
    particle_mutation_lik_function_c = NULL,
    save_history                = NULL,
    particle_history            = NULL,
    # Utility functions
    util = list(
      exp_mean              = exp_mean,
      calculate_ess         = calculate_ess,
      calc_sampling_weights = calc_sampling_weights,
      make_closure          = make_closure,
      parameter_features    = parameter_features,
      make_folder_path      = make_folder_path,
      weighted_quantile     = weighted_quantile
    ),
    # particle_node methods: -----
    # r6 Constructor
    initialize = function(parameters, n_particles, end_T, particle_mutation_lik_function_c, save_history = FALSE){
      self$parameters                   <- parameters
      self$n_particles                  <- as.integer(n_particles)
      self$end_T                        <- as.integer(end_T)
      self$lik_t                        <- rep(NA_real_, end_T)
      self$particle_ess                 <- rep(NA_real_, end_T)
      self$current_particles            <- matrix(nrow = n_particles, ncol = 1)
      self$w                            <- rep(0.0, n_particles)
      self$save_history                 <- save_history
      
      
      self$particle_mutation_lik_function_c <- particle_mutation_lik_function_c
      
    },
    
    # Function that advances the particles forward from t-1 to t and returns the likelihood
    pf_one_step = function(t) {
      temp_out <- self$particle_mutation_lik_function_c(
        parameters = self$parameters,
        particles  = self$current_particles,
        weights    = self$w,
        t          = t)
      
      purrr::walk2(c('current_particles', 'w', 'lik'), temp_out, ~assign(.x, .y, envir = self))
      
      if(is.null(self$particle_history) & self$save_history) {
        
        self$particle_history <- replicate(self$end_T, 
                                           matrix(NA_real_, 
                                                  nrow = nrow(self$current_particles), 
                                                  ncol = ncol(self$current_particles) + 1, 
                                                  dimnames = list(NULL, c(colnames(self$current_particles), 'w'))), simplify = FALSE)
        
      }
      
      
      # # Record value in list
      if(self$save_history) self$particle_history[[t]][] <- cbind(self$current_particles, self$w)
      
      # Sanity checks on output
      stopifnot(all(!is.na(self$w)))
      stopifnot(!is.na(self$lik))
      
      if(is.infinite(self$lik) | all(is.infinite(self$w))) return(-Inf)
      assertthat::assert_that(length(self$particle_ess) >= t)
      self$particle_ess[t]     <- self$util$calculate_ess(self$w, log_x = T, underflow = T)
      
      return(self$lik)
    },
    # Main particle filter function. T_range is a range of t values (e.g. 1:100). 
    run_pf      = function(T_range) {
      
      for (t in T_range) {
        assertthat::assert_that(length(self$lik_t) >= t)
        self$lik_t[t] <- self$pf_one_step(t)
        # If our marginal likelihood is -Inf, the whole thing is -Inf
        if(is.infinite(self$lik_t[t])) return(-Inf)
      }
      
      # Garbage collection
      invisible(gc(FALSE))
      #if(class(try(self$lik_t[T_range])) == "try-error") stop("sdfasdf")
      return(sum(self$lik_t[T_range]))
    }, 
    
    forecast_pf = function(fcast_start, fcast_windows, fcast_extract_n, forecast_mutation_lik_function_c) {
      current_particles <- self$current_particles
      weights           <- self$w
      fcast_range       <- 1:max(fcast_windows)
      

      fcast_dens_mat    <- matrix(NA_real_, nrow = 1, 
                                  ncol = length(fcast_windows), 
                                  dimnames = list(NULL, fcast_windows))
      
      if(!is.null(fcast_extract_n)) {
        fcast_sample_mat  <- matrix(NA_real_, nrow = fcast_extract_n,
                                    ncol =  length(fcast_windows),
                                    dimnames = list(NULL, fcast_windows))
      } else {
        fcast_sample_mat  <- NULL
      }

      
      for (horizon in fcast_range) {
        fcast_out <- forecast_mutation_lik_function_c(
          self$parameters, 
          current_particles,
          weights, 
          t = fcast_start,
          horizon = horizon, 
          fcast_window = fcast_windows,
          extract_n = fcast_extract_n)
      
      weights           <- NULL
      current_particles <- fcast_out$particles
      
      if(horizon %in% fcast_windows) {
        fcast_dens_mat[, as.character(horizon)]   <- fcast_out$fcast_dens
        if(!is.null(fcast_extract_n)) {
         # stop(length(fcast_out$fcast_sample))
          fcast_sample_mat[, as.character(horizon)] <- fcast_out$fcast_sample
        }
      }
      
      }
      
    return(list(fcast_dens = fcast_dens_mat, fcast_sample = fcast_sample_mat))
  },
      
      
    
    # Useful for SMC^2 to get a filtered distribution of latent variables without returning all of them
    subsample_latent_var_hist = function(T_range, n_x, extract_var) {
      if(is.null(self$particle_history)) stop("Error: Particle history not available")
      
      subsample <- purrr::map(self$particle_history[T_range], function(particles) {
        w <- self$util$calc_sampling_weights(particles[, 'w'])
        a <- sample.int(length(w), n_x, replace=TRUE, w)
        particles[a, extract_var, drop=FALSE]
      })
      
      purrr::set_names(subsample, T_range)
    }, 
    
    subsample_latent_var_current = function(n_x, extract_var) {
      w <- self$util$calc_sampling_weights(self$w)
      a <- sample.int(length(w), n_x, replace=TRUE, w)
      self$current_particles[a, extract_var, drop=FALSE]
    }
  ))