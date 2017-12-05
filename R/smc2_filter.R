#' Run a standard vanilla SMC^2 particle filter.
#'
#' @param cluster_object A cluster object either created by the snow/parallel package or the Rhpc package. Needs to have been initialized by initialize_cluster_environment().
#' @param parameters A matrix of parameters, where each row is a sample from the prior.
#' @param n_particles Number of particles to use in the filter
#' @param t_start Time period to start with in the data. If greater than 1, parameters are assumed to be drawn from \eqn{p(\theta | y_{1:tstart}}.
#' @param end_T Time period to end with in the data. This means the output posterior will end at end_T.
#' @param ess_threshold The threshold at which the parameter population is rejuvenated. Default is .5.
#' @param prior_function Function that evaluates the density of the priors
#' @param parameter_proposal_function M-H proposal for PMCMC step. \eqn{\theta_new \sim J(. | \theta_old)}
#' @param parameter_proposal_density Evaluates \eqn{J(\theta_new | \theta_old)}.
#' @param pmcmc_mh_steps Number of steps in the rejuvenation step. Defaults to 5.
#'
#' @details  This is a particle filter that maintains "state" on the class. We need to evaluate the likelihood at each T for SMC^2, so we need a way of maintaining particles in the cluster 
#' while we check the sample of particles for each likelihood. This is slightly more involved than batch sampling, but it allows for online filtering.
#' @return
#' @export
#' @include particle_node.R
#' @include utilities.R
#' @include pmcmc.R
#' 
#' @examples
smc2_particle_filter <- function(
  cluster_object,
  parameters,
  prior_function,
  particle_mutation_lik_function,
  parameter_proposal_function,
  parameter_proposal_density,
  n_particles,
  start_T,
  end_T,
  ess_threshold                  = 0.5,
  pmcmc_mh_steps                 = 10,
  prev_results                   = NULL,
  extract_variables              = NULL,
  extract_moments                = c(1,2,3,4),
  extract_quantiles              = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
  extract_n                      = 10,
  forecast_mutation_lik_function = NULL,
  forecast_settings              = list(
    window                          = c(1),
    moments                         = c(1,2,3,4),
    quantiles                       = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99),
    forecast_n                      = NULL
    ),
  save_steps                     = NULL,
  save_prefix                    = NULL) {

  if (start_T > 1) {
    T_range <- (1:(start_T - 1))
  } else {
    T_range <- NULL
  }
  
  # Send functions/data to cluster. If a function for forecasting is provided, it is exported to the cluster as well
  initialize_cluster_environment(cluster_object, particle_mutation_lik_function, forecast_mutation_lik_function)
  
  # Setup the initial likelihood of each parameter value
  lik <- initialize_remote_particle_node(cluster_object = cluster_object, parameters, T_range, end_T, n_particles)
  if(is.null(T_range)) {
    lik <- rep(0.0, nrow(parameters))
  }
  # Keep track of the parameter history in a list
  parameter_history <- purrr::set_names(list(parameters), (start_T - 1))
  
  # Incremental parameter weights
  w <- rep(0.0, nrow(parameters))
  
  
  # Keep track of ess history
  ess <- rep(NA_real_, end_T)
  
  # Placeholder that holds extracted latent variable information
  extracted_moments_df   <- NULL
  extracted_quantiles_df <- NULL
  
  # Placeholders that hold forecasting information
  fcast_dens_df <- dplyr::tbl_df(NULL)
  fcast_moments_df <- dplyr::tbl_df(NULL)
  fcast_quantiles_df <- dplyr::tbl_df(NULL)

  run_forecasts <- !is.null(forecast_mutation_lik_function)

  # Now run a for-loop between start_T and end_T ====
  for (t in start_T:end_T) {
    pf_output <-  run_remote_particle_node(cluster_object = cluster_object, t_cycle = t, extract_variables = extract_variables,
                                             extract_n = extract_n)
    
    log_weights <- pf_output$lik
    w           <- w + log_weights
    lik         <- lik + log_weights

    # If we extracted any particles, weight them by parameters and condense them
    if(!is.null(extract_variables)) {
      suppressPackageStartupMessages({
        library(dplyr)
        library(data.table)
        library(dtplyr)})

      dens_w <- exp(w - max(w))

      extracted_vars <- pf_output$x_sample %>%
        as.data.frame %>%
        tbl_df %>%
        mutate(dens_w = rep(dens_w, each = extract_n)) %>%
        tidyr::gather(-dens_w, key = "extracted_var", value = "value")

      # Calculate moments
      
      extracted_moments_df_t <- purrr::map_df(extract_moments, 
                                         ~group_by(extracted_vars, extracted_var) %>% 
                                           summarize(value = sum(dens_w * value^.x) / sum(dens_w)) %>% 
                                           mutate(moment = .x, t = t))
      extracted_moments_df      <- dplyr::bind_rows(extracted_moments_df, extracted_moments_df_t)
      
      
      extracted_quantiles_df_t <- 
        plyr::ddply(extracted_vars, "extracted_var", 
                    function(x, q) weighted_quantile(x$value, x$dens_w), 
                    q = extract_quantiles) %>%
        dplyr::mutate(t = t) %>%
        dplyr::select(extracted_var, t, everything())
      
      extracted_quantiles_df    <- bind_rows(extracted_quantiles_df, extracted_quantiles_df_t)
    }

    ess[t]    <- calculate_ess(w, log_x = T, underflow = T)
    message_t <- paste0("t = ", t, " (ESS = ", format(ess[t]), ")")
    message(message_t)

    # Rejuvenation stage ======
    if(ess[t] < ess_threshold * length(lik)) {
            message("Rejuvenating at t = ", t, " (ESS = ", format(ess[t]), ")")

            pmcmc_out <- pmcmc_step_parallel_standard(cluster_object, 
                                                      parameters                       = parameters, 
                                                      t                                = t,
                                                      log_parameter_weights            = w,
                                                      log_parameter_likelihood         = lik, 
                                                      prior_function                   = prior_function,
                                                      parameter_proposal_function      = parameter_proposal_function,
                                                      parameter_proposal_density       = parameter_proposal_density,
                                                      n_particles                      = n_particles,
                                                      end_T                            = end_T,
                                                      pmcmc_mh_steps                   = pmcmc_mh_steps,
                                                      # No tempering via coefficient, so set these to 1.
                                                      tempering_coefficient_resampling = 1, 
                                                      tempering_coefficient_mh         = 1
            )
            parameters          <- pmcmc_out$parameters
            parameter_history   <- c(parameter_history, purrr::set_names(list(parameters), t))
            last_acceptance     <- median(pmcmc_out$acceptance_rates)
            
      
            lik          <- pmcmc_out$log_parameter_likelihood
            w            <- rep(0.0, nrow(parameters))
            t_last_pmcmc <- t
    }

    # Run forecasts ==========================
    if(run_forecasts) {
      forecast_out <- run_forecasts_particle_node(
        cluster_object, 
        fcast_start = t, 
        fcast_windows = forecast_settings$window, 
        fcast_extract_n = forecast_settings$forecast_n)

      # Average out forecasting density for each horizon
      dens_w       <- exp(w - max(w))
      fcast_dens_t <- 
        purrr::array_branch(forecast_out$fcast_dens, 2) %>%
        purrr::map(function(l_dens) {
          ld_max <- max(l_dens)
          log(sum((exp(l_dens - ld_max) * dens_w))) - log(sum(dens_w)) + ld_max
        }) %>%
        as.data.frame %>%
        purrr::set_names(forecast_settings$window) %>%
        mutate(t = t)

      fcast_dens_df <- dplyr::bind_rows(fcast_dens_df, fcast_dens_t)

      # Calculate quantiles and moments of forecasts for each horizon (if fcast_n > 0)
      if(!is.null(forecast_settings$forecast_n)) {
        extracted_y          <- forecast_out$fcast_sample %>%
          as.data.frame %>%
          tbl_df %>%
          mutate(dens_w = rep(dens_w, each = extract_n)) %>%
          tidyr::gather(-dens_w, key = "horizon", value = "value") %>%
          mutate(horizon = factor(horizon, levels = as.character(forecast_settings$window), ordered = TRUE))

        fcast_moments_df_t   <- purrr::map_df(forecast_settings$moments, 
          ~group_by(extracted_y, horizon) %>%
          summarize(value = sum(dens_w * value^.x) / sum(dens_w)) %>% 
          mutate(moment = .x, t = t))

        fcast_moments_df     <- dplyr::bind_rows(fcast_moments_df, fcast_moments_df_t)
      
      
        fcast_quantiles_df_t <- 
          plyr::ddply(extracted_y, "horizon", 
                    function(x, q) weighted_quantile(x$value, x$dens_w), 
                    q = extract_quantiles) %>%
          dplyr::mutate(t = t) %>%
          dplyr::select(horizon, t, everything())
      
        fcast_quantiles_df   <- dplyr::bind_rows(fcast_quantiles_df, fcast_quantiles_df_t)
      }
    }
    
    # Save temporary output ========================
    if(!is.null(save_steps)) {
      if(t %% save_steps == 0) {
        temp_output <- list(
          parameters            = parameters,
          w                     = w,
          parameter_array       = parameter_history,
          log_likelihood        = as.double(lik),
          ess                   = ess
        )
        
        if(!is.null(extract_variables)) {
          temp_output$extracted_moments_df      <- extracted_moments_df  
          temp_output$extracted_quantiles_df    <- extracted_quantiles_df
        }

        if(!is.null(fcast_dens_df)) {
          temp_output$fcast_dens_df <- fcast_dens_df
        }

        if(!is.null(fcast_moments_df)) {
          temp_output$fcast_moments_df      <- fcast_moments_df  
          temp_output$fcast_quantiles_df    <- fcast_quantiles_df
        }

        message("Writing out at T = ", t)
        saveRDS(temp_output, paste0(save_prefix, "_", t, ".rds"), compress = FALSE)

      }
    }
  }

  output <- list(
    parameters            = parameters,
    parameter_array       = parameter_history,
    log_likelihood        = as.double(lik),
    ess                   = ess
  )
  
  if(!is.null(extract_variables)) {
    output$extracted_moments_df      <- extracted_moments_df  
    output$extracted_quantiles_df    <- extracted_quantiles_df
  }

  if(!is.null(fcast_dens_df)) {
    output$fcast_dens_df <- fcast_dens_df
  }

  if(!is.null(fcast_moments_df)) {
    output$fcast_moments_df      <- fcast_moments_df  
    output$fcast_quantiles_df    <- fcast_quantiles_df
  }
  
  return(output)
}
                           
