
#' Run a density-tempered Bayesian inference algorithm
#'
#' @param cluster_object A cluster object compatible with SNOW
#' @param parameters A matrix of parameter samples, with each row a new sample.
#' @param particle_mutation_lik_function 
#' @param prior_function A function that evaluates the prior density of the parameter vector
#' @param parameter_proposal_function For PMCMC; an M-H proposal kernel
#' @param parameter_proposal_density For PMCMC; the density of the M-H proposal
#' @param pmcmc_step_parallel A PMCMC function to be used for rejuventation
#' @param gamma_threshold For PMCMC; the ESS threshold below which we resample
#' @param pmcmc_mh_steps For PMCMC; the number of times the PMCMC kernel should be run.
#' @param n_particles The number of particles to use in the filter
#' @param end_T The number of iterations to go
#' @param grid_steps The size of the grid search when tempering the likelihood.
#' @param check_lik_iter A vector of integers at which to check the likelihood MC error.
#' @param acceptance_threshold If the mean acceptance drops below this threshold, we automatically check the likelihood and adjust the particles accordingly. Currently not implemented.
#' @param adaptive_grid If set above zero, the grid step size will double if the change in ESS while searching for Xi is below this threshold. Meant to speed up grid search
#' @param temp_output_path A string that denotes a file path (appended with 000_rds) to write each intermediate stage to. Useful if there is a crash or something
#' @param xi_start_value 
#' @param iter_start_value 
#'
#' @return An list of output from the particle filter.
#' @export
#' 
#' @include particle_node.R
#' @include pmcmc.R
#' @include utilities.R
#' 
#'
#'
density_tempered_pf <- function(
  cluster_object,
  parameters,
  particle_mutation_lik_function,
  prior_function,
  parameter_proposal_function,
  parameter_proposal_density,
  n_particles,
  end_T,
  gamma_threshold      = 0.5,
  pmcmc_mh_steps       = 5,
  grid_steps           = 1e-6,
  adaptive_grid        = 0,
  check_lik_iter       = FALSE,
  acceptance_threshold = 0.2,
  pmcmc_step_parallel  = pmcmc_step_parallel_standard,
  temp_output_path     = NULL,
  xi_start_value       = 0.0,
  iter_start_value     = 1
) {
        library(magrittr)
        n_parameters <- nrow(parameters)

        # Sanity Checks
        # This is a top-level algorithm, so we initialize the cluster
        initialize_cluster_environment(cluster_object, particle_mutation_lik_function)

        prior_likelihoods <- initialize_remote_particle_node(cluster_object = cluster_object, parameters, 1:end_T, end_T, n_particles)
        
        #Set up the loop
        iter                <- iter_start_value
        xi                  <- xi_start_value

        parameter_history   <- array(parameters, dim = c(dim(parameters), 1))
        likelihood_history  <- array(prior_likelihoods, dim = c(length(prior_likelihoods), 1))

        current_parameters  <- parameters
        sampling_likelihood <- prior_likelihoods

        ess_history         <- n_parameters

        log_lik_var_hist                <- double(length(check_lik_iter))
        names(log_lik_var_hist)         <- check_lik_iter

        #This is the main loop that exits only when the density-tempering power equals 1
        while(xi[iter] < 1) {
                iter     <- iter + 1
                message("Picking new value for Xi...")

                # Grid search algorithm ----
                # Here we move up in increments of step_size for our tempering exponent Xi
                # until the effective sample size (ESS) of our sampled parameters drops below a
                # threshold (i.e. (ESS)/(n_parameters) < gamma_threshold).
                #xi_diff    <- 0
                step_size  <- grid_steps
                xi_diff    <- step_size
                ess        <- calculate_ess(sampling_likelihood * xi_diff)
                ess_prev   <- ess
                
                while(ess / n_parameters > gamma_threshold & (xi_diff + xi[iter - 1]) < 1) {
                        xi_diff <-  xi_diff + step_size
                        ess     <-  calculate_ess(sampling_likelihood * xi_diff)
                        if(ess_prev - ess < adaptive_grid) {
                                step_size <- step_size * 2
                        }
                }
                xi[iter] <- min(xi_diff + xi[iter - 1], 1.0)

                message(paste0('Xi set to ', format(xi[iter]), ' (difference = ', format(xi_diff), ')'))
                message(paste0('ESS = ', format(ess)))

                # Rejuvenate the parameters with a PMCMC step----
                        pmcmc_out                <-
                        pmcmc_step_parallel(
                                cluster_object             = cluster_object,
                                parameters                 = current_parameters,
                                t                          = end_T,
                                end_T                      = end_T,
                                log_parameter_weights      = sampling_likelihood,
                                log_parameter_likelihood   = sampling_likelihood,
                                prior_function             = prior_function,
                                parameter_proposal_function = parameter_proposal_function,
                                parameter_proposal_density  = parameter_proposal_density,
                                n_particles                = n_particles,
                                pmcmc_mh_steps             = pmcmc_mh_steps,
                                tempering_coefficient_resampling = xi_diff,
                                tempering_coefficient_mh         = xi[iter]
                        )
                #Process the output from PMCMC. This includes updating the parameters,
                #log likelihood, and recording the history.

                current_parameters                 <- pmcmc_out$parameters
                sampling_likelihood                <- pmcmc_out$log_parameter_likelihood
                acceptance_rates                   <- pmcmc_out$acceptance_rates

                #History
                parameter_history                  <- abind(parameter_history, current_parameters)
                likelihood_history                 <- abind(likelihood_history, sampling_likelihood)
                ess_history                        <- c(ess_history, ess)

                # We also check likelihood at user-specified iterations ----
                # This is to diagnose our Monte Carlo error from the particle filter. If we have too
                # high of an error in our likelihood estimates, then we need to up the number of particles
                check_lik_var <- (iter %in% check_lik_iter)
                if(check_lik_var) {
                        message("Checking likelihood variance with estimated mean at iter=", iter, "...")

                        test_parameters <- matrix(
                            rep(colMeans(current_parameters), n_parameters),
                            nrow = n_parameters,
                            byrow = T)

                        #Reattach parameter names if they exist
                        colnames(test_parameters) <- colnames(parameters)
                        
                        test_log_likelihoods <- initialize_remote_particle_node(cluster_object = cluster_object, parameters, 1:end_T, end_T, n_particles, pn_list_name = "test_list")
                        log_lik_var <- var(test_log_likelihoods)

                        message("Variance in log-likelihood is: ",
                                format(log_lik_var),
                                " (tempered = ",
                                format(log_lik_var * xi[iter]^2),
                                ")")

                        log_lik_var_hist[as.character(iter)] <- log_lik_var

               }

                if(!is.null(temp_output_path)) {
                    saveRDS(
                        list(
                            parameters            = current_parameters,
                            parameter_array       = parameter_history,
                            ess                   = ess_history,
                            log_likelihood        = likelihood_history[, dim(likelihood_history)[2]],
                            xi                    = xi,
                            log_lik_var           = log_lik_var_hist
                        ),
                        paste0(temp_output_path, "_", sprintf("%03d", iter), ".rds")
                    )
                }


        }
        output <-
            list(
                        parameters            = current_parameters,
                        parameter_array       = parameter_history,
                        ess                   = ess_history,
                        log_likelihood        = likelihood_history[, dim(likelihood_history)[2]],
                        xi                    = xi,
                        log_lik_var           = log_lik_var_hist
                )

        return(output)

}
