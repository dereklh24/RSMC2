#' Determine whether or not to accept proposed parameter values
#'
#' @param prior_func A function that evaluates the prior log-density of a given parameter
#' @param current_parameters The current value of the parameters
#' @param new_parameters The new (proposed) values of the parameters
#' @param current_log_likelihood The current log-likelihoods of the parameters
#' @param new_log_likelihood The new log-likelihoods of the parameters
#' @param jumping_distribution A function that evaluates the log-density of the M-H proposal kernel
#' @param parameter_features_list Contains the weighted mean and covariance of the current parameters.
#' @param iter The iteration number (passed on to the proposal sampler and density).
#'
#' @importFrom assertthat assert_that noNA are_equal
#'
#' @include utilities.R
#'
#' @return A vector of boolean values indicating which new parameters to accept.
#'
#' @examples
mh_acceptance <- function(prior_func,
                          current_parameters,
                          new_parameters,
                          current_log_likelihood,
                          new_log_likelihood,
                          jumping_distribution,
                          parameter_features_list,
                          iter)
{
        #Sanity checks--------------------------------------
        assertthat::assert_that(assertthat::noNA(new_log_likelihood))
        assertthat::assert_that(assertthat::noNA(current_log_likelihood))

        #Calculate prior values
        current_prior                                               <- prior_func(current_parameters)
        new_prior                                                   <- prior_func(new_parameters)

        assertthat::assert_that(assertthat::noNA(current_prior))
        assertthat::assert_that(assertthat::noNA(new_prior))

        prior_ratio                                                 <- new_prior - current_prior
        assertthat::assert_that(assertthat::noNA(prior_ratio))

        #We need to assert that our current prior log density is not -Inf
        #We should never be at a parameter value with a prior value of 0
        assertthat::assert_that(sum(is.infinite(current_prior)) == 0)


        #Calculate proposal densities of the parameters
        jump_current_to_new <-
                jumping_distribution(
                        current_parameters,
                        new_parameters,
                        parameter_features = parameter_features_list,
                        t                  = t,
                        iter               = iter
                )
        jump_new_to_current <-
                jumping_distribution(
                        new_parameters,
                        current_parameters,
                        parameter_features = parameter_features_list,
                        t                  = t,
                        iter               = iter
                )

        assertthat::assert_that(assertthat::noNA(jump_current_to_new))
        assertthat::assert_that(assertthat::noNA(jump_new_to_current))

        jump_ratio                                                  <- jump_new_to_current - jump_current_to_new

        assertthat::assert_that(assertthat::noNA(jump_ratio))

        #Ratio of the likelihoods of the two parameters
        lik_ratio                                                   <- new_log_likelihood - current_log_likelihood

        #If both the new_log_likelihood and current_log_likelihood are -Inf, then
        #we reject the new proposal
        lik_ratio[new_log_likelihood == -Inf]                       <- -Inf
        assertthat::assert_that(assertthat::noNA(lik_ratio))

        #Calculate the probability of acceptance (ceiling at 1)
        acceptance_probs                                            <- prior_ratio + lik_ratio + jump_ratio

        #So we will get NA values if we have +/- Inf in the prior and -/+ Inf in the likelihood
        #To reconcile this, we will first defer to the prior. If there is a -Inf in the new_prior,
        #we will automatically reject it.
        acceptance_probs[new_prior == -Inf]                         <- -Inf

        #If we somehow have -Inf for our starting value in the prior, we should move
        acceptance_probs[current_prior == -Inf & new_prior != -Inf] <- Inf

        acceptance_probs[prior_ratio == -Inf]                       <- -Inf

        assertthat::assert_that(assertthat::noNA(acceptance_probs))

        #Simulate Bernoulli trials to determine acceptance
        acceptance                                                  <- log(runif(length(acceptance_probs))) < acceptance_probs


        #Test if there ara NA values in the acceptance vector
        assertthat::assert_that(assertthat::noNA(acceptance))
        return(acceptance)

}

#' A generic PMCMC algorithm
#'
#' @param cluster_object A cluster object compatible with SNOW
#' @param parameters A matrix of parameter samples, with each row a new sample.
#' @param t The current time-value
#' @param log_parameter_weights The sampling weights of the parameters (y_{s:t}), where s is the time since the last PCMC step
#' @param log_parameter_likelihood The log-likelihood of each parameter y_{1:t}
#' @param prior_function A function that evaluates the prior density of the parameter vector
#' @param parameter_proposal_function An M-H proposal kernel
#' @param parameter_proposal_density The density of the M-H proposal
#' @param n_particles Number of particles in the particle filters
#' @param pmcmc_mh_steps The number of times the PMCMC kernel should be run.
#' @param tempering_coefficient_resampling A tempering coefficient to be applied to the resampling weights
#' @param tempering_coefficient_mh A tempering coefficient to be applied to the M-H acceptance likelihoods
#'
#' @include utilities.R
#' @include particle_node.R
#' @importFrom assertthat assert_that noNA
#'
#' @return A list with three elements: ```parameters```, ```log_parameter_weights```, and ```acceptance_rates```.
#' @export
#'
#' @examples
pmcmc_step_parallel_standard <- function(
        cluster_object,
        parameters,
        t,
        log_parameter_weights,  #The conditional likelhood SINCE THE LAST RESAMPLNG
        log_parameter_likelihood, #The total likelihood of the parameter
        prior_function,
        parameter_proposal_function,
        parameter_proposal_density,
        n_particles,
        start_T = 1,
        end_T,
        pmcmc_mh_steps                    = 5,
        tempering_coefficient_resampling  = 1,
        tempering_coefficient_mh          = 1
) {


        n_parameters                   <- nrow(parameters)

        #Calculate parameter features before resampling for the M -H step

        parameter_features_list        <- parameter_features(parameters, tempering_coefficient_mh * log_parameter_weights, log_features = F)


        #Temper the weights (if needed)
        sampling_weights               <-
                calc_sampling_weights(tempering_coefficient_resampling * log_parameter_weights)

        parameter_ancestors            <- sample(
                seq_along(sampling_weights),
                length(sampling_weights),
                T,
                sampling_weights
                )

        parameters                     <- parameters[parameter_ancestors, ]
        log_parameter_likelihood       <- log_parameter_likelihood[parameter_ancestors]

        parameter_features_list$resampled_cov <- cov(parameters)
        
        # Code to generate nice output about current parameter sample (before rejuvenation)----
        if(require(skimr) & require(tidyr) & require(dplyr)) {
          options(tibble.width = Inf)
          parameter_stats <- 
            skimr::skim(as.data.frame(parameters)) %>%
            tbl_df
          
          parameter_stats$tidy_value <- format(parameter_stats$value, digits = 2)
          parameter_stats$tidy_value[parameter_stats$stat == "hist"] <- parameter_stats$level[parameter_stats$stat == "hist"]
          parameter_stats$level[parameter_stats$stat == "hist"] <- ""
          
          parameter_stats_clean <- parameter_stats %>%
            tidyr::unite(col = stat_level, stat, level, sep = "") %>%
            dplyr::select(-value) %>%
            tidyr::spread(key = stat_level, value = tidy_value) %>%
            select(var, mean = mean.all, sd = sd.all, min = min.all, q25 = `quantile25%`, med= median.all, q75 = `quantile75%`, max = max.all, hist)
          message("Parameter Summary (before rejuvenation):")
          print(parameter_stats_clean )
        } else {
          message("Parameter Mean:")
          print(colMeans(parameters))
          message("Covariance Matrix:")
          print(cov(parameters))
        }

        acceptance_rates   <- rep(0.0, pmcmc_mh_steps)
        total_acceptance   <- rep(FALSE, nrow(parameters))

        # Beginning of MH steps
        for (i in 1:pmcmc_mh_steps)
        {
                loop_start <- Sys.time()
                #Propose new parameters from the old ones using the proposal
                #function
                new_parameters         <- parameter_proposal_function(
                        parameters         = parameters,
                        parameter_features = parameter_features_list,
                        t                  = t,
                        iter               = i
                )

                #Now we run a particle filter using these parameters.
                filter_time            <- Sys.time()
                new_log_likelihood <-  initialize_remote_particle_node(cluster_object = cluster_object, new_parameters, 1:t, end_T, n_particles)
#                 new_parameter_list_nodes <- cluster_split(cluster_object, lapply(seq_len(nrow(new_parameters)), function(idx) new_parameters[idx, ]))
# 
#                 if (!("Rmpi" %in% cluster_object)) {
#                   new_log_likelihood_out   <- cluster_lapply(cl          = cluster_object,
#                                                              x           = new_parameter_list_nodes,
#                                                              t_cycle     = 1:t,
#                                                              n_particles = n_particles,
#                                                              fun         = function(parameter_list, t_cycle, n_particles) {
#                                                                for (idx in seq_along(pn_list)) {
#                                                                  pn_list[[idx]]$parameters <- parameter_list[[idx]]
#                                                                }
#                                                                
#                                                               # pn_list <<- purrr::map(parameter_list, ~particle_node$new(.x, n_particles, end_T, particle_mutation_function_c, likelihood_function_c))
#                                                                liks <- purrr::map_dbl(pn_list, ~.x$run_pf(t_cycle))
#                                                                return(list(liks = liks, mem_used = pryr::mem_used()))
#                                                              }) 
#                   
#                 } else {
#                   mpi.scatter.Robj2slave(new_parameter_list_nodes)
#                   t_cycle     = 1:t
#                   mpi.bcast.Robj2slave(t_cycle)
#                   new_log_likelihood_out <- mpi.remote.exec({
#                     for (idx in seq_along(pn_list)) {
#                       pn_list[[idx]]$parameters <- new_parameter_list_nodes[[idx]]
#                     }
#                     
#                     liks <- purrr::map_dbl(pn_list, ~.x$run_pf(t_cycle))
#                     
#                     list(liks = liks, mem_used = pryr::mem_used())
#                   }, simplify = FALSE)
# 
#                   if("try-error" %in% purrr::map_chr(new_log_likelihood_out, class)) {
#                     stop(paste0("Remote error (first encountered is displayed)\\n", purrr::detect(new_log_likelihood_out, ~class(.x) == "try-error")))
#                   }
#                 }
# #browser()
#                 new_log_likelihood <- abind::abind(purrr::map(new_log_likelihood_out, "liks"), along = 1)

                filter_time          <- Sys.time() - filter_time
                #Now we decide to accept or reject particles based on their
                #prior, likelihood, and proposal weights. By default, the
                #proposals are assumed to be symmetric
                acceptance          <-
                        mh_acceptance(
                                prior_func              = prior_function,
                                current_parameters      = parameters,
                                new_parameters          = new_parameters,
                                current_log_likelihood  = tempering_coefficient_mh * log_parameter_likelihood,
                                new_log_likelihood      = tempering_coefficient_mh * new_log_likelihood,
                                jumping_distribution    = parameter_proposal_density,
                                parameter_features_list = parameter_features_list,
                                iter                    = i
                        )

                assertthat::assert_that(assertthat::noNA(acceptance))

                total_acceptance <- total_acceptance | acceptance
                #Write out the new parameters
                parameters[acceptance]                  <- new_parameters[acceptance]
                #Write out new log weights
                log_parameter_likelihood[acceptance]    <- new_log_likelihood[acceptance]
              

                acceptance_rates[i] <- mean(acceptance)
                message(
                  "Iteration ", sprintf("%02d", i), ": Acceptance = ", 
                  format(acceptance_rates[i]), 
                  " (filter_time = ", format(filter_time, digits = 4))
                  #"; mem_used = ", format(mem_usage, digits = 4), " MB)")
                
               # message("Loop time = ", format(Sys.time() - loop_start))

        }
        message("Total acceptance = ", mean(total_acceptance))

        message("Running particle filter for all parameters...")

        
        log_parameter_likelihood <- initialize_remote_particle_node(cluster_object = cluster_object, parameters, 1:t, end_T, n_particles)
        
        # log_parameter_likelihood[!total_acceptance] <- no_acceptance_liks
        
        # Code to generate nice output about current parameter sample (after rejuvenation)----
        if(require(skimr) & require(tidyr) & require(dplyr)) {
          options(tibble.width = Inf)
          parameter_stats <- 
            skimr::skim(as.data.frame(parameters)) %>%
            tbl_df
          
          parameter_stats$tidy_value <- format(parameter_stats$value, digits = 2)
          parameter_stats$tidy_value[parameter_stats$stat == "hist"] <- parameter_stats$level[parameter_stats$stat == "hist"]
          parameter_stats$level[parameter_stats$stat == "hist"] <- ""
          
          parameter_stats_clean <- parameter_stats %>%
            tidyr::unite(col = stat_level, stat, level, sep = "") %>%
            dplyr::select(-value) %>%
            tidyr::spread(key = stat_level, value = tidy_value) %>%
            select(var, mean = mean.all, sd = sd.all, min = min.all, q25 = `quantile25%`, med= median.all, q75 = `quantile75%`, max = max.all, hist)
          message("Parameter Summary (after rejuvenation):")
          print(parameter_stats_clean )
        } else {
          message("Parameter Mean:")
          print(colMeans(parameters))
          message("Covariance Matrix:")
          print(cov(parameters))
        }
        
        return(list(
                parameters               = parameters,
                log_parameter_likelihood = log_parameter_likelihood,
                acceptance_rates         = acceptance_rates
        ) )
}

#' Perform a M-H jump based on a truncated-normal proposal distribution
#'
#' @param x Starting value
#' @param mu Mean of the proposal. Usually should be the starting value.
#' @param sigma Covariance of the proposal distribution
#' @param lower A list length M of functions that take as input an M-1 - length vector and return the corresponding lower bound
#' @param upper A list length M of functions that take as input an M-1 - length vector and return the corresponding upper bound
#'
#' @return A length(x) vector at the new position
#' @export
#'
r_nltmvtnorm_jump <- function(x, mu, sigma, lower, upper, fix = c()) {
  A_list      = lapply(
    seq_len(length(mu)), 
    function(j) {
      A                 = solve(sigma[-j, -j], sigma[-j, j])
      sigma_j           = sigma[j,j] - sigma[j, -j] %*% A
      return(list(A = A, sigma = sqrt(sigma_j)))
    }
  )
  
  for (j in seq_along(x)) {
    if(j %in% fix) {
      next
    }
    mean_j  = mu[j] + crossprod(A_list[[j]]$A, x[-j] - mu[-j])
    l       = (lower[[j]](x[-j]) - mean_j) / A_list[[j]]$sigma
    u       = (upper[[j]](x[-j]) - mean_j) / A_list[[j]]$sigma
    
    x[j]    = qnorm(runif(1, min = pnorm(l), max = pnorm(u))) * A_list[[j]]$sigma + mean_j 
  }
  
  return(x)
}

#' Calculate the log-density of a truncated multivariate jump
#'
#' @param x Initial value
#' @param x_2 New value
#' @param mu Mean of the proposal. Usually the initial value.
#' @param sigma Covariance of the proposal
#' @param lower A list length M of functions that take as input an M-1 - length vector and return the corresponding lower bound
#' @param upper A list length M of functions that take as input an M-1 - length vector and return the corresponding upper bound
#'
#' @return The log-density of the jump
#' @export
#'
d_nltmvtnorm_jump <- function(x, x_2, mu, sigma, lower, upper, fix = c()) {
  A_list      = lapply(
    seq_len(length(mu)), 
    function(j) {
      A                 = solve(sigma[-j, -j], sigma[-j, j])
      sigma_j           = sigma[j,j] - sigma[j, -j] %*% A
      return(list(A = A, sigma = sqrt(sigma_j)))
    }
  )
  
  lik <- 0
  
  x_current <- x
  
  for (j in seq_along(x)) {
    if(j %in% fix) {
      next
    }
    mean_j  = mu[j] + crossprod(A_list[[j]]$A, x_current[-j] - mu[-j])
    l       = (lower[[j]](x_current[-j])) 
    u       = (upper[[j]](x_current[-j])) 
    
    lik             = lik + 
      dnorm(x_2[j], mean = mean_j, sd = A_list[[j]]$sigma, log = T) - 
      log(pnorm(u) - pnorm(l))
    
    x_current[j]    = x_2[j]
  }
  
  return(lik)
}
