#' Title
#'
#' @param cluster_object
#' @param particle_mutation_function
#' @param likelihood_function
#' @param experiment_data
#' @param model_omega
#'
#' @return
#' @export
#' @include particle_node.R
#' @include utilities.R
#' @importFrom compiler cmpfun
#'
#' @examples
initialize_cluster_environment_psock <- function(
  cluster_object,
  ## particle_mutation_lik_function,
  pfilter_spec,
  forecast_mutation_lik_function = NULL
) {


  message("Exporting data to cluster....")
  export_time <- Sys.time()

  objects_to_export                 <- 'pfilter_spec'
  if("pomp" %in% class(pfilter_spec)) {
    ;
  } else {
    pfilter_spec  <- compiler::cmpfun(pfilter_spec)

    if(!is.null(forecast_mutation_lik_function)) {
      forecast_mutation_lik_function_c <- compiler::cmpfun(forecast_mutation_lik_function)
      objects_to_export                <- c(objects_to_export, 'forecast_mutation_lik_function_c')
    }
   }

  parallel::clusterExport(cluster_object, objects_to_export, envir = environment())
  message("Done exporting data.... (time = ", format(Sys.time() - export_time), ")")
  #This function is just to setup the environment
  return(NULL)
}

initialize_remote_particle_node_psock <- function(cluster_object, parameters, t_cycle = NULL, end_T, 
                                                n_particles, pn_list_name = "pn_list", save_history = FALSE) {
  # Create list of list of parameters
  parameter_list_nodes <- cluster_split(cluster_object, lapply(seq_len(nrow(parameters)), function(idx) parameters[idx, ]))

  ## browser()
  prior_likelihood_out <- parallel::parLapply(cluster_object, parameter_list_nodes, function(parameter_list, n_particles, t_cycle, end_T, save_history, pn_list_name) {
    # Setup particle_node R6 objects
    if("pomp" %in% class(pfilter_spec)) {
      assign(pn_list_name,
             purrr::map(parameter_list, pomp_node$new, n_particles, end_T, pfilter_spec,
                        save_history),
             envir = .GlobalEnv)
    } else {
      assign(pn_list_name,
            purrr::map(parameter_list, particle_node$new, n_particles, end_T, pfilter_spec,
                        save_history),
            envir = .GlobalEnv)
    }
    
    # Run them in parallel (if supplied with a range)
    if(!is.null(t_cycle)) {
      return(purrr::map_dbl(get(pn_list_name, envir = .GlobalEnv), ~.x$run_pf(t_cycle)))
    } else{
      return(NULL)
    }
  }, n_particles, t_cycle, end_T, save_history, pn_list_name)
  return(abind::abind(prior_likelihood_out, along = 1))
}

run_remote_particle_node_psock <- function(cluster_object, t_cycle, pn_list_name = "pn_list", extract_variables = NULL, extract_n = 10) {
  x <- clusterCall(cl = cluster_object, fun = function(t_cycle, extract_variables, extract_n, pn_list_name) {
    purrr::map(get(pn_list_name, envir = .GlobalEnv), function(pn) {
      lik      <- pn$run_pf(t_cycle)
      if(!is.null(extract_variables)) {
        x_sample <- pn$subsample_latent_var_current(extract_n, extract_variables)
      } else {
        x_sample <- NULL
      }
      return(list(lik = lik, x_sample = x_sample))
    })
  }, t_cycle = t_cycle, extract_variables = extract_variables, extract_n = extract_n, pn_list_name = pn_list_name)
  
  x <- purrr::flatten(x)
  lik      <- abind::abind(purrr::map(x, "lik"), along = 1)
  if(!is.null(extract_variables)) {
    x_sample <- abind::abind(purrr::map(x, "x_sample"), along = 1)
  } else {
    x_sample <- NULL
  }
  return(list(lik = lik, x_sample = x_sample))
}

run_forecasts_particle_node_psock <- function(cluster_object, fcast_start, fcast_windows, fcast_extract_n = NULL, pn_list_name = "pn_list") {
  x <- parallel::clusterCall(cl = cluster_object, fun = function(fcast_start, fcast_windows, fcast_extract_n, pn_list_name) {
    purrr::map(get(pn_list_name, envir = .GlobalEnv), function(pn) {
      pn$forecast_pf(fcast_start, fcast_windows, fcast_extract_n, forecast_mutation_lik_function_c)
    })
  }, fcast_start = fcast_start, fcast_windows = fcast_windows, fcast_extract_n = fcast_extract_n, pn_list_name = pn_list_name)

  x          <- purrr::flatten(x)
  fcast_dens <- abind::abind(purrr::map(x, "fcast_dens"), along = 1)

  if(!is.null(fcast_extract_n)) {
    fcast_sample <- abind::abind(purrr::map(x, "fcast_sample"), along = 1)
  } else {
    fcast_sample <- NULL
  }

  return(list(fcast_dens = fcast_dens, fcast_sample = fcast_sample))
}

### This will be developed into some sort of S3 Methods thing later
initialize_cluster_environment <- function(cluster_object, pfilter_spec, forecast_mutation_lik_function = NULL) {
  if ("SOCKcluster" %in% class(cluster_object)) {
    initialize_cluster_environment_psock(cluster_object, pfilter_spec, forecast_mutation_lik_function)
  } else if (cluster_object == "Rmpi") {
    initialize_cluster_environment_mpi(cluster_object, pfilter_spec, forecast_mutation_lik_function)
  } else {
    stop("Invalid cluster detected")
  }
}

initialize_remote_particle_node <- function(cluster_object, parameters, t_cycle = NULL, end_T, 
                                            n_particles, pn_list_name = "pn_list", save_history = FALSE) {
  if ("SOCKcluster" %in% class(cluster_object)) {
    initialize_remote_particle_node_psock(cluster_object = cluster_object, parameters, t_cycle, end_T, 
                                          n_particles, pn_list_name, save_history)
  } else if (cluster_object == "Rmpi") {
    initialize_remote_particle_node_mpi(cluster_object = cluster_object, parameters, t_cycle, end_T, 
                                        n_particles, pn_list_name, save_history)
  } else {
    stop("Invalid cluster detected")
  }
}

run_remote_particle_node <- function(cluster_object, t_cycle, pn_list_name = "pn_list", extract_variables = NULL, extract_n = 10) {
  if ("SOCKcluster" %in% class(cluster_object)) {
    run_remote_particle_node_psock(cluster_object, t_cycle, pn_list_name = pn_list_name, extract_variables = extract_variables, extract_n = extract_n)
  } else if (cluster_object == "Rmpi") {
    run_remote_particle_node_mpi(cluster_object, t_cycle, pn_list_name = "pn_list", extract_variables = extract_variables, extract_n = extract_n)
  } else {
    stop("Invalid cluster detected")
  }
}

run_forecasts_particle_node <- function(cluster_object, fcast_start, fcast_windows, fcast_extract_n = NULL, pn_list_name = "pn_list") {
  if ("SOCKcluster" %in% class(cluster_object)) {
    run_forecasts_particle_node_psock(cluster_object, fcast_start, fcast_windows, fcast_extract_n, pn_list_name)
  } else if (cluster_object == "Rmpi") {
    run_forecasts_particle_node_mpi(cluster_object, fcast_start, fcast_windows, fcast_extract_n, pn_list_name)
  } else {
    stop("Invalid cluster detected")
  }
}

cluster_split <- function(cl, seq) {
  if ("SOCKcluster" %in% class(cl)) {
    parallel::clusterSplit(cl, seq)
  } else if (cluster_object == "Rmpi") {
    cluster_split_mpi(cl, seq)
  } else {
    stop("Invalid cluster detected")
  }
}


