mpi_remote_exec_error_checking <- function(...) {
  x <- mpi.remote.exec(..., simplify=FALSE)

  if("try-error" %in% purrr::map_chr(x, class)) {
    stop(paste0("Remote error (first encountered is displayed): ", purrr::detect(x, ~class(.x) == "try-error")))
  } else {
    return(x)
  }
}


particleSS_cl_obj <- 
  c(
    'particle_node',
    'exp_mean',              
    'calculate_ess',         
    'calc_sampling_weights', 
    'make_closure',          
    'parameter_features',    
    'make_folder_path',
    'particle_mutation_function_c',
    'likelihood_function_c',
    'particle_gibbs_function',
    'particle_node_forecast', 
    'forecast_mutation_function_c', 
    'forecast_sampling_function_c'
    
  )


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
initialize_cluster_environment_mpi <- function(
  cluster_object,
  particle_mutation_lik_function
) {
  if(!require(Rmpi)) {
    stop("Rmpi needs to be installed and an MPI cluster available to use.")
  }
  
  # The following objects are needed for particle filtering
  particle_filtering_objects <- c(
    'particle_mutation_lik_function_c'
  )
  
  particle_mutation_lik_function_c   <- compiler::cmpfun(particle_mutation_lik_function)
  
  message("Exporting data to cluster....")
  export_time <- Sys.time()
  
  for(object_str in particle_filtering_objects) {
    obj     <- list(objname=object_str,obj=get(object_str))
    mpi.bcast.cmd(cmd=.tmpRobj <- mpi.bcast.Robj(comm=1),
                  rank=0, comm=1)
    mpi.bcast.Robj(obj, rank=0, comm=1)
    mpi.bcast.cmd(cmd=assign(.tmpRobj$objname,.tmpRobj$obj), rank=0, comm=1)
  
  }
  

  message("Done exporting data.... (time = ", format(Sys.time() - export_time), ")")
  
  #This function is just to setup the environment
  return(NULL)
  
}


initialize_remote_particle_node_mpi <- function(cluster_object, parameters, t_cycle = NULL, end_T, 
                                            n_particles, pn_list_name = "pn_list", save_history = FALSE) {
  if(!require(Rmpi)) {
    stop("Rmpi needs to be installed and an MPI cluster available to use.")
  }
  # Create list of list of parameters
  parameter_list_nodes <- cluster_split_mpi(cluster_object, lapply(seq_len(nrow(parameters)), function(idx) parameters[idx, ]))

  # Send arguments across MPI Cluster
  mpi.bcast.Robj2slave(t_cycle)
  mpi.bcast.Robj2slave(n_particles)
  mpi.bcast.Robj2slave(end_T)
  mpi.bcast.Robj2slave(pn_list_name)
  mpi.bcast.Robj2slave(save_history)

  mpi.scatter.Robj2slave(parameter_list_nodes)
  # Initialize and run filter
  prior_likelihood_out <- mpi_remote_exec_error_checking({
    # Setup particle_node R6 objects
    assign(pn_list_name, purrr::map(
      parameter_list_nodes, 
      particle_node$new, 
      n_particles, 
      end_T, 
      particle_mutation_lik_function_c,
      save_history))
    # Run them in parallel (if supplied with a range)
    if(!is.null(t_cycle)) {
     purrr::map_dbl(get(pn_list_name), ~.x$run_pf(t_cycle))
    } else{
      NULL
    }
  })
  return(abind::abind(prior_likelihood_out, along = 1))
}

run_remote_particle_node_mpi <- function(cluster_object, t_cycle, pn_list_name = "pn_list", extract_variables = NULL, extract_n = 10) {
  if(!require(Rmpi)) {
    stop("Rmpi needs to be installed and an MPI cluster available to use.")
  }
  mpi.bcast.Robj2slave(t_cycle)
  mpi.bcast.Robj2slave(extract_variables)
  mpi.bcast.Robj2slave(extract_n)
  
  x <- mpi_remote_exec_error_checking({
    purrr::map(get(pn_list_name, envir = .GlobalEnv), function(pn) {
      lik      <- pn$run_pf(t_cycle)
      if(!is.null(extract_variables)) {
        x_sample <- pn$subsample_latent_var_current(extract_n, extract_variables)
      } else {
        x_sample <- NULL
      }
      list(lik = lik, x_sample = x_sample)
    })
  })
  x <- purrr::flatten(x)
  lik      <- abind::abind(purrr::map(x, "lik"), along = 1)
  if(!is.null(extract_variables)) {
    x_sample <- abind::abind(purrr::map(x, "x_sample"), along = 1)
  } else {
    x_sample <- NULL
  }
  return(list(lik = lik, x_sample = x_sample))
}


#' Determine how tasks will be divided on a cluster
#'
#' @param cl The cluster object
#' @param seq A list or array 
#'
#' @return A list as long as \code{cl} with each object in \code{seq} distributed as they would be on a call to cluster_lapply.
#' @export
#' @importFrom purrr map
#' @importFrom parallel splitIndices
#'
#' @details This distribution differs depending on whether you use a SNOW cluster or an RMPI cluster. Also, you cannot assume this distribution if using a load-balanced function
#' @examples
#' cluster_split(cl, 1:10)
#' [[1]]
#' [1] 1 2 3

#' [[2]]
#' [1] 4 5 6 7

#' [[3]]
#' [1]  8  9 10
#'
cluster_split_mpi <- function(cl, seq) {
  if(!require(Rmpi)) {
    stop("Rmpi needs to be installed and an MPI cluster available to use.")
  }
  
  cluster_size <- mpi.comm.size() - 1
  index_map    <- parallel::splitIndices(length(seq), cluster_size)
  purrr::map(index_map, ~ seq[.x])
  
}