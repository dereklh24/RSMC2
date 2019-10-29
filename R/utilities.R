# Various utility functions 

#' Find mean likelihood value from vector of log-likelihoods
#'
#' @param x Vector of log-likelihoods
#'
#' @return The mean of the exponents of \code{x}
#' @export
#'
#' @examples
#' exp_mean(c(log(2),log(4),log(6)))
#' ## [1] 4
exp_mean <- function(x, underflow = TRUE, return_log = FALSE, log_weights = NULL) {
	if(underflow) {
		max_x      <- max(x)
		if(!is.null(log_weights))  {
		  w <- calc_sampling_weights(log_weights)
		  exp_x_plus <- sum(w * exp(x - max_x)) / sum(w)
		} else {
		  exp_x_plus <- mean(exp(x - max_x))
		}
		if(return_log) return(log(exp_x_plus) + max_x)
		else return(exp_x_plus * exp(max_x))
	}
	else {
		if(return_log) return(log(mean(exp(x))))
		else return(mean(exp(x)))
	}
}

exp_var <- function(x, underflow = TRUE, return_log = FALSE) {
        if(underflow) {
                max_x      <- max(x)
                exp_x_plus <- var(exp(x - max_x))
                if(return_log) return(log(exp_x_plus) + max_x * 2)
                else return(exp_x_plus * exp(max_x * 2))
        }
        else {
                if(return_log) return(log(var(exp(x))))
                else return(var(exp(x)))
        }
}


#' Calculate ESS from vector of likelihood weights
#'
#' @param x Vector of weights
#'
#' @return The ESS of those weights
#' @export
#'
#' @examples
#' calculate_ess(c(.2,.001,.01,.5))
#' ## [1] 1.742569
calculate_ess <- function(x, log_x = TRUE, underflow = TRUE) {
        #Prevent underflow
        if(log_x) x_scaled <- exp(x - as.integer(underflow) * max(x)  )
        else x_scaled <- x / (1 + (as.integer(underflow) * max(x) - 1))

        ess <- sum(x_scaled)^2 / sum(x_scaled^2)

        if(is.nan(ess)) ess <- 0
        return(ess)
}

#' Calculate sampling weights from log-likelihoods, compensating for underflow
#'
#' @param x Vector of log-weights
#' @param underflow Default TRUE. Subtracts the maximum of x from x before exponentiating to avoid underflows
#'
#' @return A vector of sampling weights. Note that these are not necessarily normalized to sum to one, but still work in sample()
#' @export
#'
#' @examples
#' calc_sampling_weights(c(-3,-4,-5,-2.5))
#' # [1] 0.6065307 0.2231302 0.0820850 1.0000000
#'
calc_sampling_weights <- function(x, underflow = TRUE){
        return(exp(x - as.integer(underflow) * max(x)))
}

#' Make a closure that is entirely self-contained
#'
#' @param ...fun The function
#' @param ... Arguments to pass to \code{purrr:partial}
#'
#' @return
#' @export
#' @importFrom purrr partial
#' @importFrom compiler cmpfun
#' @details This is a simple wrapper function to the \code{partial} function from purrr. It forces evaluation of the function and stores it into an enclosing environment. This allows for preserving the function even when passed to the cluster. It also byte-compiles the function.
#' @examples
#'
make_closure <- function(...fun, ...) {
        ...fun_cmp <- compiler::cmpfun(...fun)
        compiler::cmpfun(purrr::partial(.f = ...fun_cmp, .lazy = FALSE, ...))
}


#' Calculate the weighted mean and covariance from a sample of parameters
#'
#' @param parameters A matrix of parameters, each row a sample
#' @param log_weights The weights of the sample
#' @param log_features Boolean indicating whether to calculate features for the log of parameters
#'
#' @return A list with the weighted mean, weighted covariance, parameter values, and log weights
#' @export
#'
#' @examples
#' 
#' 
parameter_features <- function(parameters, log_weights, log_features = F) {
        weights      <- exp(log_weights - max(log_weights))
        total_weight <- sum(weights)

        theta_mean    <- mapply(
                function(row_idx, weight) {
                        weight * parameters[row_idx, ]
                },
                1:nrow(parameters),
                weights,
                SIMPLIFY = F
        )
        theta_mean %<>%
                do.call(what=rbind) %>%
                apply(2, sum) %>%
                magrittr::divide_by(total_weight)

        theta_cov <-
                mapply(
                        function(row_idx, weight) {
                                weight * (parameters[row_idx, ] - theta_mean) %*% t(parameters[row_idx, ] - theta_mean)
                        },
                        row_idx = 1:nrow(parameters),
                        weight = weights,
                        SIMPLIFY = F
                )

        theta_cov <- Reduce(`+`, theta_cov) %>%
                magrittr::divide_by(total_weight)

        rownames(theta_cov) <- colnames(theta_cov)

        return_list <-
                list(
                        mean = theta_mean,
                        cov  = theta_cov,
                        parameters = parameters,
                        log_weights = log_weights
                )

        if(log_features) {
                return_list$log_features <- parameter_features(log(parameters), log_weights, log_features = F)
        }


        return(return_list)
}

#' Make a folder path
#'
#' @param path_name Path of the folder
#' @param suppress_warnings If the folder is already there, suppress the warning.
#'
#' @return True or False
#' @export
#'
#' @examples
make_folder_path <- function(path_name, suppress_warnings = F){
  #Purpose: Takes the path of a folder, checks if it exists. If it does not, creates all
  #necessary sub-folders.

  #Input: path_name for the folder.

  #Output: Returns T for successful creation of a folder or if folder already existed.
        #Errors out if the designated folder is a file or if a folder cannot be successfully
        #created.

  isDir <- file.info(path_name)$isdir
  isDir[is.na(isDir)] <- "NA"

  success <- F
  switch(isDir,
         "TRUE"  = {
                        if(!suppress_warnings) warning(paste(path_name, "is an already existing directory."))
                        success <- T
                 },
         "FALSE" = stop(paste(path_name, "already exists, but is not a directory.")),
         "NA"    = {success <- dir.create(path_name, recursive=T);
                    if(!success) stop(paste(path_name, "was not successfully created.",
                                       "Check permissions along the path and that the path is valid."))
                   }
  )

  return(success)
}


#' Calculate the quantiles of a distribution with weights 
#'
#' @param x Vector of values
#' @param w Vector of weights
#' @param q Vector of quantiles
#'
#' @return
#' @export
#'
#' @examples
weighted_quantile <- function(x, w, q = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)) {
  i          <- order(x)
  q          <- sort(q)
  out        <- weighted_quantile_cpp(x[i], w[i], q)
  names(out) <- q
  return(out)
}


# This takes an IXJXlength(theta) array and casts it to an N_theta x length(theta) matrix
parameter_array_to_matrix <- function(params) {
  array(params, dim = c(dim(params)[1] * dim(params)[2], dim(params)[3]), dimnames = list(NULL, dimnames(params)[[3]]))
}

parameter_matrix_to_array <- function(params, I, J) {
  array(params, dim = c(I, J, ncol(params)), dimnames = list(NULL, NULL, dimnames(params)[[2]]))
}

# Calculates a biased estimate of the weighted covariance (but bias goes to zero as obs increases)
#' Title
#'
#' @param x 
#' @param l_w 
#'
#' @return
#' @export
#'
#' @examples
weighted_cov <- function(x, l_w) {
  ws   <- (exp(l_w - max(l_w)))

  if(sum(ws) == 0) {
    return(0)
  }

 # lm1  <- log(sum(x * ws)) - log(sum(ws))

  m1  <- sum(x * ws) / sum(ws)
  m2  <- sum(x^2 * ws) / sum(ws)
  
  m2 - m1^2
  
  #list(m1 = m1, m2= m2, v = m2 - m1^2)
}

