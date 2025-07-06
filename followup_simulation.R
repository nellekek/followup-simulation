library(dplyr)
library(survival)
library(survminer)
library(tidyr)
library(purrr)
library(boot)

set.seed(123)

start_time <- Sys.time()    # To track elapsed time 

N <- 100      # Individuals per trial
K <- 1000     # Number of simulations per parameter combination
B <- 1000     # Number of bootstrap replicates

# Parameter grids
lambda_x_vals <- c(1, 1/2, 1/3, 1/5, 1/15)
lambda_c2_vals <- c(1/20, 1/30, 1/25, 1/15, 1/10)
ab_list <- list(c(2, 4), c(2, 5), c(2, 4), c(2, 4), c(2,5))

# Store all results
results_list <- list()

# Simulation loop
parameter_combo_index <- 1 

# Loop over parameter combinations 
for (i in 1:5) { 
  message(paste("Parameter combo:", i))
  
  lambda_x <- lambda_x_vals[i]
  lambda_c2 <- lambda_c2_vals[i]
  ab <- ab_list[[i]]
      
      a <- ab[1]
      b <- ab[2]
      
      # Calculate true follow-up time for this (a,b,lambda_c2)
      
      # ### True fup: median of C_obs (uncomment to use this method) ###
      # theta_exp <- log(2) / (lambda_c2)                              #
      #                                                                #
      # if (theta_exp < a) {                                           #
      #   true_fup <- theta_exp                                        #
      # } else {                                                       #
      #   # Define the function to solve numerically                   #
      #   g <- function(theta) {                                       #
      #     ((b - theta) / (b - a)) * exp(-(lambda_c2) * theta) - 0.5  #
      #   }                                                            #
      #   # Solve in [a, b]                                            #
      #   true_fup <- uniroot(g, lower = a, upper = b)$root            #
      # }                                                              #
      # ################################################################
      
      # ### True fup: median of T (uncomment to use this method) #################
      # theta_exp <- log(2) / (lambda_x + lambda_c2)                             #
      #                                                                          #
      # if (theta_exp < a) {                                                     #
      #   true_fup <- theta_exp                                                  #
      # } else {                                                                 #
      #   # Define the function to solve numerically                             #
      #   g <- function(theta) {                                                 #
      #     ((b - theta) / (b - a)) * exp(-(lambda_c2 + lambda_x) * theta) - 0.5 #
      #   }                                                                      #
      #   # Solve in [a, b]                                                      #
      #   true_fup <- uniroot(g, lower = a, upper = b)$root                      #
      # }                                                                        #
      # ##########################################################################

      # Store MSE, bias, coverage for each method
      method_mse <- c()
      method_bias <- c()
      method_coverage <- c()
      method_names <- c()
      
      # Storage for each simulation run
      method_estimates <- list(obs_time = c(),
                               cens_time = c(),
                               end_of_study = c(),
                               known_func = c(),
                               reverse_km = c())
      
      censoring_percents <- c()
      
      ci_coverage <- list(obs_time = c(),
                          cens_time = c(),
                          end_of_study = c(),
                          known_func = c(),
                          reverse_km = c())
      
      # Generate K datasets 
      for (k in 1:K) {
        message(paste("Dataset:", k))
        
        delta <- time <- X <- C1 <- C2 <- numeric(N)
        
        # Generate data for N individuals 
        for (n in 1:N) {
          X[n] <- rexp(1, lambda_x)
          C1[n] <- runif(1, a, b)
          C2[n] <- rexp(1, lambda_c2)
          
          if (X[n] == min(X[n], C1[n], C2[n])) {
            delta[n] <- 1
            time[n] <- X[n]
          } else {
            delta[n] <- 0
            time[n] <- min(C1[n], C2[n])
          }
        }
        
        censoring_percents[k]<- mean(delta == 0) * 100      # Check
        
        df <- data.frame(ID = 1:N, time, delta, X, C1, C2)
        
        # Method 1: Observation time
        method_estimates$obs_time[k] <- median(df$time)
        
        # Method 2: Censored time
        method_estimates$cens_time[k] <- median(df$time[df$delta == 0])
        
        # Method 3: End-of-study method
        method_estimates$end_of_study[k] <- median(C1)
        
        # Method 4: Known function time
        kft <- ifelse(df$delta == 0, df$time, C1)
        method_estimates$known_func[k] <- median(kft)
        
        # Method 5: Reverse KM
        survfit_km <- survfit(Surv(time = df$time, event = 1 - df$delta) ~ 1)
        med_km <- surv_median(survfit_km)$median
        method_estimates$reverse_km[k] <- ifelse(is.null(med_km), NA, med_km)
        
        # Bootstrap CI 
        # Loop over all methods
        for (method in names(method_estimates)) {
          # Define the function to compute the estimator on a dataset
          theta_hat_fun <- function(data, indices) {
            boot_data <- data[indices, ]
            if (method == "obs_time") {
              return(median(boot_data$time))
            } else if (method == "cens_time") {
              return(median(boot_data$time[boot_data$delta == 0]))
            } else if (method == "end_of_study") {
              return(median(boot_data$C1))
            } else if (method == "known_func") {
              kft <- ifelse(boot_data$delta == 0, boot_data$time, boot_data$C1)
              return(median(kft))
            } else if (method == "reverse_km") {
              survfit_km <- survfit(Surv(time = boot_data$time, event = 1 - boot_data$delta) ~ 1)
              return(surv_median(survfit_km)$median)
            }
          }
          
          # Apply bootstrapping
          boot_obj <- boot(data = df, statistic = theta_hat_fun, R = B)
          
          # Standard error and normal CI
          theta_hat <- method_estimates[[method]][k]
          se_hat <- sd(boot_obj$t, na.rm = TRUE)
          
          ci_low <- theta_hat - 1.96 * se_hat
          ci_high <- theta_hat + 1.96 * se_hat
          
          # Check if true value is within CI
          covered <- (true_fup >= ci_low) & (true_fup <= ci_high)
          ci_coverage[[method]] <- c(ci_coverage[[method]], as.integer(covered))
        }
        
        # Time estimation (print every 5 iterations)
        if (k %% 5 == 0 || k == 1) {
          elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
          avg_time <- elapsed / k
          remaining <- avg_time * (K - k)
          cat(sprintf("Param combo %d | Iteration %d/%d | Elapsed: %.1fs | Est. time left: %.1fs\n",
                      i, k, K, elapsed, remaining))
        }
      }
      
      # Calculate average censoring percentages
      avg_censoring_percent <- mean(censoring_percents)
      
      # Calculate bias, MSE, coverage per method
      for (method in names(method_estimates)) {
        estimates <- method_estimates[[method]]
        estimates <- estimates[!is.na(estimates)]
        mse <- mean((estimates - true_fup)^2)
        bias <- mean(estimates - true_fup)
        coverage <- 100 * mean(ci_coverage[[method]], na.rm = TRUE)
        
        method_mse <- c(method_mse, mse)
        method_bias <- c(method_bias, bias)
        method_coverage <- c(method_coverage, coverage)
        method_names <- c(method_names, method)
      }
      
      # Store all results
      results_list[[parameter_combo_index]] <- data.frame(
        lambda_x = lambda_x,
        lambda_c2 = lambda_c2,
        a = a,
        b = b,
        censoring_percent = avg_censoring_percent,
        method = method_names,
        bias = method_bias,
        mse = method_mse,
        coverage = method_coverage
      )
      parameter_combo_index <- parameter_combo_index + 1
}

# Combine results
simulation_results <- bind_rows(results_list)
