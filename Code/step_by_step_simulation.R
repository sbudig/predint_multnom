# Step-by-Step Simulation Script
# Run this line-by-line to inspect intermediate results
# All variables will be available in your R session for inspection

# =============================================================================
# SETUP
# =============================================================================

library(here)
library(dplyr)
library(VGAM)
library(MCMCpack)
library(mvtnorm)
library(cmdstanr)

setwd(here())

# Load source functions
source("Code/predint_mult_source.R")

# Compile Stan models
mod_gamma <- cmdstanr::cmdstan_model("Code/dirichlet_multinomial_gamma.stan")
mod_cauchy <- cmdstanr::cmdstan_model("Code/dirichlet_multinomial_cauchy.stan")

# =============================================================================
# DEFINE TEST SETTINGS
# =============================================================================

# Simulation parameters
df_sim_settings <- expand.grid(
  k = c(5, 20),
  n = c(10, 50),
  m = c(1),
  phi = c(1.01, 5),
  B = 10000,
  quant_min = 0.01,
  quant_max = 10,
  alpha = 0.05,
  dispersion = "afroz",
  alternative = "both",
  tol = 1e-3,
  n_bisec = 60
)

# Probability vector (true probabilities for each category)
m_true_pi_vec <- matrix(c(0.1, 0.3, 0.6), ncol = 3)

# Method flags (enable only what you want to test)
l_methods <- list(
  bisection = FALSE,
  bisection_asym = FALSE,
  Bonferroni = FALSE,
  Bonf_bisec = FALSE,
  pointwise = FALSE,
  percentile_bt1 = FALSE,
  mvn = FALSE,
  percentile_bt_gemini = FALSE,
  nonparametric_bt = FALSE,  
  bayesian_mcmc_gamma = TRUE, 
  bayesian_mcmc_cauchy = TRUE 
)

# store original parameter settings dataframe
df_sim_settings_og <- df_sim_settings

# future cluster size set to historical cluster size, when m=1
# that means equal cluster size for historical and future samples is used
df_sim_settings$m[df_sim_settings$m == 1] <- df_sim_settings$n[df_sim_settings$m == 1]

# generates a list for each setting, where each element consists of the 
# historical data
# each element of the list is: $\bm{x}_k=(x_k1,...,x_kC)$
l_x_k_vec <- lapply(apply(df_sim_settings, 1, function(x)
  apply(m_true_pi_vec, 1, function(y)
    f_rdirmultinom(
      phi = as.numeric(x[4]),
      k = as.numeric(x[1]),
      n = as.numeric(x[2]),
      v_pi = y
    ), simplify = FALSE)), function(z)
      cbind(as.data.frame(do.call(rbind, z)),
            "level" = as.factor(rep(
              seq(1:nrow(m_true_pi_vec)),
              each = (nrow(as.data.frame(do.call(
                rbind, z
              ))) / nrow(m_true_pi_vec))
            ))))

# generates a list for each setting, where each element consists of the 
# future validation sample
# $\bm{y}=(y_1,...,y_C)$
l_y_vec <- lapply(apply(df_sim_settings, 1, function(x)
  apply(m_true_pi_vec, 1, function(y)
    f_rdirmultinom(
      phi = as.numeric(x[4]),
      k = 1,
      n = as.numeric(x[3]),
      v_pi = y
    ), simplify = FALSE)), function(z)
      cbind(as.data.frame(do.call(rbind, z)),
            "level" = as.factor(rep(
              seq(1:nrow(m_true_pi_vec)),
              each = (nrow(as.data.frame(do.call(
                rbind, z
              ))) / nrow(m_true_pi_vec))
            ))))

# track if there is any column that consists of just zeroes
v_historical_zerocounts <-
  unlist(sapply(l_x_k_vec, f_countzeroes))

# Iterate through all lists and check whether any category has just zeroes
# If so: replace random row of that category with a 1
# this prevents extremely large standard errors and estimates
l_x_k_vec <- lapply(l_x_k_vec, function(df) {
  for (col in colnames(df)) {
    if (all(df[[col]] == 0)) {
      zero_indices <- which(df[[col]] == 0)
      random_row <- sample(zero_indices, 1)
      df[[col]][random_row] <- 1
    }
  }
  return(df)  # Return the modified dataframe
})

# take each element (historical data sets for each seeting) of l_x_k_vec 
# and fit a multinomial model the elememt
l_models_x <- lapply(l_x_k_vec, function(df) {
  f_fit_mult(df_i = df, modeltype = "multinomial")
})

# For each model element of the list l_models_x, 
# errors and warnings are extraced and stored in the a new list
l_errwarn <- lapply(l_models_x, `[[`,2)

# write the errors into a dataframe
df_err_warn <- f_check_err_warn(l_errwarn)

# extract modelfit list
l_models_x <- lapply(l_models_x, `[[`,1)

#The respective dispersion parameter is estimated for each model 
# in l_models_x.  NA is returned for missing models.
l_phi_hat <- mapply(function(model, disp) {
  f_estimate_phi(model = model, dispersion = as.character(disp))
}, l_models_x, df_sim_settings$dispersion, SIMPLIFY = FALSE)

# for each model of l_models_x the probabilities for each category
# are estimated. $\hat{\bm{\pi}}=(\hat{\pi}_1,...,\hat{\pi}_C)
l_pi_hat_vec <- lapply(l_models_x, function(model) {
  if (!is.null(model)) {
    f_calc_pi(model)
  } else {
    return(NULL)
  }
})

# for each model element of the list l_models_x the total number of observations
# across all historical studies is calculated.
# $N_{hist}=\sum_{k=1}^{K}n_k$
l_N_hist <- sapply(l_x_k_vec, function(df)
  sum(rowSums(df[grep("^V", names(df))])))

# Remove l_mult_dat to free up memory (not removed, need for later methods)
# rm(l_x_k_vec)
# gc()

# Expected count for each category
# $\hat{\bm{y}}=(\bm{y}_1,...,\bm{y}_C$
l_y_hat_vec <- mapply(
  FUN = f_calc_y_hat,
  pi_hat = l_pi_hat_vec,
  m = as.list(df_sim_settings$m),
  SIMPLIFY = FALSE
)

# Variance of future observation for each setting
l_var_y_vec <- mapply(
  FUN = f_calc_var_y,
  phi_hat = l_phi_hat,
  m = as.list(df_sim_settings$m),
  pi_hat = l_pi_hat_vec,
  SIMPLIFY = FALSE
)

# Variance of the estimated future observation for each setting
l_var_y_hat_vec <- mapply(
  FUN = f_calc_var_y_hat,
  phi_hat = l_phi_hat,
  m = as.list(df_sim_settings$m),
  pi_hat = l_pi_hat_vec,
  N_hist = l_N_hist,
  SIMPLIFY = FALSE
)

# Prediction standard error for each category for each setting
l_se_pred_vec <- mapply(
  FUN = f_calc_se_pred,
  var_y = l_var_y_vec,
  var_y_hat = l_var_y_hat_vec,
  SIMPLIFY = FALSE
)

# change all estimated \hat{\phi} that are <=1 to 1.01
l_phi_hat <- f_replace_phi_below_one(l_phi_hat)

# phi needs to be adjusted when close to cluster size n or m
l_phi_hat_m <- list()
l_phi_hat_n <- list()
l_c_phi_hat_m <- vector()
l_c_phi_hat_n <- vector()

for(i in seq_along(l_phi_hat)){
  if(df_sim_settings$m[i] < l_phi_hat[[i]][1]){
    l_phi_hat_m[[i]] <- df_sim_settings$m[i]*0.975
    l_c_phi_hat_m[[i]] <- 1
  } else {
    l_phi_hat_m[[i]] <- l_phi_hat[[i]]
    l_c_phi_hat_m[[i]] <- 0
  }
  if(df_sim_settings$n[i] < l_phi_hat[[i]][1]){
    l_phi_hat_n[[i]] <- df_sim_settings$n[i]*0.975
    l_c_phi_hat_n[[i]] <- 1
  } else {
    l_phi_hat_n[[i]] <- l_phi_hat[[i]]
    l_c_phi_hat_n[[i]] <- 0
  }
}

if (l_methods$bisection |
    l_methods$bisection_asym |
    l_methods$Bonf_bisec |
    l_methods$percentile_bt1 |
    l_methods$percentile_bt_gemini) {
  
# Define a helper function to run for each setting
run_bootstrap_iterations <- function(setting_index) {
  B <- df_sim_settings$B[setting_index]
  
  # Pre-allocate lists to store results from each bootstrap iteration
  pi_star_hat_b_list <- vector("list", B)
  phi_star_hat_b_list <- vector("list", B)
  y_star_hat_b_list <- vector("list", B)
  se_pred_star_b_list <- vector("list", B)
  y_star_b_list <- vector("list", B) # For bisection methods
  
  # Get parameters for the current setting
  current_phi_n <- l_phi_hat_n[[setting_index]]
  current_phi_m <- l_phi_hat_m[[setting_index]]
  current_pi <- l_pi_hat_vec[[setting_index]]
  current_n <- df_sim_settings$n[setting_index]
  current_k <- df_sim_settings$k[setting_index]
  current_m <- df_sim_settings$m[setting_index]
  current_disp <- as.character(df_sim_settings$dispersion[setting_index])
  current_N_hist <- l_N_hist[setting_index]
  
  # Return NULL if initial pi estimation failed
  # Corrected code with more robust check
  if(is.null(current_pi) || any(is.na(current_pi)) || is.null(current_phi_n) || is.na(current_phi_n)) {
    # If the initial estimates are invalid, we cannot proceed with this simulation setting.
    print(l_pi_hat_vec)
    print(l_phi_hat_m)
    return(NULL)
  }
  
  for (b in 1:B) {
    
    # Initialize a variable to hold the result for this iteration
    x_star_b <- NULL
    count_error_mbt <- 0
    # This loop will continue until x_star_b is successfully assigned a value
    while (is.null(x_star_b)) {
      
      tryCatch({
        # 1. Generate ONE historical bootstrap sample
        x_star_b <- as.data.frame(f_rdirmultinom(
          phi = current_phi_n,
          k = current_k,
          n = current_n,
          v_pi = current_pi
        ))
      }, error = function(e) {
        # This block is executed ONLY if an error occurs in the code above.
        # We print a message to the console to stay informed that a retry is happening.
        message("Main bootstrap x: An error in data generation was caught. Retrying...")
        count_error_mbt <- count_error_mbt + 1
        
        if (count_error_mbt > 10) {
          message("Main bootstrap x: Ten errors occurred and were caught. Retrying with phi*0.99")
          current_phi_n <- current_phi_n*0.99
        }
        # By not assigning anything to x_star_b here, the variable remains NULL,
        # ensuring the while-loop runs again.
      })
    }
    
    # Handle all-zero columns
    for (col in colnames(x_star_b)) {
      if (all(x_star_b[[col]] == 0)) {
        x_star_b[[sample(nrow(x_star_b), 1), col]] <- 1
      }
    }
    
    # 2. Fit model to the bootstrap sample
    model_x_star <- f_fit_mult(df_i = x_star_b, modeltype = "multinomial")
    
    if (!is.null(model_x_star[[1]])) {
      # 3. Calculate bootstrap statistics
      pi_star_hat_b <- f_calc_pi(model_x_star[[1]])
      phi_star_hat_b <- f_estimate_phi(model = model_x_star[[1]], dispersion = current_disp)
      
      pi_star_hat_b_list[[b]] <- pi_star_hat_b
      phi_star_hat_b_list[[b]] <- phi_star_hat_b
      
      y_star_hat_b_list[[b]] <- f_calc_y_hat(pi_star_hat_b, current_m)
      se_pred_star_b_list[[b]] <- f_calc_se_pred_star_internal(phi_star_hat_b, pi_star_hat_b, current_m, current_N_hist)
      
      # 4. Generate ONE future bootstrap sample (if needed by methods)
      if (l_methods$bisection || l_methods$bisection_asym || l_methods$Bonf_bisec || l_methods$percentile_bt1) {
        
        y_star_b <- NULL
        count_error_mbty <- 0
        
        # This loop will continue until y_star_b is successfully assigned a value
        while (is.null(y_star_b)) {
          tryCatch({
            y_star_b <- as.data.frame(f_rdirmultinom(
              phi = current_phi_m,
              v_pi = current_pi,
              k = 1,
              n = current_m
            ))
          }, error = function(e) {
            # This block is executed ONLY if an error occurs in the code above.
            # We print a message to the console to stay informed that a retry is happening.
            message("Main bootstrap y: An error in data generation was caught. Retrying...")
            count_error_mbty <- count_error_mbty + 1
            
            if (count_error_mbty > 10) {
              message("Main bootstrap y: Ten errors occurred and were caught. Retrying with phi*0.99")
              current_phi_m <- current_phi_m*0.99
            }
            
            # By not assigning anything to x_star_b here, the variable remains NULL,
            # ensuring the while-loop runs again.
          })
        }
        
        y_star_b_list[[b]] <- y_star_b
      }
    }
    # The large objects 'x_star_b' and 'model_x_star' are discarded at the end of the loop
  }
  
  # Return a list of all calculated bootstrap results for this setting
  return(list(
    pi_star_hat_b_vec = pi_star_hat_b_list,
    phi_star_hat_b = phi_star_hat_b_list,
    y_star_hat_b_vec = y_star_hat_b_list,
    se_pred_star_b_vec = se_pred_star_b_list,
    y_star_b_vec = y_star_b_list
  ))
}

# Run the iterative bootstrap process for each parameter setting
bootstrap_results <- lapply(seq_len(nrow(df_sim_settings)), run_bootstrap_iterations)

# Unpack the results into the list structures your existing code expects
l_pi_star_hat_b_vec <- lapply(bootstrap_results, `[[`, "pi_star_hat_b_vec")
l_phi_star_hat_b <- lapply(bootstrap_results, `[[`, "phi_star_hat_b")
l_y_star_hat_b_vec <- lapply(bootstrap_results, `[[`, "y_star_hat_b_vec")
l_se_pred_star_b_vec <- lapply(bootstrap_results, `[[`, "se_pred_star_b_vec")
l_y_star_b_vec <- lapply(bootstrap_results, `[[`, "y_star_b_vec")

}

# Explicitly call garbage collection after the big loop
gc()

# create simulation result dataframe
df_sim_res <- cbind(
  df_sim_settings_og, # original simiulation settings dataframe
  data.frame(
    sim = 1,
    C = ncol(m_true_pi_vec), #number of categories
    minpn = min(m_true_pi_vec) * df_sim_settings$n,
    minpnk = min(m_true_pi_vec) * df_sim_settings$n * df_sim_settings$k,
    minp = min(m_true_pi_vec),
    maxp = max(m_true_pi_vec),
    phi_hat = unlist(l_phi_hat),
    czercounts = v_historical_zerocounts,
    c_disphat_m = l_c_phi_hat_m,
    c_disphat_n = l_c_phi_hat_n
  ),
  df_err_warn
)

# use Pointwise Normal Approximation to compute prediction intervals
if (l_methods$pointwise){
  
  # calculate vector of quantiles depending on respective alphas
  v_q_norm <- qnorm((1-df_sim_settings$alpha/2))
  
  # calculate pointwise prediction intervals 
  l_pred_int_pointwise <- mapply(FUN = f_calc_prediction_interval,
                                  y_hat = l_y_hat_vec,
                                  q = v_q_norm,
                                  se_pred = l_se_pred_vec,
                                  alternative = as.list(df_sim_settings$alternative),
                                  SIMPLIFY = FALSE)
  
  
  l_res_pointwise<- lapply(seq_along(l_pred_int_pointwise), function(i) {
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_pointwise[[i]]
    
    # Exclude non-category columns (e.g., 'level') from df_counts
    category_cols <- setdiff(names(df_counts), "level")
    
    # Extract counts for all categories
    counts <- as.numeric(df_counts[1, category_cols])
    
    # Ensure counts and intervals have the same length
    if (length(counts) != nrow(df_intervals)) {
      stop(paste("Mismatch in counts and intervals at index", i))
    }
    
    # Extract lower and upper intervals
    lower <- df_intervals$lower
    upper <- df_intervals$upper
    
    # Check if all counts are within the prediction intervals
    result <- all(counts >= lower & counts <= upper)
    
    return(result)
    
  })
  
  # Check whether whether future obs are higher, lower or within the intervals
  l_cwi_pointwise <- mapply(FUN = f_compare_with_intervals, 
                            data = l_y_vec, 
                            intervals = l_pred_int_pointwise, SIMPLIFY = FALSE)
  
  df_res_pointwise <- as.data.frame(do.call(rbind, l_res_pointwise))
  # rename the columns
  colnames(df_res_pointwise) <- paste0("coverage_pointwise")
  
  # Convert the list of padded vectors into a data frame
  df_cwi_pointwise <- as.data.frame(do.call(rbind, l_cwi_pointwise))
  # rename the columns
  colnames(df_cwi_pointwise) <- paste0("pointwise_", colnames(df_cwi_pointwise))
  
  df_width_stats_pointwise <- do.call(rbind, lapply(l_pred_int_pointwise, f_calc_interval_width_stats))
  colnames(df_width_stats_pointwise) <- paste0("pointwise_", colnames(df_width_stats_pointwise))
  
  df_sim_res <- cbind(df_sim_res, df_res_pointwise, df_cwi_pointwise, df_width_stats_pointwise)
  
}

gc()

# use Bonferroni Adjustment to compute prediction intervals
if (l_methods$Bonferroni){
  
  # list with elements of number of categories
  l_n_cat <- lapply(l_y_hat_vec, function(x) length(x))
  
  # calculate bonferroni adjusted q
  l_q_norm <- mapply(f_adjust_q, 
                      n = l_n_cat,
                      alpha = df_sim_settings$alpha,
                      SIMPLIFY = FALSE)
  
  # calculate prediction intervals with bonferroni adjusted q
  l_pred_int_bonf <- mapply(FUN = f_calc_prediction_interval,
                            y_hat = l_y_hat_vec,
                            q = l_q_norm,
                            se_pred = l_se_pred_vec,
                            alternative = as.list(df_sim_settings$alternative),
                            SIMPLIFY = FALSE)
  
  
  l_res_bonf <- lapply(seq_along(l_pred_int_bonf), function(i) {
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_bonf[[i]]
    
    # Exclude non-category columns (e.g., 'level') from df_counts
    category_cols <- setdiff(names(df_counts), "level")
    
    # Extract counts for all categories
    counts <- as.numeric(df_counts[1, category_cols])
    
    # Ensure counts and intervals have the same length
    if (length(counts) != nrow(df_intervals)) {
      stop(paste("Mismatch in counts and intervals at index", i))
    }
    
    # Extract lower and upper intervals
    lower <- df_intervals$lower
    upper <- df_intervals$upper
    
    # Check if all counts are within the prediction intervals
    result <- all(counts >= lower & counts <= upper)
    
    return(result)
    
  })
  
  # Check whether whether future obs are higher, lower or within the intervals
  l_cwi_bonferoni <- mapply(FUN = f_compare_with_intervals, 
                            data = l_y_vec, 
                            intervals = l_pred_int_bonf, SIMPLIFY = FALSE)
  
  df_res_bonferroni <- as.data.frame(do.call(rbind, l_res_bonf))
  # rename the columns
  colnames(df_res_bonferroni) <- paste0("coverage_bonferroni")
  
  # Convert the list of padded vectors into a data frame
  df_cwi_bonferroni <- as.data.frame(do.call(rbind, l_cwi_bonferoni))
  # rename the columns
  colnames(df_cwi_bonferroni) <- paste0("bonferroni_", colnames(df_cwi_bonferroni))
  
  df_width_stats_bonf <- do.call(rbind, lapply(l_pred_int_bonf, f_calc_interval_width_stats))
  colnames(df_width_stats_bonf) <- paste0("bonferroni_", colnames(df_width_stats_bonf))
  
  df_sim_res <- cbind(df_sim_res, df_res_bonferroni, df_cwi_bonferroni, df_width_stats_bonf)
  
  }

gc()

# use Symmetrical Bisection to compute prediction intervals
if (l_methods$bisection){
  
  # calibrate q with bisection
  l_alpha_calib_sym <- mapply(FUN = f_bisection,
                              y_star_hat = l_y_star_hat_b_vec,
                              se_pred = l_se_pred_star_b_vec,
                              y_star = l_y_star_b_vec,
                              alternative = as.list(df_sim_settings$alternative),
                              quant_min = as.list(df_sim_settings$quant_min),
                              quant_max = as.list(df_sim_settings$quant_max),
                              n_bisec = as.list(df_sim_settings$n_bisec),
                              tol = as.list(df_sim_settings$tol),
                              alpha = as.list(df_sim_settings$alpha),
                              MoreArgs = list(traceplot = FALSE),
                              SIMPLIFY = FALSE)
  
  # calculate prediction intervals with calibrated q
  l_pred_int_bisection_sym <- mapply(FUN = f_calc_prediction_interval,
                                      y_hat = l_y_hat_vec,
                                      q = l_alpha_calib_sym,
                                      se_pred = l_se_pred_vec,
                                      alternative = as.list(df_sim_settings$alternative),
                                      SIMPLIFY = FALSE)
  
  
  l_res_bisection_sym <- lapply(seq_along(l_pred_int_bisection_sym), function(i) {
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_bisection_sym[[i]]
    
    # Exclude non-category columns (e.g., 'level') from df_counts
    category_cols <- setdiff(names(df_counts), "level")
    
    # Extract counts for all categories
    counts <- as.numeric(df_counts[1, category_cols])
    
    # Ensure counts and intervals have the same length
    if (length(counts) != nrow(df_intervals)) {
      stop(paste("Mismatch in counts and intervals at index", i))
    }
    
    # Extract lower and upper intervals
    lower <- df_intervals$lower
    upper <- df_intervals$upper
    
    # Check if all counts are within the prediction intervals
    result <- all(counts >= lower & counts <= upper)
    
    return(result)
    
  })
  
  # Check whether whether future obs are higher, lower or within the intervals
  l_cwi_bisection_sym <- mapply(FUN = f_compare_with_intervals, 
                                data = l_y_vec, 
                                intervals = l_pred_int_bisection_sym, SIMPLIFY = FALSE)
  
  df_res_bisection_sym <- as.data.frame(do.call(rbind, l_res_bisection_sym))
  # rename the columns
  colnames(df_res_bisection_sym) <- paste0("coverage_bisection")
  
  # Convert the list of padded vectors into a data frame
  df_cwi_bisection_sym <- as.data.frame(do.call(rbind, l_cwi_bisection_sym))
  # rename the columns
  colnames(df_cwi_bisection_sym) <- paste0("bisection_", colnames(df_cwi_bisection_sym))
  
  df_width_stats_bisection <- do.call(rbind, lapply(l_pred_int_bisection_sym, f_calc_interval_width_stats))
  colnames(df_width_stats_bisection) <- paste0("bisection_", colnames(df_width_stats_bisection))
  
  df_sim_res <- cbind(df_sim_res, df_res_bisection_sym, df_cwi_bisection_sym, df_width_stats_bisection)
  
}

gc()

# use Asymmetrical Bisection to compute prediction intervals
if (l_methods$bisection_asym){
  
  # calibrate q with bisection for lower border
  l_alpha_calib_asym_lower <- mapply(FUN = f_bisection,
                                      y_star_hat = l_y_star_hat_b_vec,
                                      se_pred = l_se_pred_star_b_vec,
                                      y_star = l_y_star_b_vec,
                                      quant_min = as.list(df_sim_settings$quant_min),
                                      quant_max = as.list(df_sim_settings$quant_max),
                                      n_bisec = as.list(df_sim_settings$n_bisec),
                                      tol = as.list(df_sim_settings$tol),
                                      alpha = as.list(df_sim_settings$alpha/2),
                                      MoreArgs = list(traceplot = FALSE,
                                                      alternative = "lower"),
                                      SIMPLIFY = FALSE)
  
  # calibrate q with bisection for upper boarder
  l_alpha_calib_asym_upper <- mapply(FUN = f_bisection,
                                      y_star_hat = l_y_star_hat_b_vec,
                                      se_pred = l_se_pred_star_b_vec,
                                      y_star = l_y_star_b_vec,
                                      quant_min = as.list(df_sim_settings$quant_min),
                                      quant_max = as.list(df_sim_settings$quant_max),
                                      n_bisec = as.list(df_sim_settings$n_bisec),
                                      tol = as.list(df_sim_settings$tol),
                                      alpha = as.list(df_sim_settings$alpha/2),
                                      MoreArgs = list(traceplot = FALSE,
                                                      alternative = "upper"),
                                      SIMPLIFY = FALSE)
  
  l_alpha_calib_asym_both <- Map(c, l_alpha_calib_asym_lower, l_alpha_calib_asym_upper)
  
  # calculate prediction intervals with calibrated q
  l_pred_int_bisection_asym <- mapply(FUN = f_calc_prediction_interval,
                                      y_hat = l_y_hat_vec,
                                      q = l_alpha_calib_asym_both,
                                      se_pred = l_se_pred_vec,
                                      alternative = as.list(df_sim_settings$alternative),
                                      SIMPLIFY = FALSE)
  
  l_res_bisection_asym <- lapply(seq_along(l_pred_int_bisection_asym), function(i) {
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_bisection_asym[[i]]
    
    # Exclude non-category columns (e.g., 'level') from df_counts
    category_cols <- setdiff(names(df_counts), "level")
    
    # Extract counts for all categories
    counts <- as.numeric(df_counts[1, category_cols])
    
    # Ensure counts and intervals have the same length
    if (length(counts) != nrow(df_intervals)) {
      stop(paste("Mismatch in counts and intervals at index", i))
    }
    
    # Extract lower and upper intervals
    lower <- df_intervals$lower
    upper <- df_intervals$upper
    
    # Check if all counts are within the prediction intervals
    result <- all(counts >= lower & counts <= upper)
    
    return(result)
    
  })
  
  # Check whether whether future obs are higher, lower or within the intervals
  l_cwi_bisection_asym <- mapply(FUN = f_compare_with_intervals, 
                                  data = l_y_vec, 
                                  intervals = l_pred_int_bisection_asym, SIMPLIFY = FALSE)
  
  df_res_bisection_asym <- as.data.frame(do.call(rbind, l_res_bisection_asym))
  # rename the columns
  colnames(df_res_bisection_asym) <- paste0("coverage_bisection_asym")
  
  # Convert the list of padded vectors into a data frame
  df_cwi_bisection_asym <- as.data.frame(do.call(rbind, l_cwi_bisection_asym))
  # rename the columns
  colnames(df_cwi_bisection_asym) <- paste0("bisection_asym_", colnames(df_cwi_bisection_asym))
  
  df_width_stats_bisection_asym <- do.call(rbind, lapply(l_pred_int_bisection_asym, f_calc_interval_width_stats))
  colnames(df_width_stats_bisection_asym) <- paste0("bisection_asym_", colnames(df_width_stats_bisection_asym))
  
  df_sim_res <- cbind(df_sim_res, df_res_bisection_asym, df_cwi_bisection_asym, df_width_stats_bisection_asym)
  
}

gc()

# use Bonferroni with bisection method to compute prediction intervals
if (l_methods$Bonf_bisec){
  
  # list with elements of number of categories
  l_n_cat <- lapply(l_y_hat_vec, function(x) length(x))
  
  l_alpha_calib_bonf_bisec <- mapply(FUN = f_single_bisec,
                                      n_cat = l_n_cat,
                                      y_star_hat = l_y_star_hat_b_vec,
                                      se_pred = l_se_pred_star_b_vec,
                                      y_star = l_y_star_b_vec,
                                      quant_min = as.list(df_sim_settings$quant_min),
                                      quant_max = as.list(df_sim_settings$quant_max),
                                      n_bisec = as.list(df_sim_settings$n_bisec),
                                      tol = as.list(df_sim_settings$tol),
                                      alpha = as.list(df_sim_settings$alpha),
                                      alternative = as.list(df_sim_settings$alternative),
                                      MoreArgs = list(traceplot = FALSE),
                                      SIMPLIFY = FALSE)
  
  # calculate prediction intervals with bonferroni adjusted q
  l_pred_int_bonf_bisec <- mapply(FUN = f_calc_prediction_interval,
                                  y_hat = l_y_hat_vec,
                                  q = l_alpha_calib_bonf_bisec,
                                  se_pred = l_se_pred_vec,
                                  alternative = as.list(df_sim_settings$alternative),
                                  SIMPLIFY = FALSE)
  
  
  l_res_bonf_bisec <- lapply(seq_along(l_pred_int_bonf_bisec), function(i) {
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_bonf_bisec[[i]]
    
    # Exclude non-category columns (e.g., 'level') from df_counts
    category_cols <- setdiff(names(df_counts), "level")
    
    # Extract counts for all categories
    counts <- as.numeric(df_counts[1, category_cols])
    
    # Ensure counts and intervals have the same length
    if (length(counts) != nrow(df_intervals)) {
      stop(paste("Mismatch in counts and intervals at index", i))
    }
    
    # Extract lower and upper intervals
    lower <- df_intervals$lower
    upper <- df_intervals$upper
    
    # Check if all counts are within the prediction intervals
    result <- all(counts >= lower & counts <= upper)
    
    return(result)
    
  })
  
  # Check whether whether future obs are higher, lower or within the intervals
  l_cwi_bonf_bisec <- mapply(FUN = f_compare_with_intervals, 
                              data = l_y_vec, 
                              intervals = l_pred_int_bonf_bisec, SIMPLIFY = FALSE)
  
  df_res_bonf_bisec <- as.data.frame(do.call(rbind, l_res_bonf_bisec))
  # rename the columns
  colnames(df_res_bonf_bisec) <- paste0("coverage_bonf_bisec")
  
  # Convert the list of padded vectors into a data frame
  df_cwi_bonf_bisec <- as.data.frame(do.call(rbind, l_cwi_bonf_bisec))
  # rename the columns
  colnames(df_cwi_bonf_bisec) <- paste0("bonf_bisec_", colnames(df_cwi_bonf_bisec))
  
  df_width_stats_bonf_bisec <- do.call(rbind, lapply(l_pred_int_bonf_bisec, f_calc_interval_width_stats))
  colnames(df_width_stats_bonf_bisec) <- paste0("bonf_bisec_", colnames(df_width_stats_bonf_bisec))
  
  df_sim_res <- cbind(df_sim_res, df_res_bonf_bisec, df_cwi_bonf_bisec, df_width_stats_bonf_bisec)
  
}

gc()

# Frank Aufschrieb
# y_star: bt_fut_y_star
# y_star_hat: bt_hist_y_star_hat
# se: bt_hist_se
if (l_methods$percentile_bt1) {
  
  l_df_zc <- mapply(FUN = f_calc_zc_main,
                    l_y_star = l_y_star_b_vec,
                    l_y_star_hat = l_y_star_hat_b_vec,
                    l_se = l_se_pred_star_b_vec,
                    MoreArgs = list(n_cat = ncol(m_true_pi_vec)),
                    SIMPLIFY = FALSE)
  
  l_scs_zc <- lapply(l_df_zc, function(df){
    abs(SCSrank(as.matrix(df))$conf.int)
  })
  
  
  
  # calculate prediction intervals with rectangular simultaneous confidence set
  l_pred_int_percentile_bt1 <- mapply(FUN = f_calc_prediction_interval,
                                      y_hat = l_y_hat_vec,
                                      q = l_scs_zc,
                                      se_pred = l_se_pred_vec,
                                      alternative = as.list(df_sim_settings$alternative),
                                      SIMPLIFY = FALSE)
  
  
  l_res_percentile_bt1 <- lapply(seq_along(l_pred_int_percentile_bt1), function(i) {
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_percentile_bt1[[i]]
    
    # Exclude non-category columns (e.g., 'level') from df_counts
    category_cols <- setdiff(names(df_counts), "level")
    
    # Extract counts for all categories
    counts <- as.numeric(df_counts[1, category_cols])
    
    # Ensure counts and intervals have the same length
    if (length(counts) != nrow(df_intervals)) {
      stop(paste("Mismatch in counts and intervals at index", i))
    }
    
    # Extract lower and upper intervals
    lower <- df_intervals$lower
    upper <- df_intervals$upper
    
    # Check if all counts are within the prediction intervals
    result <- all(counts >= lower & counts <= upper)
    
    return(result)
    
  })
  
  # Check whether whether future obs are higher, lower or within the intervals
  l_cwi_percentile_bt1 <- mapply(FUN = f_compare_with_intervals,
                                  data = l_y_vec,
                                  intervals = l_pred_int_percentile_bt1, SIMPLIFY = FALSE)
  
  df_res_percentile_bt1 <- as.data.frame(do.call(rbind, l_res_percentile_bt1))
  # rename the columns
  colnames(df_res_percentile_bt1) <- paste0("coverage_percentile_bt1")
  
  # Convert the list of padded vectors into a data frame
  df_cwi_percentile_bt1 <- as.data.frame(do.call(rbind, l_cwi_percentile_bt1))
  # rename the columns
  colnames(df_cwi_percentile_bt1) <- paste0("percentile_bt1_", colnames(df_cwi_percentile_bt1))
  
  df_width_stats_percentile_bt1 <- do.call(rbind, lapply(l_pred_int_percentile_bt1, f_calc_interval_width_stats))
  colnames(df_width_stats_percentile_bt1) <- paste0("percentile_bt1_", colnames(df_width_stats_percentile_bt1))
  
  df_sim_res <- cbind(df_sim_res, df_res_percentile_bt1, df_cwi_percentile_bt1, df_width_stats_percentile_bt1)
  
}

gc()

# use Multivariate Normal Approximation to compute prediction intervals
if (l_methods$mvn){
  
  # Calculate the MVN-adjusted quantile for each setting
  l_q_mvn <- mapply(FUN = f_calc_q_mvn,
                    phi_hat = l_phi_hat,
                    pi_hat = l_pi_hat_vec,
                    m = as.list(df_sim_settings$m),
                    N_hist = l_N_hist,
                    alpha = df_sim_settings$alpha,
                    SIMPLIFY = FALSE)
  
  # calculate prediction intervals with MVN-adjusted q
  l_pred_int_mvn <- mapply(FUN = f_calc_prediction_interval,
                            y_hat = l_y_hat_vec,
                            q = l_q_mvn,
                            se_pred = l_se_pred_vec,
                            alternative = as.list(df_sim_settings$alternative),
                            SIMPLIFY = FALSE)
  
  
  l_res_mvn <- lapply(seq_along(l_pred_int_mvn), function(i) {
    # Handle cases where interval calculation failed
    if (is.null(l_pred_int_mvn[[i]])) return(NA)
    
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_mvn[[i]]
    
    category_cols <- setdiff(names(df_counts), "level")
    counts <- as.numeric(df_counts[1, category_cols])
    
    if (length(counts) != nrow(df_intervals)) {
      stop(paste("Mismatch in counts and intervals at index", i))
    }
    
    lower <- df_intervals$lower
    upper <- df_intervals$upper
    
    # Check if all counts are within the prediction intervals
    all(counts >= lower & counts <= upper)
  })
  
  # Check whether future obs are higher, lower or within the intervals
  l_cwi_mvn <- mapply(FUN = f_compare_with_intervals, 
                      data = l_y_vec, 
                      intervals = l_pred_int_mvn, SIMPLIFY = FALSE)
  
  df_res_mvn <- as.data.frame(do.call(rbind, l_res_mvn))
  colnames(df_res_mvn) <- "coverage_mvn"
  
  df_cwi_mvn <- as.data.frame(do.call(rbind, l_cwi_mvn))
  colnames(df_cwi_mvn) <- paste0("mvn_", colnames(df_cwi_mvn))
  
  df_width_stats_mvn <- do.call(rbind, lapply(l_pred_int_mvn, f_calc_interval_width_stats))
  colnames(df_width_stats_mvn) <- paste0("mvn_", colnames(df_width_stats_mvn))
  
  df_sim_res <- cbind(df_sim_res, df_res_mvn, df_cwi_mvn, df_width_stats_mvn)
  
}


# use Percentile Bootstrap (Gemini Corrected Method) to compute prediction intervals
if (l_methods$percentile_bt_gemini) {
  
  # This part remains the same: calculate the matrix of pivotal quantities (standardized residuals)
  # where each row is a bootstrap replication and each column is a category.
  l_df_zc <- mapply(FUN = f_calc_zc_main,
                    l_y_star = l_y_star_b_vec,
                    l_y_star_hat = l_y_star_hat_b_vec,
                    l_se = l_se_pred_star_b_vec,
                    MoreArgs = list(n_cat = ncol(m_true_pi_vec)),
                    SIMPLIFY = FALSE)
  
  # --- CORRECTED LOGIC STARTS HERE ---
  # Instead of SCSrank, we now find the quantile of the maximum absolute residuals.
  l_q_percentile <- lapply(seq_along(l_df_zc), function(i) {
    
    df_zc <- l_df_zc[[i]]
    
    # Handle cases with errors or insufficient data
    if(is.null(df_zc) || any(is.na(df_zc)) || nrow(df_zc) < 2) {
      return(NA)
    }
    
    # For each bootstrap replication (row), find the maximum absolute standardized residual.
    max_abs_residuals <- apply(abs(df_zc), 1, max)
    
    # The critical value 'q' is the (1-alpha) quantile of this distribution of maximums.
    q_val <- quantile(max_abs_residuals, probs = 1 - df_sim_settings$alpha[i], na.rm = TRUE)
    
    return(q_val)
  })
  # --- CORRECTED LOGIC ENDS HERE ---
  
  
  # Calculate prediction intervals with the new single critical value q
  l_pred_int_percentile_bt_gemini <- mapply(FUN = f_calc_prediction_interval,
                                            y_hat = l_y_hat_vec,
                                            q = l_q_percentile, # Use the correctly calculated q
                                            se_pred = l_se_pred_vec,
                                            alternative = as.list(df_sim_settings$alternative),
                                            SIMPLIFY = FALSE)
  
  
  # This part remains the same: Calculate coverage using the new intervals
  l_res_percentile_bt_gemini <- lapply(seq_along(l_pred_int_percentile_bt_gemini), function(i) {
    if (is.null(l_pred_int_percentile_bt_gemini[[i]])) return(NA)
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_percentile_bt_gemini[[i]]
    
    category_cols <- setdiff(names(df_counts), "level")
    counts <- as.numeric(df_counts[1, category_cols])
    
    if (length(counts) != nrow(df_intervals)) {
      stop(paste("Mismatch in counts and intervals at index", i))
    }
    
    lower <- df_intervals$lower
    upper <- df_intervals$upper
    
    all(counts >= lower & counts <= upper)
  })
  
  # This part also remains the same
  l_cwi_percentile_bt_gemini <- mapply(FUN = f_compare_with_intervals,
                                        data = l_y_vec,
                                        intervals = l_pred_int_percentile_bt_gemini, SIMPLIFY = FALSE)
  
  df_res_percentile_bt_gemini <- as.data.frame(do.call(rbind, l_res_percentile_bt_gemini))
  colnames(df_res_percentile_bt_gemini) <- paste0("coverage_percentile_bt_gemini")
  
  df_cwi_percentile_bt_gemini <- as.data.frame(do.call(rbind, l_cwi_percentile_bt_gemini))
  colnames(df_cwi_percentile_bt_gemini) <- paste0("percentile_bt_gemini_", colnames(df_cwi_percentile_bt_gemini))
  
  df_width_stats_percentile_bt_gemini <- do.call(rbind, lapply(l_pred_int_percentile_bt_gemini, f_calc_interval_width_stats))
  colnames(df_width_stats_percentile_bt_gemini) <- paste0("percentile_bt_gemini_", colnames(df_width_stats_percentile_bt_gemini))
  
  df_sim_res <- cbind(df_sim_res, df_res_percentile_bt_gemini, df_cwi_percentile_bt_gemini, df_width_stats_percentile_bt_gemini)
  
}

gc()

# use Non-Parametric Bootstrap to compute prediction intervals
if (l_methods$nonparametric_bt) {
  
  # Helper function to perform a single bootstrap iteration for one simulation setting.
  # This function encapsulates the entire process for one bootstrap sample, from resampling
  # to the calculation of the pivotal quantity, to avoid storing large intermediate objects.
  perform_bootstrap_iteration <- function(df_hist, sim_setting, N_hist_i) {
    
    # 1. Generate one resampled dataset
    # We set B=1 to get a single bootstrap sample.
    x_star_b <- f_bt_mult_dat_nonparametric(df_hist = df_hist, B = 1)[[1]]
    
    # Correct for columns with only zeros in the bootstrap sample, which would prevent model fitting.
    # A single count is added to a randomly selected zero-entry row in that column.
    for (col in colnames(x_star_b)) {
      if (all(x_star_b[[col]] == 0)) {
        zero_indices <- which(x_star_b[[col]] == 0)
        random_row <- sample(zero_indices, 1)
        x_star_b[[col]][random_row] <- 1
      }
    }
    
    # 2. Fit model to the resampled dataset
    model_x_star <- f_fit_mult(df_i = x_star_b, modeltype = "multinomial")
    
    # If model fitting fails for the bootstrap sample, return NULL.
    if (is.null(model_x_star[[1]])) return(NULL)
    
    # 3. Estimate parameters from the bootstrap sample
    pi_star_hat <- f_calc_pi(model_x_star[[1]])
    phi_star_hat <- f_estimate_phi(model = model_x_star[[1]], dispersion = as.character(sim_setting$dispersion))
    
    # Handle cases where parameter estimation might fail.
    if (is.null(pi_star_hat) || is.null(phi_star_hat)) return(NULL)
    
    # Adjust the estimated dispersion parameter phi. These functions expect lists.
    phi_star_hat <- f_replace_phi_below_one(list(phi_star_hat))[[1]]
    phi_star_hat <- f_replace_phi_above_val(phi = list(phi_star_hat), val = sim_setting$m)[[1]]
    
    if (is.na(phi_star_hat)) return(NULL)
    
    # 4. Generate a future observation y_star based on the bootstrap estimates.
    # A while-loop with tryCatch is used to handle potential errors in f_rdirmultinom.
    y_star <- NULL
    
    count_error_npbt <- 0
    
    while (is.null(y_star)) {
      tryCatch({
        y_star <- f_rdirmultinom(phi = phi_star_hat, v_pi = pi_star_hat, k = 1, n = sim_setting$m)
      }, error = function(e) {
        message("np_bootstrap: An error in data generation was caught. Retrying...")
        count_error_npbt <- count_error_npbt + 1
        
        if (count_error_npbt > 10) {
          message("np_bootstrap: Ten errors occurred and were caught. Retrying with phi*0.99")
          phi_star_hat <- phi_star_hat*0.99
        }
      })
    }
    y_star <- as.data.frame(y_star)
    
    # 5. Calculate predicted value and standard error for the bootstrap sample.
    y_star_hat <- f_calc_y_hat(pi_star_hat, sim_setting$m)
    
    # We call f_calc_se_pred_star with lists of length 1 to get the SE for this single bootstrap sample.
    se_pred_star <- f_calc_se_pred_star_internal(phi_star_hat, pi_star_hat, sim_setting$m, N_hist_i)
    
    # 6. Calculate and return the pivotal quantity for this bootstrap sample.
    zc <- (y_star - y_star_hat) / se_pred_star
    return(as.numeric(zc))
  }
  
  # Number of simulation settings to iterate over
  num_settings <- nrow(df_sim_settings)
  
  # List to store the critical values from SCSrank for each setting
  l_scs_zc_np <- vector("list", num_settings)
  
  # Loop over each simulation setting
  for (i in 1:num_settings) {
    
    # Extract data for the current simulation setting
    df_hist_i <- l_x_k_vec[[i]]
    sim_setting_i <- df_sim_settings[i, ]
    N_hist_i <- l_N_hist[[i]]
    B_i <- sim_setting_i$B # Number of bootstrap replicates
    
    # Replicate the bootstrap iteration B times for the current setting.
    # 'replicate' returns a matrix where each column is a pivotal quantity vector from one iteration.
    # This is far more memory-efficient than storing all intermediate objects.
    zc_matrix <- tryCatch({
      replicate(B_i, perform_bootstrap_iteration(df_hist_i, sim_setting_i, N_hist_i))
    }, error = function(e) {
      warning(paste("Bootstrap replication failed for setting", i, ":", e$message))
      return(NULL)
    })
    
    # Continue to next setting if bootstrap replication failed
    if (is.null(zc_matrix)) {
      l_scs_zc_np[[i]] <- NULL
      next
    }
    
    # Filter out NULL results from failed single iterations and transpose the matrix
    # so that each row represents one bootstrap observation.
    zc_list <- as.list(data.frame(zc_matrix))
    valid_zc <- zc_list[!sapply(zc_list, is.null)]
    
    if (length(valid_zc) < 2) { # Need at least 2 valid replicates for SCSrank
      l_scs_zc_np[[i]] <- NULL
      next
    }
    
    df_zc <- do.call(rbind, valid_zc)
    
    # Calculate the critical value using SCSrank
    l_scs_zc_np[[i]] <- tryCatch({
      abs(SCSrank(as.matrix(df_zc))$conf.int)
    }, error = function(e) {
      warning(paste("SCSrank failed for setting", i, ":", e$message))
      return(NULL)
    })
  }
  
  # The rest of the logic remains the same as in your original code.
  # It uses the calculated critical values to construct prediction intervals and assess coverage.
  
  l_pred_int_nonparametric_bt <- mapply(FUN = f_calc_prediction_interval,
                                        y_hat = l_y_hat_vec,
                                        q = l_scs_zc_np,
                                        se_pred = l_se_pred_vec,
                                        alternative = as.list(df_sim_settings$alternative),
                                        SIMPLIFY = FALSE)
  
  # Calculate coverage
  l_res_nonparametric_bt <- lapply(seq_along(l_pred_int_nonparametric_bt), function(i) {
    if (is.null(l_pred_int_nonparametric_bt[[i]])) return(NA)
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_nonparametric_bt[[i]]
    category_cols <- setdiff(names(df_counts), "level")
    counts <- as.numeric(df_counts[1, category_cols])
    if (length(counts) != nrow(df_intervals)) stop(paste("Mismatch in counts and intervals at index", i))
    lower <- df_intervals$lower
    upper <- df_intervals$upper
    all(counts >= lower & counts <= upper)
  })
  
  l_cwi_nonparametric_bt <- mapply(FUN = f_compare_with_intervals,
                                    data = l_y_vec,
                                    intervals = l_pred_int_nonparametric_bt, SIMPLIFY = FALSE)
  
  df_res_nonparametric_bt <- as.data.frame(do.call(rbind, l_res_nonparametric_bt))
  colnames(df_res_nonparametric_bt) <- "coverage_nonparametric_bt"
  
  df_cwi_nonparametric_bt <- as.data.frame(do.call(rbind, l_cwi_nonparametric_bt))
  colnames(df_cwi_nonparametric_bt) <- paste0("nonparametric_bt_", colnames(df_cwi_nonparametric_bt))
  
  df_width_stats_nonparametric_bt <- do.call(rbind, lapply(l_pred_int_nonparametric_bt, f_calc_interval_width_stats))
  colnames(df_width_stats_nonparametric_bt) <- paste0("nonparametric_bt_", colnames(df_width_stats_nonparametric_bt))
  
  df_sim_res <- cbind(df_sim_res, df_res_nonparametric_bt, df_cwi_nonparametric_bt, df_width_stats_nonparametric_bt)
  
}

# use Bayesian Hierarchical MCMC Model (Gamma Prior)
if (l_methods$bayesian_mcmc_gamma) {
  
  # l_pred_int_bayesian_gamma is now a list of lists.
  # Each element is list(mean_centered = df, marginal = df, scs_rank = df)
  l_pred_int_bayesian_gamma <- mapply(
    FUN = f_calc_pi_bayesian_mcmc,
    df_hist = l_x_k_vec,
    m = as.list(df_sim_settings$m),
    alpha = as.list(df_sim_settings$alpha),
    MoreArgs = list(stan_model = mod_gamma),
    SIMPLIFY = FALSE
  )
  
  # Extract the individual lists of data frames
  l_pred_int_bayesian_gamma_mean <- lapply(l_pred_int_bayesian_gamma, `[[`, "mean_centered")
  l_pred_int_bayesian_gamma_marg <- lapply(l_pred_int_bayesian_gamma, `[[`, "marginal")
  l_pred_int_bayesian_gamma_scs  <- lapply(l_pred_int_bayesian_gamma, `[[`, "scs_rank")
  
  # --- 1. Bayesian MCMC Gamma (Mean Centered) ---
  l_res_bayesian_gamma_mean <- lapply(seq_along(l_pred_int_bayesian_gamma_mean), function(i) {
    if (is.null(l_pred_int_bayesian_gamma_mean[[i]])) return(NA)
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_bayesian_gamma_mean[[i]]
    category_cols <- setdiff(names(df_counts), "level")
    counts <- as.numeric(df_counts[1, category_cols])
    if (length(counts) != nrow(df_intervals)) stop(paste("Mismatch at index", i))
    all(counts >= df_intervals$lower & counts <= df_intervals$upper)
  })
  l_cwi_bayesian_gamma_mean <- mapply(FUN = f_compare_with_intervals, data = l_y_vec, intervals = l_pred_int_bayesian_gamma_mean, SIMPLIFY = FALSE)
  df_res_bayesian_gamma_mean <- as.data.frame(do.call(rbind, l_res_bayesian_gamma_mean))
  colnames(df_res_bayesian_gamma_mean) <- "coverage_bayesian_mcmc_gamma_mean"
  df_cwi_bayesian_gamma_mean <- as.data.frame(do.call(rbind, l_cwi_bayesian_gamma_mean))
  colnames(df_cwi_bayesian_gamma_mean) <- paste0("bayesian_mcmc_gamma_mean_", colnames(df_cwi_bayesian_gamma_mean))
  df_width_stats_bayesian_gamma_mean <- do.call(rbind, lapply(l_pred_int_bayesian_gamma_mean, f_calc_interval_width_stats))
  colnames(df_width_stats_bayesian_gamma_mean) <- paste0("bayesian_mcmc_gamma_mean_", colnames(df_width_stats_bayesian_gamma_mean))
  df_sim_res <- cbind(df_sim_res, df_res_bayesian_gamma_mean, df_cwi_bayesian_gamma_mean, df_width_stats_bayesian_gamma_mean)
  
  # --- 2. Bayesian MCMC Gamma (Marginal) ---
  l_res_bayesian_gamma_marg <- lapply(seq_along(l_pred_int_bayesian_gamma_marg), function(i) {
    if (is.null(l_pred_int_bayesian_gamma_marg[[i]])) return(NA)
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_bayesian_gamma_marg[[i]]
    category_cols <- setdiff(names(df_counts), "level")
    counts <- as.numeric(df_counts[1, category_cols])
    if (length(counts) != nrow(df_intervals)) stop(paste("Mismatch at index", i))
    all(counts >= df_intervals$lower & counts <= df_intervals$upper)
  })
  l_cwi_bayesian_gamma_marg <- mapply(FUN = f_compare_with_intervals, data = l_y_vec, intervals = l_pred_int_bayesian_gamma_marg, SIMPLIFY = FALSE)
  df_res_bayesian_gamma_marg <- as.data.frame(do.call(rbind, l_res_bayesian_gamma_marg))
  colnames(df_res_bayesian_gamma_marg) <- "coverage_bayesian_mcmc_gamma_marg"
  df_cwi_bayesian_gamma_marg <- as.data.frame(do.call(rbind, l_cwi_bayesian_gamma_marg))
  colnames(df_cwi_bayesian_gamma_marg) <- paste0("bayesian_mcmc_gamma_marg_", colnames(df_cwi_bayesian_gamma_marg))
  df_width_stats_bayesian_gamma_marg <- do.call(rbind, lapply(l_pred_int_bayesian_gamma_marg, f_calc_interval_width_stats))
  colnames(df_width_stats_bayesian_gamma_marg) <- paste0("bayesian_mcmc_gamma_marg_", colnames(df_width_stats_bayesian_gamma_marg))
  df_sim_res <- cbind(df_sim_res, df_res_bayesian_gamma_marg, df_cwi_bayesian_gamma_marg, df_width_stats_bayesian_gamma_marg)
  
  # --- 3. Bayesian MCMC Gamma (SCS Rank) ---
  l_res_bayesian_gamma_scs <- lapply(seq_along(l_pred_int_bayesian_gamma_scs), function(i) {
    if (is.null(l_pred_int_bayesian_gamma_scs[[i]])) return(NA)
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_bayesian_gamma_scs[[i]]
    category_cols <- setdiff(names(df_counts), "level")
    counts <- as.numeric(df_counts[1, category_cols])
    if (length(counts) != nrow(df_intervals)) stop(paste("Mismatch at index", i))
    all(counts >= df_intervals$lower & counts <= df_intervals$upper)
  })
  l_cwi_bayesian_gamma_scs <- mapply(FUN = f_compare_with_intervals, data = l_y_vec, intervals = l_pred_int_bayesian_gamma_scs, SIMPLIFY = FALSE)
  df_res_bayesian_gamma_scs <- as.data.frame(do.call(rbind, l_res_bayesian_gamma_scs))
  colnames(df_res_bayesian_gamma_scs) <- "coverage_bayesian_mcmc_gamma_scs"
  df_cwi_bayesian_gamma_scs <- as.data.frame(do.call(rbind, l_cwi_bayesian_gamma_scs))
  colnames(df_cwi_bayesian_gamma_scs) <- paste0("bayesian_mcmc_gamma_scs_", colnames(df_cwi_bayesian_gamma_scs))
  df_width_stats_bayesian_gamma_scs <- do.call(rbind, lapply(l_pred_int_bayesian_gamma_scs, f_calc_interval_width_stats))
  colnames(df_width_stats_bayesian_gamma_scs) <- paste0("bayesian_mcmc_gamma_scs_", colnames(df_width_stats_bayesian_gamma_scs))
  df_sim_res <- cbind(df_sim_res, df_res_bayesian_gamma_scs, df_cwi_bayesian_gamma_scs, df_width_stats_bayesian_gamma_scs)
}

# use Bayesian Hierarchical MCMC Model (Cauchy Prior)
if (l_methods$bayesian_mcmc_cauchy) {
  
  # l_pred_int_bayesian_cauchy is now a list of lists.
  l_pred_int_bayesian_cauchy <- mapply(
    FUN = f_calc_pi_bayesian_mcmc,
    df_hist = l_x_k_vec,
    m = as.list(df_sim_settings$m),
    alpha = as.list(df_sim_settings$alpha),
    MoreArgs = list(stan_model = mod_cauchy),
    SIMPLIFY = FALSE
  )
  
  # Extract the individual lists of data frames
  l_pred_int_bayesian_cauchy_mean <- lapply(l_pred_int_bayesian_cauchy, `[[`, "mean_centered")
  l_pred_int_bayesian_cauchy_marg <- lapply(l_pred_int_bayesian_cauchy, `[[`, "marginal")
  l_pred_int_bayesian_cauchy_scs  <- lapply(l_pred_int_bayesian_cauchy, `[[`, "scs_rank")
  
  # --- 1. Bayesian MCMC Cauchy (Mean Centered) ---
  l_res_bayesian_cauchy_mean <- lapply(seq_along(l_pred_int_bayesian_cauchy_mean), function(i) {
    if (is.null(l_pred_int_bayesian_cauchy_mean[[i]])) return(NA)
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_bayesian_cauchy_mean[[i]]
    category_cols <- setdiff(names(df_counts), "level")
    counts <- as.numeric(df_counts[1, category_cols])
    if (length(counts) != nrow(df_intervals)) stop(paste("Mismatch at index", i))
    all(counts >= df_intervals$lower & counts <= df_intervals$upper)
  })
  l_cwi_bayesian_cauchy_mean <- mapply(FUN = f_compare_with_intervals, data = l_y_vec, intervals = l_pred_int_bayesian_cauchy_mean, SIMPLIFY = FALSE)
  df_res_bayesian_cauchy_mean <- as.data.frame(do.call(rbind, l_res_bayesian_cauchy_mean))
  colnames(df_res_bayesian_cauchy_mean) <- "coverage_bayesian_mcmc_cauchy_mean"
  df_cwi_bayesian_cauchy_mean <- as.data.frame(do.call(rbind, l_cwi_bayesian_cauchy_mean))
  colnames(df_cwi_bayesian_cauchy_mean) <- paste0("bayesian_mcmc_cauchy_mean_", colnames(df_cwi_bayesian_cauchy_mean))
  df_width_stats_bayesian_cauchy_mean <- do.call(rbind, lapply(l_pred_int_bayesian_cauchy_mean, f_calc_interval_width_stats))
  colnames(df_width_stats_bayesian_cauchy_mean) <- paste0("bayesian_mcmc_cauchy_mean_", colnames(df_width_stats_bayesian_cauchy_mean))
  df_sim_res <- cbind(df_sim_res, df_res_bayesian_cauchy_mean, df_cwi_bayesian_cauchy_mean, df_width_stats_bayesian_cauchy_mean)
  
  # --- 2. Bayesian MCMC Cauchy (Marginal) ---
  l_res_bayesian_cauchy_marg <- lapply(seq_along(l_pred_int_bayesian_cauchy_marg), function(i) {
    if (is.null(l_pred_int_bayesian_cauchy_marg[[i]])) return(NA)
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_bayesian_cauchy_marg[[i]]
    category_cols <- setdiff(names(df_counts), "level")
    counts <- as.numeric(df_counts[1, category_cols])
    if (length(counts) != nrow(df_intervals)) stop(paste("Mismatch at index", i))
    all(counts >= df_intervals$lower & counts <= df_intervals$upper)
  })
  l_cwi_bayesian_cauchy_marg <- mapply(FUN = f_compare_with_intervals, data = l_y_vec, intervals = l_pred_int_bayesian_cauchy_marg, SIMPLIFY = FALSE)
  df_res_bayesian_cauchy_marg <- as.data.frame(do.call(rbind, l_res_bayesian_cauchy_marg))
  colnames(df_res_bayesian_cauchy_marg) <- "coverage_bayesian_mcmc_cauchy_marg"
  df_cwi_bayesian_cauchy_marg <- as.data.frame(do.call(rbind, l_cwi_bayesian_cauchy_marg))
  colnames(df_cwi_bayesian_cauchy_marg) <- paste0("bayesian_mcmc_cauchy_marg_", colnames(df_cwi_bayesian_cauchy_marg))
  df_width_stats_bayesian_cauchy_marg <- do.call(rbind, lapply(l_pred_int_bayesian_cauchy_marg, f_calc_interval_width_stats))
  colnames(df_width_stats_bayesian_cauchy_marg) <- paste0("bayesian_mcmc_cauchy_marg_", colnames(df_width_stats_bayesian_cauchy_marg))
  df_sim_res <- cbind(df_sim_res, df_res_bayesian_cauchy_marg, df_cwi_bayesian_cauchy_marg, df_width_stats_bayesian_cauchy_marg)
  
  # --- 3. Bayesian MCMC Cauchy (SCS Rank) ---
  l_res_bayesian_cauchy_scs <- lapply(seq_along(l_pred_int_bayesian_cauchy_scs), function(i) {
    if (is.null(l_pred_int_bayesian_cauchy_scs[[i]])) return(NA)
    df_counts <- l_y_vec[[i]]
    df_intervals <- l_pred_int_bayesian_cauchy_scs[[i]]
    category_cols <- setdiff(names(df_counts), "level")
    counts <- as.numeric(df_counts[1, category_cols])
    if (length(counts) != nrow(df_intervals)) stop(paste("Mismatch at index", i))
    all(counts >= df_intervals$lower & counts <= df_intervals$upper)
  })
  l_cwi_bayesian_cauchy_scs <- mapply(FUN = f_compare_with_intervals, data = l_y_vec, intervals = l_pred_int_bayesian_cauchy_scs, SIMPLIFY = FALSE)
  df_res_bayesian_cauchy_scs <- as.data.frame(do.call(rbind, l_res_bayesian_cauchy_scs))
  colnames(df_res_bayesian_cauchy_scs) <- "coverage_bayesian_mcmc_cauchy_scs"
  df_cwi_bayesian_cauchy_scs <- as.data.frame(do.call(rbind, l_cwi_bayesian_cauchy_scs))
  colnames(df_cwi_bayesian_cauchy_scs) <- paste0("bayesian_mcmc_cauchy_scs_", colnames(df_cwi_bayesian_cauchy_scs))
  df_width_stats_bayesian_cauchy_scs <- do.call(rbind, lapply(l_pred_int_bayesian_cauchy_scs, f_calc_interval_width_stats))
  colnames(df_width_stats_bayesian_cauchy_scs) <- paste0("bayesian_mcmc_cauchy_scs_", colnames(df_width_stats_bayesian_cauchy_scs))
  df_sim_res <- cbind(df_sim_res, df_res_bayesian_cauchy_scs, df_cwi_bayesian_cauchy_scs, df_width_stats_bayesian_cauchy_scs)
}

df_sim_res

  l_pred_int_bayesian_cauchy <- mapply(
    FUN = f_calc_pi_bayesian_mcmc,
    df_hist = l_x_k_vec,
    m = as.list(df_sim_settings$m),
    alpha = as.list(df_sim_settings$alpha),
    MoreArgs = list(stan_model = mod_cauchy),
    SIMPLIFY = FALSE
  )

df_hist <- l_x_k_vec[[1]]
m <- df_sim_settings$m[[1]]
alpha <- df_sim_settings$alpha[[1]]
stan_model = mod_cauchy

# This function now calculates THREE types of intervals from a single MCMC run:
# 1. mean_centered: The original simultaneous interval based on max deviation from the mean.
# 2. marginal: Pointwise (non-simultaneous) intervals using simple quantiles.
# 3. scs_rank: A simultaneous interval using the SCSrank method on the posterior samples.
#
# @return A list containing three data frames: $mean_centered, $marginal, $scs_rank


# Prepare data for Stan
category_cols <- grep("^V", names(df_hist))
if (length(category_cols) == 0) {
  category_cols <- 1:(ncol(df_hist)-1)
}
X_hist <- as.matrix(df_hist[, category_cols, drop = FALSE])

stan_data <- list(
  K = nrow(X_hist),
  C = ncol(X_hist),
  n_k = rowSums(X_hist),
  X = X_hist,
  m = m
)

# Fit the model using MCMC
fit <- suppressMessages(
  stan_model$sample(
    data = stan_data,
    chains = 4,
    iter_warmup = 1250,
    iter_sampling = 2500,
    refresh = 0
  )
)

# Extract posterior predictive samples
post_pred_samples <- fit$draws("y_pred", format = "matrix")

# --- 1. Method 1: Mean-Centered Simultaneous Interval (Original) ---
pi_mean_centered <- f_construct_simultaneous_pi(post_pred_samples, alpha)

# --- 2. Method 2: Marginal Intervals (Pointwise) ---
lower_marg <- apply(post_pred_samples, 2, quantile, probs = alpha / 2, na.rm = TRUE)
upper_marg <- apply(post_pred_samples, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
pi_marginal <- data.frame(lower = lower_marg, upper = upper_marg)

# --- 3. Method 3: SCSrank Simultaneous Interval ---
pi_scs_rank <- tryCatch({
  # Use conf.level = 1 - alpha
  conf_matrix <- SCSrank(post_pred_samples, conf.level = 1 - alpha)$conf.int
  # convert matrix to data.frame
  as.data.frame(conf_matrix)
}, error = function(e) {
  # Return an NA-filled data frame if SCSrank fails
  warning(paste("SCSrank failed on posterior samples:", e$message))
  data.frame(lower = rep(NA, ncol(post_pred_samples)),
              upper = rep(NA, ncol(post_pred_samples)))
})
# Ensure column names match the other data frames
colnames(pi_scs_rank) <- c("lower", "upper")


# Clean up Stan files
suppressWarnings(file.remove(fit$output_files()))

# Return all three intervals in a named list
return(list(
  mean_centered = pi_mean_centered,
  marginal = pi_marginal,
  scs_rank = pi_scs_rank
))


unique(df_tot_covprob_long$total_sim)

df_tot_covprob_long %>%
filter(total_sim > 300)


df_tot_eqt_long %>%
  filter(base_method == "bayesian_mcmc_cauchy_scs")

df_tot_eqt_long %>%
  filter(base_method == "bonferroni")

df_test <- df_asymmetry %>%
  filter(base_method == "bayesian_mcmc_cauchy_scs")
