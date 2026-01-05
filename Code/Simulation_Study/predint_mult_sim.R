# Setup -------------------------------------------------------------------

library(dplyr)
library(VGAM)
library(MCMCpack)
library(future.apply)
#library(MCPAN) # function is called in source
library(mvtnorm)
library(cmdstanr)
library(here)


# Run this in your R console
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# cmdstanr::install_cmdstan()

#set working directory
#setwd(".\\Code_and_Data")
setwd(here())
#setwd(r"(C:\Users\Budig\Google_Drive\Uni\Phd\03_Predint_multinomial\Code\Code_and_Data)")
# Functions
source(here("Code", "predint_mult_source.R"))

mod_gamma <- cmdstan_model(here("Code", "dirichlet_multinomial_gamma.stan"), compile = TRUE)
mod_cauchy <- cmdstan_model(here("Code", "dirichlet_multinomial_cauchy.stan"), compile = TRUE)
mod_beta <- cmdstan_model(here("Code", "dirichlet_multinomial_beta_rho.stan"), compile = TRUE)

# Settings ----------------------------------------------------------------

# phi: Dispersion parameter
# k: total number of clusters
# n: cluster size of historical data
# m: cluster size of future data

#Parameters for simulation
df_sim_settings_og <- expand.grid(
  k = c(5, 20, 100),
  n = c(10, 50, 100, 500),
  m = c(1),
  phi = c(1.01, 5, 8),
  B = 10000,
  quant_min = 0.01,
  quant_max = 10,
  alpha = 0.05,
  dispersion = "afroz",
  alternative = "both",
  tol = 1e-3,
  n_bisec = 60
)

#Parameters for simulation
# df_sim_settings_og <- expand.grid(
#   k = c(5, 20, 100),
#   n = c(10, 50, 100),
#   m = c(1),
#   phi = c(1.01, 5, 8),
#   B = 100,
#   quant_min = 0.01,
#   quant_max = 10,
#   alpha = 0.05,
#   dispersion = "afroz",
#   alternative = "both",
#   tol = 1e-3,
#   n_bisec = 60
# )

# Create a list of method flags
# l_methods <- list(
#   bisection = TRUE,
#   bisection_asym = TRUE,
#   Bonferroni = TRUE,
#   Bonf_bisec = TRUE,
#   pointwise = TRUE,
#   percentile_bt1 = TRUE,
#   mvn = TRUE,
#   percentile_bt_gemini = TRUE,
#   nonparametric_bt = TRUE, 
#   bayesian_mcmc = TRUE
# )

# Create a list of method flags
l_methods <- list(
  bisection = TRUE,
  bisection_asym = TRUE,
  Bonferroni = TRUE,
  Bonf_bisec = TRUE,
  pointwise = TRUE,
  percentile_bt1 = TRUE,
  mvn = TRUE,
  percentile_bt_gemini = TRUE,
  nonparametric_bt = FALSE,  
  bayesian_mcmc_gamma = FALSE, 
  bayesian_mcmc_cauchy = TRUE,
  bayesian_mcmc_beta = TRUE
)

# three categories
#l_props <- list(
#   matrix(c(0.33, 0.33, 0.33), ncol = 3, byrow = TRUE),
#   matrix(c(0.01, 0.01, 0.98), ncol = 3, byrow = TRUE),
#   matrix(c(0.25, 0.01, 0.74), ncol = 3, byrow = TRUE),
#   matrix(c(0.49, 0.02, 0.49), ncol = 3, byrow = TRUE),
#   matrix(c(0.25, 0.25, 0.5), ncol = 3, byrow = TRUE),
#   matrix(c(0.1, 0.3, 0.6), ncol = 3, byrow = TRUE)
# )

# further three categories

l_props <- list(
  matrix(c(0.33, 0.33, 0.33), ncol = 3, byrow = TRUE),
  matrix(c(0.01, 0.01, 0.98), ncol = 3, byrow = TRUE),
  matrix(c(0.25, 0.01, 0.74), ncol = 3, byrow = TRUE),
  matrix(c(0.49, 0.02, 0.49), ncol = 3, byrow = TRUE),
  matrix(c(0.25, 0.25, 0.5), ncol = 3, byrow = TRUE),
  matrix(c(0.1, 0.3, 0.6), ncol = 3, byrow = TRUE),
  matrix(c(0.02, 0.03, 0.95), ncol = 3, byrow = TRUE),
  matrix(c(0.05, 0.05, 0.9), ncol = 3, byrow = TRUE),
  matrix(c(0.05, 0.1, 0.85), ncol = 3, byrow = TRUE),
  matrix(c(0.05, 0.15, 0.8), ncol = 3, byrow = TRUE),
  matrix(c(0.1, 0.2, 0.7), ncol = 3, byrow = TRUE),
  matrix(c(0.05, 0.35, 0.65), ncol = 3, byrow = TRUE)

)

# # Five categories
# l_props <- list(
#   matrix(c(0.20, 0.20, 0.20, 0.20, 0.20), ncol = 5, byrow = TRUE),
#   matrix(c(0.30, 0.30, 0.20, 0.10, 0.10), ncol = 5, byrow = TRUE),
#   matrix(c(0.44, 0.22, 0.11, 0.11, 0.11), ncol = 5, byrow = TRUE),
#   matrix(c(0.50, 0.30, 0.10, 0.05, 0.05), ncol = 5, byrow = TRUE),
#   matrix(c(0.45, 0.27, 0.18, 0.08, 0.01), ncol = 5, byrow = TRUE),
#   matrix(c(0.70, 0.10, 0.10, 0.05, 0.05), ncol = 5, byrow = TRUE),
#   matrix(c(0.80, 0.10, 0.05, 0.04, 0.01), ncol = 5, byrow = TRUE),
#   matrix(c(0.10, 0.10, 0.20, 0.30, 0.30), ncol = 5, byrow = TRUE),
#   matrix(c(0.11, 0.11, 0.11, 0.22, 0.44), ncol = 5, byrow = TRUE),
#   matrix(c(0.05, 0.05, 0.10, 0.30, 0.50), ncol = 5, byrow = TRUE)
# )
# 
# # ten categories
# l_props <- list(
#   matrix(c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1), ncol = 10, byrow = TRUE),
#   matrix(c(0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2), ncol = 10, byrow = TRUE),
#   matrix(c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3), ncol = 10, byrow = TRUE),
#   matrix(c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.4), ncol = 10, byrow = TRUE),
#   matrix(c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.5), ncol = 10, byrow = TRUE),
#   matrix(c(0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 0.05, 0.05, 0.1, 0.6), ncol = 10, byrow = TRUE),
#   matrix(c(0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 0.05, 0.05, 0.35, 0.35), ncol = 10, byrow = TRUE),
#   matrix(c(0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.2, 0.2, 0.2), ncol = 10, byrow = TRUE),
#   matrix(c(0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 0.2, 0.2, 0.2, 0.2), ncol = 10, byrow = TRUE),
#   matrix(c(0.025, 0.025, 0.025, 0.025, 0.05, 0.05, 0.1, 0.2, 0.2, 0.3), ncol = 10, byrow = TRUE)
# )

# Split up into multiple Iterations to get intermediate results
#print(paste("Cores", availableCores()))


# if (availableCores() == 12){
#   sim_iterations <- 33
#   df_sim_settings_og$nsim <- 12
# } else if (availableCores() == 32){
#   sim_iterations <- 32
#   df_sim_settings_og$nsim <- 32
# } else {
#   sim_iterations <- 18
#   df_sim_settings_og$nsim <- 12
# }

# print(paste("Number of Iterations", sim_iterations))
# 
# lk <- length(unique(df_sim_settings_og$k))

# df_sim_settings_og$nsim <- 1 
# system.time(f_global_sim(l_props[[1]], df_sim_settings_og[1,], l_methods))

#27 settings 3 cats, 3334.16s 

# 
# cl <- parallelly::makeClusterPSOCK(
#   "biostat-32c",
#   user = "biostat",
#   rscript = "Rscript", 
#   verbose = TRUE,
#   homogeneous = FALSE
# )


# df_sim_settings_og$nsim <- 6
# 
# worker_hosts <- c("biostat-32a", "biostat-32b", "biostat-32c")
# cores_per_worker <- 2



# worker_hosts <- c("biostat-32a", "biostat-32b")
# cores_per_worker <- 32
# 
# df_sim_settings_og$nsim <- length(worker_hosts)*cores_per_worker
# 
# 
# 
# 
# all_cores_as_workers <- rep(worker_hosts, each = cores_per_worker)
# 
# cl <- parallelly::makeClusterPSOCK(
#   all_cores_as_workers,
#   user = rep("biostat", length(all_cores_as_workers)),
#   rscript = "Rscript", 
#   verbose = TRUE,
#   homogeneous = FALSE,
#   setup_strategy = "parallel" 
# )
# 
# 
# plan(cluster, workers = cl)


# parallelly::stopCluster(cl)
# parallel::stopCluster(cl)
# stopCluster(cl)
# 
# plan(sequential)

plan(multisession)

df_sim_settings_og$nsim <- availableCores()

l_props <- rotate_list(l_props, 1)

n_chunks <- 4

chunk_indices <- cut(seq_len(nrow(df_sim_settings_og)),
                     breaks = n_chunks,
                     labels = FALSE)

list_of_setting_chunks <- split(df_sim_settings_og, chunk_indices)

for (j in 1:unique(ceiling(1000 / df_sim_settings_og$nsim))) {
  # Now, loop through your probability vectors as before
  for (i in 1:length(l_props)) {
    print(Sys.time())
    print(paste("Processing prop set", i))
    
    system.time(for (j in 1:length(list_of_setting_chunks)) {
      # Get the current chunk of settings for this iteration
      current_settings_chunk <- list_of_setting_chunks[[j]]
      
      print(paste(
        "--- Starting Main Iteration for Chunk",
        j,
        "of",
        n_chunks,
        "---"
      ))
      print(paste(
        "This chunk contains",
        nrow(current_settings_chunk),
        "setting rows."
      ))
      
      print(current_settings_chunk)
      
      # Run the simulation function on the current chunk of settings
      system.time(f_global_sim(l_props[[i]], current_settings_chunk, l_methods))
    })
  }
}

# for (j in 1:30) {
#   print(paste("Main Iteration", j))  
# for (i in 1:(nrow(df_sim_settings_og)/lk)) {
#   df_sim_settings <- df_sim_settings_og[c((i*lk-(lk-1)):(i*lk)),]
#   system.time(for (x in 1:sim_iterations) {
#     print(paste("Iteration", x))
#     print(df_sim_settings)
#     system.time(for (i in 1:length(l_props)) {
#       f_global_sim(l_props[[i]], df_sim_settings, l_methods)
#     })
#   })
# }
# }

