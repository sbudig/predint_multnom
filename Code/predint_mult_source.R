# Functions ---------------------------------------------------------------

## expit function ---------------------------------------------------------

expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

## change lprop start------------------------------------------------------

# To change ordering of the original list of probability vectors
rotate_list <- function(input_list, start_index) {
  
  # Get the total number of elements in the list
  n <- length(input_list)

  # Check if the index is a single numeric value
  if (!is.numeric(start_index) || length(start_index) != 1) {
    stop("Error: 'start_index' must be a single number.")
  }
  
  # Check if the index is an integer
  if (start_index %% 1 != 0) {
    stop("Error: 'start_index' must be an integer.")
  }
  
  # Check if the index is within the valid range of the list
  if (start_index < 1 || start_index > n) {
    stop(paste("Error: 'start_index' must be between 1 and", n))
  }
  
  if (start_index == 1) {
    return(input_list)
  }
  
  # Create the new order of indices for the rotation
  new_order <- c( (start_index:n), (1:(start_index - 1)) )
  
  # Subset the list using the newly created rotational order
  rotated_list <- input_list[new_order]
  
  # Return the rotated list
  return(rotated_list)
}


## Sample from Dirichlet-multinomial distribution -------------------------

# k: a single integer; number of multinomial-dirichlet vectors to draw

# n: single integer or vector of integers, cluster size
#  the multinomial sample size for the individual multinomial vectors to generate
#  if vector: recycled to have length n

# v_pi_c: vector of probabilities, length (v_pi_c) must be = no. of categories
#  pi should sum to 1, pi is rescaled such that it really sums to 1 

# catnam: name of categories, should equal length v_pi vector

f_rdirmultinom <- function(v_pi, phi, k, n, catnam = NULL) {
  
  v_size <- rep(n, length.out = k)
  v_props <- v_pi / sum(v_pi)
  asum <- (phi - n) / (1 - phi)
  v_alpha <- asum * v_props
  m_propdir <- rdirichlet(n = k, alpha = v_alpha)
  
  m_dat <- t(apply(
    X = cbind(v_size, m_propdir),
    MARGIN = 1,
    FUN = function(x) {
      rmultinom(n = 1,
                size = x[1],
                prob = x[-1])
    }
  ))
  
  if (!is.null(catnam)) {
    colnames(m_dat) <- catnam
  }
  
  return(m_dat)
}

# count the number of zeroes per group and category
f_countzeroes <- function(df_dat){
  df_dat %>% group_by(level)  %>%
    summarise(across(everything(), mean))  %>% 
    summarise(zerocounts=sum(.==0))
}

# Fit a multinomial model function ------------------------------------------
# fits a multinomial or dirichlet-multinomial model to the given dataset
# for pearson, afroz, fletcher and deviance, also just a multinomial model
# is fitted, as these dispersion parameters are calculated afterwards,
# based on the multinomial model outcome
f_fit_mult <-
  function(df_i,
           modeltype = c("multinomial",
                         "pearson",
                         "afroz",
                         "fletcher",
                         "deviance",
                         "dirmultinomial")) {
    
    modeltype <- match.arg(modeltype)
    
    if ("level" %in% names(df_i)) {
      t_formula <- paste0("cbind(",
                          paste0(colnames(df_i[2:(length(df_i) - 1)]),
                                 collapse = ","),",V1) ~ 1")
    } else {
      t_formula <- paste0("cbind(",
                          paste0(colnames(df_i[2:(length(df_i))]),
                                 collapse = ","),",V1) ~ 1")  
    }
    
    
    warning_text <- NULL
    
    switch(
      modeltype,
      multinomial = ,
      pearson = ,
      afroz = ,
      fletcher = ,
      deviance = {
        modeloutput <-tryCatch( withCallingHandlers(
          {
            vglm(
              formula = t_formula,
              family = multinomial,
              model = TRUE,
              data = df_i
            )},
          warning = function(w){
            warning_text <<- c(warning_text, conditionMessage(w))
            invokeRestart("muffleWarning")
          }),
          error = function(e) {
            warning_text <<-  trimws(e)
            modeloutput <- NULL
          }
        )
        
        
      },
      dirmultinomial = {
        modeloutput <-tryCatch(withCallingHandlers({
          vglm(
            formula = t_formula,
            family = dirmultinomial,
            model = TRUE,
            data = df_i,
            maxit = 100
          )},
          warning = function(w){
            warning_text <<- c(warning_text, conditionMessage(w))
            invokeRestart("muffleWarning")
          }),
          error = function(e) {
            warning_text <<-   trimws(e)
            modeloutput <- NULL
          }
        )
        
        
      }
    )
    
    return(list(modeloutput,unique(warning_text)))
  }

# check if there is any error or warning in model fitting process ---------

f_check_err_warn <- function(l_errwarn){
  
  # "NAs found in the working weights variable 'wz'"
  checkerrnas <- unlist(lapply(l_errwarn, function(i) {
    if(any(is.null(i))){0
    } else if (any(grepl("NAs", i, fixed = TRUE))) {1} else (0)
  }))
  
  # "Some elements in the working weights variable 'wz' are not finite"
  checkerrfin <- unlist(lapply(l_errwarn, function(i) {
    if(any(is.null(i))){0
    } else if (any(grepl("finite", i, fixed = TRUE))) {1} else (0)
  }))
  
  # check if there is a working wz warning
  # "diagonal elements of the working weights variable 'wz' have been replaced by 1.819e-12"
  checkwarwz <- unlist(lapply(l_errwarn, function(i) {
    if(any(is.null(i))){0
    } else if (
      any(grepl("diagonal", i, fixed = TRUE))) {1} else (0)
  }))
  
  # check if there is a convergence warning
  # "convergence not obtained in 30 IRLS iterations"
  checkwarcon <- unlist(lapply(l_errwarn, function(i) {
    if(any(is.null(i))){0
    } else if (
      any(grepl("convergence", i, fixed = TRUE))) {1} else {0}
  }))
  
  # check if there are any deleted columns due to zero counts
  # "Deleted 2 columns of the response matrix due to zero counts"
  checkwarzc <- unlist(lapply(l_errwarn, function(i) {
    if(any(is.null(i))){0
    } else if (
      any(grepl("Deleted", i, fixed = TRUE))) {1} else (0)
  }))
  
  # check if there is any other warning
  checkwarany <- unlist(lapply(l_errwarn, function(i) {
    if(any(is.null(i))){0
    } else if (
      any(grepl("diagonal", i, fixed = TRUE)) | 
      any(grepl("convergence", i, fixed = TRUE)) |
      any(grepl("Error", i, fixed = TRUE)) |
      any(grepl("Deleted", i, fixed = TRUE))) {0} else (1)
  }))
  
  return(data.frame(cerrnas = checkerrnas,
                    cerrfin = checkerrfin,
                    cwarwz = checkwarwz,
                    cwarcon = checkwarcon,
                    cwarzc = checkwarzc,
                    cwarany = checkwarany
  ))
  
}

# Estimate Dispersion Parameter ------------------------------------------
# estimates a dispersion parameter based on the model. 
# multinomial means no dispersion is estimated and its set to 1
# pearson, afroz, farrington and deviance estimators can be calculated
# when the underlying model is a multinomial model
# dirmultinomial gives the dispersion parameter, when Dirichlet-multinomial
# model was used on the data
f_estimate_phi <-
  function(model,
           dispersion = c("multinomial",
                          "pearson",
                          "afroz",
                          "farrington",
                          "deviance",
                          "dirmultinomial")) {
    disptype <- match.arg(dispersion)
    
    if (is.null(model)) {
      return(NA)
    }
    
    if (disptype == "dirmultinomial") {
      icc <-
        expit(model@coefficients[length(model@misc$predictors.names)][[1]])
      phi <- 1 + (unique(model@extra$n2) - 1) * icc
      return(phi)
    }
    
    y_vglm <- model@y
    y_vglm <- matrix(y_vglm, ncol = ncol(y_vglm))
    N_vglm <- length(y_vglm)
    w_vglm <- as.vector(model@prior.weights)
    
    pi_vglm <- model@fitted.values
    pi_vglm <- as.matrix(pi_vglm, ncol = ncol(pi_vglm))
    
    n_vglm <- model@misc$n
    
    s.bar_vglm <- sum((y_vglm - pi_vglm)/pi_vglm/(N_vglm - n_vglm))
    
    X2_vglm <- sum((y_vglm*w_vglm - pi_vglm*w_vglm) ^ 2 / (pi_vglm*w_vglm))
    
    phi.Pearson <- X2_vglm / model@df.residual
    
    # Switch-case for dispersion type
    phi <- switch(
      disptype,
      multinomial = 1,
      pearson = phi.Pearson,
      afroz = phi.Pearson / (1 + s.bar_vglm),
      farrington = phi.Pearson - (N_vglm - n_vglm) * s.bar_vglm / model@df.residual,
      deviance = model@criterion$deviance / model@df.residual
    )
    return(phi)
  }

# Function to calculate probabilities for multinomial model ---------------
# calculates the probabilities for each category c
f_calc_pi <- function(mult_model) {
  # Extract the intercept coefficients
  coefs <- unname(coef(mult_model))
  
  # Compute the exponential of the intercepts
  exp_coefs <- exp(coefs)
  
  # Calculate probabilities using the softmax function
  probs <- exp_coefs / (1 + sum(exp_coefs))
  
  # Add the probability for the first category 
  prob_first <- 1 - sum(probs)
  
  # Combine the probabilities for all categories
  return(c(prob_first, probs))
}

# Calculate the expected counts with pi and cluster size ------------------
# takes probability pi and cluster size m to calculate the expected counts
f_calc_y_hat <- function(pi_hat, m){
  y_hat <- pi_hat * m
  return(y_hat)
}

# Calculate variance of future observation --------------------------------
f_calc_var_y <- function(phi_hat, m, pi_hat){
  var_y <- phi_hat*m*pi_hat*(1-pi_hat)
  return(var_y)
}

# Calculate variance of the estimated future observation ------------------
f_calc_var_y_hat <- function(phi_hat, m, pi_hat, N_hist){
  var_y_hat <- phi_hat*m^2*pi_hat*(1-pi_hat)/N_hist
  return(var_y_hat)
}

# Calculate standard error for prediciton ---------------------------------
f_calc_se_pred <- function(var_y, var_y_hat){
  se_pred <- sqrt(var_y+var_y_hat)
  return(se_pred)
}

# Replace estimated phis that are below 1 ---------------------------------
f_replace_phi_below_one <- function(phi) {
  phi <- lapply(phi, function(x) {
   if (is.numeric(x)){
    if (x <= 1) {
      x <- 1.01
      return(x)
    } else {
      return(x)
    }
   } else {
     return(x)
   }
  })
  return(phi)
}

# Replace estimated phis that are above a specific value-------------------
f_replace_phi_above_val <- function(phi, val) {
  phi <- lapply(phi, function(x) {
    if (is.numeric(x)){
    if (x >= val) {
      return(val*0.975)
    } else {
      return(x)
    }
    } else {
      return(x)
      }
  })
  return(phi)
}


# Internal function to calculate standard error of the historical  --------

f_calc_se_pred_star_internal <- function(phi_hat, pi_hat_vec, m, N_hist){
  
  # variance of fut. random variable
  var_y_star <- m * phi_hat * pi_hat_vec * (1-pi_hat_vec)
  
  # variance of fut. expectation
  var_y_star_hat <- (m^2 * phi_hat * pi_hat_vec * (1-pi_hat_vec))/N_hist
  
  # prediction SE
  se_pred <- sqrt(var_y_star + var_y_star_hat)
  
  
  return(se_pred)
}

# Calculate standard error of the historical bootstraps -------------------

f_calc_se_pred_star <- function(l_phi_hat, l_pi_hat_vec, m, N_hist){
  
  l_se_pred_star <- mapply(FUN = f_calc_se_pred_star_internal,
                           phi_hat = l_phi_hat,
                           pi_hat_vec = l_pi_hat_vec,
                           MoreArgs = list(m = m, 
                                           N_hist = N_hist),
                           SIMPLIFY = FALSE)
  
  return(l_se_pred_star)
}

# Calculate the prediciton intervals --------------------------------------

f_calc_prediction_interval <- function(y_hat, q, se_pred, alternative)
{
  
  if (is.vector(q)){
    
    if (alternative == "both") {
      if (length(q) == 1) {
        lower <- y_hat - q * se_pred
        upper <- y_hat + q * se_pred
        
        out <- data.frame(lower, upper)
        return(out)
      }
      
      if (length(q) == 2) {
        lower <- y_hat - q[1] * se_pred
        upper <- y_hat + q[2] * se_pred
        
        out <- data.frame(lower, upper)
        return(out)
      }
      
    }
    
    if (alternative == "lower") {
      lower <- y_hat - q * se_pred
      
      out <- data.frame(lower)
      return(out)
    }
    
    if (alternative == "upper") {
      upper <- y_hat + q * se_pred
      
      out <- data.frame(upper)
      return(out)
    }
  } else if (is.matrix(q)) {
    
    out <- data.frame()
    
    for (i in 1:nrow(q)) {
      
      lower <- y_hat[i] - q[i,1] * se_pred[i]
      upper <- y_hat[i] + q[i,2] * se_pred[i]
      
      out <- rbind(out, data.frame(lower, upper))
      
    }
    return(out)
    
  } else {
    return(NULL)
  }
  
}

# Extract whether future obs are higher, lower or within the interv -------

f_compare_with_intervals <- function(data, intervals) {
  result <- sapply(names(data)[-which(names(data) == "level")], function(col) {
    value <- data[[col]]
    lower <- intervals$lower[which(names(data) == col)]
    upper <- intervals$upper[which(names(data) == col)]
    
    if (value < lower) {
      -1
    } else if (value > upper) {
      1
    } else {
      0
    }
  })
  return(result)
}

# Calculate adjusted q ----------------------------------------------------
# adjust q depending on number of n (Bonferroni adjustment)
f_adjust_q <- function(n, alpha) {
  qnorm(1 - (alpha/n)/ 2)
}

# Function to assess the coverage probability -----------------------------
f_coverage_prob <- function(y_star_hat,
                            se_pred,
                            q,
                            y_star,
                            alternative){
  
  # y_star_hat has to be a list
  if(!is.list(y_star_hat)){
    stop("!is.list(y_star_hat)")
  }
  
  # se_pred has to be a list
  if(!is.list(se_pred)){
    stop("!is.list(se_pred)")
  }
  
  # y_star has to be a list
  if(!is.list(y_star)){
    stop("!is.list(y_star)")
  }
  
  # all three lists have to have the same length
  if(length(unique(c(length(y_star_hat), length(se_pred), length(y_star)))) != 1){
    stop("length(unique(c(length(y_star_hat), length(se_pred), length(y_star)))) != 1")
  }
  
  # q needs to be numeric
  if(!is.numeric(q)){
    stop("!is.numeric(q)")
  }
  
  # alternative must be defined
  if(isTRUE(alternative!="both" && alternative!="lower" && alternative!="upper")){
    stop("alternative must be either both, lower or upper")
  }
  
  
  # q times se_pred
  q_se_list <- mapply(FUN = function(x,y){x*y},
                      x = se_pred,
                      MoreArgs = list(y=q),
                      SIMPLIFY=FALSE)
  
  if(alternative=="both"){
    
    # Lower border
    lower_list <- mapply(FUN = function(x,y){x-y},
                         x = y_star_hat,
                         y = q_se_list,
                         SIMPLIFY=FALSE)
    
    # Upper border
    upper_list <- mapply(FUN = function(x,y){x+y},
                         x = y_star_hat,
                         y = q_se_list,
                         SIMPLIFY=FALSE)
    
  }
  
  if(alternative=="lower"){
    
    # Lower border
    lower_list <- mapply(FUN = function(x,y){x-y},
                         x = y_star_hat,
                         y = q_se_list,
                         SIMPLIFY=FALSE)
  }
  
  if(alternative=="upper"){
    
    # Upper border
    upper_list <- mapply(FUN = function(x,y){x+y},
                         x = y_star_hat,
                         y = q_se_list,
                         SIMPLIFY=FALSE)
    
  }
  
  
  # Function to check the coverage
  cover_fun <- function(lower=NULL,
                        upper=NULL,
                        y_star,
                        alternative){
    
    if(alternative=="both"){
      
      # If  all y_star are covered set output to 1
      if(all(lower <  y_star & y_star < upper)){
        return(1)
      }
      
      # If not all y_star are covered set output to 0
      if(!all(lower <  y_star & y_star < upper)){
        return(0)
      }
    }
    
    if(alternative=="lower"){
      
      # If  all y_star are covered set output to 1
      if(all(lower <  y_star)){
        return(1)
      }
      
      # If not all y_star are covered set output to 0
      if(!all(lower <  y_star )){
        return(0)
      }
    }
    
    if(alternative=="upper"){
      
      # If  all y_star are covered set output to 1
      if(all(y_star < upper)){
        return(1)
      }
      
      # If not all y_star are covered set output to 0
      if(!all(y_star < upper)){
        return(0)
      }
    }
  }
  
  
  if(alternative=="both"){
    
    cover_list <- mapply(FUN = cover_fun,
                         lower=lower_list,
                         upper=upper_list,
                         y_star=y_star,
                         MoreArgs = list(alternative="both"),
                         SIMPLIFY=FALSE)
  }
  
  if(alternative=="lower"){
    
    cover_list <- mapply(FUN = cover_fun,
                         lower=lower_list,
                         y_star=y_star,
                         MoreArgs = list(alternative="lower"),
                         SIMPLIFY=FALSE)
  }
  
  if(alternative=="upper"){
    
    cover_list <- mapply(FUN = cover_fun,
                         upper=upper_list,
                         y_star=y_star,
                         MoreArgs = list(alternative="upper"),
                         SIMPLIFY=FALSE)
  }
  
  
  cover_prob <- sum(unlist(cover_list))/length(cover_list)
  
  return(cover_prob)
  
}


# Bisection function ------------------------------------------------------
f_bisection <- function(y_star_hat,
                        se_pred,
                        y_star,
                        alternative,
                        quant_min,
                        quant_max,
                        n_bisec,
                        tol,
                        alpha,
                        traceplot=TRUE){
  
  # y_star_hat has to be a list
  if(!is.list(y_star_hat)){
    stop("!is.list(y_star_hat)")
  }
  
  # se_pred has to be a list
  if(!is.list(se_pred)){
    stop("!is.list(se_pred)")
  }
  
  # y_star has to be a list
  if(!is.list(y_star)){
    stop("!is.list(y_star)")
  }
  
  # all three lists have to have the same length
  if(length(unique(c(length(y_star_hat), length(se_pred), length(y_star)))) != 1){
    stop("length(unique(c(length(y_star_hat), length(se_pred), length(y_star)))) != 1")
  }
  
  # all elements of y_star_hat have to be of length M
  if(length(unique(sapply(y_star_hat, length))) != 1){
    stop("all elements in y_star_hat have to have the same length M")
  }
  
  # all elements of se_pred have to be of length M
  if(length(unique(sapply(se_pred, length))) != 1){
    stop("all elements in se_pred have to have the same length M")
  }
  
  # all elements of y_star have to be of length M
  if(length(unique(sapply(y_star, length))) !=1){
    stop("all elements in y_star have to have the same length M")
  }
  
  # # M has to be the same for all three lists
  # if(length(unique(unique(sapply(y_star_hat, length)), unique(sapply(se_pred, length)), unique(sapply(y_star, length)))) != 1){
  # # if(var(c(unique(sapply(y_star_hat, length)), unique(sapply(se_pred, length)), unique(sapply(y_star, length)))) != 0){
  #         stop("M differs between y_star_hat, se_pred and y_star")
  # }
  
  
  # alternative must be defined
  if(isTRUE(alternative!="both" && alternative!="lower" && alternative!="upper")){
    stop("alternative must be either both, lower or upper")
  }
  
  # quant_min needs to be numeric
  if(!is.numeric(quant_min)){
    stop("!is.numeric(quant_min)")
  }
  
  # quant_max needs to be numeric
  if(!is.numeric(quant_max)){
    stop("!is.numeric(quant_max)")
  }
  
  # n_bisec needs ot be an integer number
  if(!is.numeric(n_bisec)){
    stop("!is.numeric(n_bisec)")
  }
  
  if(!all(n_bisec == floor(n_bisec))){
    stop("!all(n_bisec == floor(n_bisec))")
  }
  
  # tolerance needs to bve a number
  if(!is.numeric(tol)){
    stop("!is.numeric(tol)")
  }
  
  # Tolerance needs to be bigger than 0
  if(!(tol > 0)){
    stop("tol needs to be bigger than 0")
  }
  
  # Tolerance needs to be small to yield accurate results
  if(tol>0.01){
    warning("The tolerance is higher than 0.01: The bisection resulds might be imprecise.")
  }
  
  if(!isTRUE(traceplot)){
    if(!(traceplot==FALSE)){
      stop("traceplote needs to be TRUE or FALSE")
    }
  }
  
  
  c_i <- vector()
  runval_i <- vector()
  
  
  # Calculate coverages for start points
  
  cover_quant_min <- f_coverage_prob(y_star_hat = y_star_hat,
                                     se_pred = se_pred,
                                     q = quant_min,
                                     y_star = y_star,
                                     alternative = alternative)
  
  cover_quant_max <- f_coverage_prob(y_star_hat = y_star_hat,
                                     se_pred = se_pred,
                                     q = quant_max,
                                     y_star = y_star,
                                     alternative = alternative)
  
  
  # if the coverage is bigger for both quant take quant_min
  if ((cover_quant_min > 1-alpha+tol)) {
    
    warning(paste("observed coverage probability of ",
                  cover_quant_min,
                  "for quant_min is bigger than 1-alpha+tol =",
                  1-alpha+tol))
    
    if(traceplot==TRUE){
      
      plot(x=quant_min,
           y=cover_quant_min-(1-alpha),
           type="p",
           pch=20,
           xlab="calibration value",
           ylab="obs. coverage - nom. coverage",
           main=paste("f(quant_min) > 1-alpha+tol"),
           ylim=c(cover_quant_min-(1-alpha)+tol, -tol))
      abline(a=0, b=0, lty="dashed")
      abline(a=tol, b=0, col="grey")
    }
    
    return(quant_min)
  }
  
  
  # if the coverage is bigger for both quant take quant_max
  
  else if ((cover_quant_max < 1-alpha-tol)) {
    
    warning(paste("observed coverage probability of",
                  cover_quant_max,
                  "for quant_max is smaller than 1-alpha-tol =",
                  1-alpha-tol))
    
    
    if(traceplot==TRUE){
      
      plot(x=quant_max,
           y=cover_quant_max-(1-alpha),
           type="p", pch=20,
           xlab="calibration value",
           ylab="obs. coverage - nom. coverage",
           main=paste("f(quant_max) < 1-alpha-tol"),
           ylim=c(cover_quant_max-(1-alpha)-tol, tol))
      abline(a=0, b=0, lty="dashed")
      abline(a=-tol, b=0, col="grey")
    }
    
    
    return(quant_max)
  }
  
  
  # run bisection
  
  else for (i in 1:n_bisec) {
    
    # Calculate midpoint
    c <- (quant_min + quant_max) / 2
    
    cprop <- f_coverage_prob(y_star_hat = y_star_hat,
                             se_pred = se_pred,
                             q = c,
                             y_star = y_star,
                             alternative = alternative)
    
    runval <- (1-alpha) - cprop
    
    # Assigning c and runval into the vectors
    c_i[i] <- c
    runval_i[i] <- runval
    # runval_i[i] <- cprop
    
    # if ((1-alpha)-tol < cprop & cprop < (1-alpha)+tol) {
    if (abs(runval) < tol) {
      
      if(traceplot==TRUE){
        
        plot(x=c_i,
             y=-runval_i,
             type="p",
             pch=20,
             xlab="calibration value",
             ylab="obs. coverage - nom. coverage",
             main=paste("Trace with", i, "iterations"))
        lines(x=c_i, y=-runval_i, type="s", col="red")
        abline(a=0, b=0, lty="dashed")
        abline(a=tol, b=0, col="grey")
        abline(a=-tol, b=0, col="grey")
      }
      
      return(c)
    }
    
    if(sign(runval)==1){
      quant_min <- c}
    
    else if(sign(runval)==-1){
      quant_max <- c}
    
    
  }
  warning('Too many iterations, but the quantile of the last step is returned')
  
  if(traceplot==TRUE){
    
    plot(x=c_i,
         y=runval_i,
         type="p",
         pch=20,
         xlab="calibration value",
         ylab="obs. coverage - nom. coverage",
         main=paste("Trace with", i, "iterations"))
    lines(x=c_i, y=runval_i, type="s", col="red")
    abline(a=0, b=0, lty="dashed")
    abline(a=tol, b=0, col="grey")
    abline(a=-tol, b=0, col="grey")
  }
  return(c)
  
} 

# bisection head function for bonferroni bisection calibration ------------

f_single_bisec <- function(n_cat,
                           y_star_hat,
                           se_pred,
                           y_star,
                           quant_min,
                           quant_max,
                           n_bisec,
                           tol,
                           alpha,
                           alternative,
                           traceplot){
  
  m_quant_calib <- c()
  
  for (i in 1:n_cat) {
    y_star_hat_t <- lapply(y_star_hat, `[`, i)
    se_pred_t <- lapply(se_pred, `[`, i)
    y_star_t <- lapply(y_star, `[`, i)
    
    quant_calib_lower_t <- f_bisection(y_star_hat = y_star_hat_t,
                                       se_pred = se_pred_t,
                                       y_star = y_star_t,
                                       quant_min = quant_min,
                                       quant_max = quant_max,
                                       n_bisec = n_bisec,
                                       tol = tol,
                                       alpha = alpha/(n_cat*2),
                                       alternative = "lower",
                                       traceplot = FALSE
    )
    
    quant_calib_upper_t <- f_bisection(y_star_hat = y_star_hat_t,
                                       se_pred = se_pred_t,
                                       y_star = y_star_t,
                                       quant_min = quant_min,
                                       quant_max = quant_max,
                                       n_bisec = n_bisec,
                                       tol = tol,
                                       alpha = alpha/(n_cat*2),
                                       alternative = "upper",
                                       traceplot = FALSE
    )
    
    m_quant_calib <- rbind(m_quant_calib, c(quant_calib_lower_t, quant_calib_upper_t))
  }
  
  return(m_quant_calib)
}

# SCSrank function form package MCPAN -------------------------------------

SCSrank <-
  function(x, conf.level=0.95, alternative="two.sided")
  {
    alternative <- match.arg(alternative, choices=c("two.sided","less","greater"))
    
    DataMatrix <- x
    N <- nrow(DataMatrix)
    k <- round(conf.level*N,0)
    RankDat <- apply(DataMatrix,2,rank)
    
    switch(alternative,
           
           "two.sided"={
             W1 <- apply(RankDat,1,max)
             W2 <- N + 1 - apply(RankDat,1,min)
             
             Wmat <- cbind(W1,W2)
             w <- apply(Wmat,1,max)
             tstar <- round(sort(w)[k],0)
             
             SCI <- function(x)
             {
               sortx <- sort(x)
               cbind(sortx[N+1-tstar],sortx[tstar])
             }
             
             SCS <- t(apply(DataMatrix,2,SCI))
           },
           
           "less"={
             W1 <- apply(RankDat,1,max)
             tstar <- round(sort(W1)[k],0)
             
             SCI <- function(x)
             {
               sortx <- sort(x)
               cbind(-Inf, sortx[tstar])
             }
             
             SCS<-t(apply(DataMatrix,2,SCI))
           },
           
           "greater"={
             W2 <- N + 1 - apply(RankDat,1,min)
             tstar <- round(sort(W2)[k],0)
             
             SCI <- function(x)
             {
               sortx <- sort(x)
               cbind(sortx[N+1-tstar], Inf)
             }
             
             SCS<-t(apply(DataMatrix,2,SCI))
             
           }
    )
    # end of switch
    
    colnames(SCS)<-c("lower","upper")
    
    attr(SCS, which="k")<-k
    attr(SCS, which="N")<-N
    OUT<-list(conf.int=SCS, conf.level=conf.level, alternative=alternative)
    return(OUT)
  }

# Calculate critical value based on MVN approximation --------------------
f_calc_q_mvn <- function(phi_hat, pi_hat, m, N_hist, alpha) {
  
  if (is.null(phi_hat) || is.null(pi_hat)) {
    return(NA)
  }
  
  C <- length(pi_hat)
  
  if (C <= 1) {
    return(qnorm(1 - alpha / 2))
  }
  
  q_mvn <- tryCatch({
    
    Sigma_pi_hat <- diag(pi_hat) - pi_hat %*% t(pi_hat)
    Sigma_pred <- phi_hat * m * (1 + m / N_hist) * Sigma_pi_hat
    
    R_pred <- cov2cor(Sigma_pred)
    
    qmvnorm(1 - alpha, corr = R_pred, tail = "both.tails")$quantile
    
  }, error = function(e) {
    return(NA)
  })
  
  return(q_mvn)
}

f_calc_zc_help <- function(y_star, y_star_hat, se){
  zc <- (y_star-y_star_hat)/se
  return(zc)
}

f_calc_zc_main <- function(l_y_star, l_y_star_hat, l_se, n_cat) {
  
  if (length(l_y_star) == 0 || length(l_y_star_hat) == 0 || length(l_se) == 0) {
    return(NULL)
  }
  
  valid_indices <- which(
    sapply(l_y_star, function(x) !is.null(x) && length(x) == n_cat) &
      sapply(l_y_star_hat, function(x) !is.null(x) && length(x) == n_cat) &
      sapply(l_se, function(x) !is.null(x) && length(x) == n_cat)
  )
  
  if (length(valid_indices) == 0) {
    return(NULL)
  }
  
  l_y_star_valid <- l_y_star[valid_indices]
  l_y_star_hat_valid <- l_y_star_hat[valid_indices]
  l_se_valid <- l_se[valid_indices]
  
  zc_list <- mapply(
    FUN = f_calc_zc_help,
    y_star = l_y_star_valid,
    y_star_hat = l_y_star_hat_valid,
    se = l_se_valid,
    SIMPLIFY = FALSE
  )

  do.call(rbind, zc_list)
  
}

# Non-parametric Bootstrap data generation --------------------------------
# Resamples the original historical data clusters with replacement
f_bt_mult_dat_nonparametric <- function(df_hist, B) {
  
  k <- nrow(df_hist)
  
  replicate(n = B, {
    resampled_indices <- sample(1:k, size = k, replace = TRUE)
    df_hist[resampled_indices, ]
  }, simplify = FALSE)
}

# Adjust a list of phi values based on a minimum and maximum value --------
f_adjust_phi_list <- function(l_phi, m_value) {
  
  if (is.null(l_phi)) {
    return(NULL)
  }
  
  lapply(l_phi, function(phi) {
    
    if (is.null(phi) || is.na(phi)) {
      return(phi)
    }

    if (phi <= 1) {
      return(1.01)
    } else if (phi > m_value) {
      return(m_value * 0.975)
    } else {
      return(phi)
    }
  })
}

#  Helper function to construct a simultaneous PI from posterior p --------
f_calc_pi_bayesian_standardized <- function(post_pred_samples, alpha) {
  
  # 1. Calculate moments of the posterior predictive distribution
  y_bar_bayes <- colMeans(post_pred_samples)
  y_sd_bayes  <- apply(post_pred_samples, 2, sd)
  
  # Handle zero variance if a category never occurs (prevent div by zero)
  y_sd_bayes[y_sd_bayes == 0] <- 1e-6 
  
  # 2. Calculate standardized residuals for every sample
  # (Sample - Mean) / SD
  # sweep subtracts mean, then sweep divides by SD
  z_scores <- sweep(post_pred_samples, 2, y_bar_bayes, "-")
  z_scores <- sweep(z_scores, 2, y_sd_bayes, "/")
  
  # 3. Find the maximum absolute standardized deviation for each sample (Simultaneous step)
  max_abs_z <- apply(abs(z_scores), 1, max)
  
  # 4. Find the critical value (1-alpha quantile of the MAX z-scores)
  q_crit_standardized <- quantile(max_abs_z, probs = 1 - alpha, na.rm = TRUE)
  
  # 5. Construct Interval: Mean +/- q * SD
  lower <- y_bar_bayes - q_crit_standardized * y_sd_bayes
  upper <- y_bar_bayes + q_crit_standardized * y_sd_bayes
  
  # Floor/Ceiling to ensure integer bounds (optional but recommended for counts)
  #lower <- floor(lower)
  #upper <- ceiling(upper)
  
  # Ensure non-negative lower bound
  lower[lower < 0] <- 0
  
  return(data.frame(lower, upper))
}

# Calculate Bayesian PI using a Hierarchical MCMC model -------------------
#
# This function now calculates THREE types of intervals from a single MCMC run:
# 1. mean_centered: The original simultaneous interval based on max deviation from the mean.
# 2. marginal: Pointwise (non-simultaneous) intervals using simple quantiles.
# 3. scs_rank: A simultaneous interval using the SCSrank method on the posterior samples.
#
# @return A list containing three data frames: $mean_centered, $marginal, $scs_rank

f_calc_pi_bayesian_mcmc <- function(df_hist, m, alpha, stan_model) {
  
  # Ensure cmdstanr is loaded
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("Package 'cmdstanr' is required for the MCMC Bayesian method.")
  }
  
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
      refresh = 0,
      show_messages = FALSE
    )
  )
  
  # Extract posterior predictive samples
  post_pred_samples <- fit$draws("y_pred", format = "matrix")
  
  # --- 1. Method 1: Standardized (Pivotal) Simultaneous Interval  ---
  pi_mean_centered <- f_calc_pi_bayesian_standardized(post_pred_samples, alpha)
  
  # --- 2. Method 2: Marginal Intervals (Pointwise) ---
  C <- ncol(post_pred_samples)
  alpha_bonf <- alpha / C
  lower_marg <- apply(post_pred_samples, 2, quantile, probs = alpha_bonf / 2, na.rm = TRUE)
  upper_marg <- apply(post_pred_samples, 2, quantile, probs = 1 - alpha_bonf / 2, na.rm = TRUE)
  pi_marginal <- data.frame(lower = lower_marg, upper = upper_marg)
  
  # --- 3. Method 3: SCSrank Simultaneous Interval ---
  pi_scs_rank <- tryCatch({
    conf_matrix <- SCSrank(post_pred_samples, conf.level = 1 - alpha)$conf.int
    as.data.frame(conf_matrix)
  }, error = function(e) {
    warning(paste("SCSrank failed on posterior samples:", e$message))
    data.frame(lower = rep(NA, ncol(post_pred_samples)),
               upper = rep(NA, ncol(post_pred_samples)))
  })

  colnames(pi_scs_rank) <- c("lower", "upper")
  
  suppressWarnings(file.remove(fit$output_files()))
  
  return(list(
    mean_centered = pi_mean_centered,
    marginal = pi_marginal,
    scs_rank = pi_scs_rank
  ))
}

# Calculate statistics for interval widths --------------------------------
#
# @param df_interval A data frame with 'lower' and 'upper' columns.
# @return A single-row data frame with average, median, sd, and quantile widths.

f_calc_interval_width_stats <- function(df_interval) {
  
  # Return NA for all stats if the input is NULL (e.g., if a method failed)
  if (is.null(df_interval) || nrow(df_interval) == 0) {
    return(data.frame(avg_width = NA,
                      median_width = NA,
                      sd_width = NA,
                      q25_width = NA,
                      q75_width = NA))
  }
  
  # Calculate the width for each category
  widths <- df_interval$upper - df_interval$lower
  
  # Calculate summary statistics, handling potential NA values
  avg_width    <- mean(widths, na.rm = TRUE)
  median_width <- median(widths, na.rm = TRUE)
  sd_width     <- sd(widths, na.rm = TRUE)
  quantiles    <- quantile(widths, probs = c(0.25, 0.75), na.rm = TRUE)
  q25_width    <- quantiles[1]
  q75_width    <- quantiles[2]
  
  # Return as a single-row data frame
  return(data.frame(avg_width, median_width, sd_width, q25_width, q75_width))
}


# Simulation function -----------------------------------------------------
# Main simulation function

f_predint_mult_sim <- function(m_true_pi_vec, df_sim_settings, l_methods,  mod_gamma, mod_cauchy, mod_beta){
  
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
  
  # Remove l_mult_dat to free up memory
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
          # This block is executed only if an error occurs in the code above.
          message("Main bootstrap x: An error in data generation was caught. Retrying...")
          count_error_mbt <- count_error_mbt + 1
          
          if (count_error_mbt > 10) {
            message("Main bootstrap x: Ten errors occurred and were caught. Retrying with phi*0.99")
            current_phi_n <- current_phi_n*0.99
          }
        })
      }
      
      # Handle all-zero columns
      for (col in colnames(x_star_b)) {
        if (all(x_star_b[[col]] == 0)) {
          x_star_b[[sample(nrow(x_star_b), 1), col]] <- 1
        }
      }
      
      # Fit model to the bootstrap sample
      model_x_star <- f_fit_mult(df_i = x_star_b, modeltype = "multinomial")
      
      if (!is.null(model_x_star[[1]])) {
        # Calculate bootstrap statistics
        pi_star_hat_b <- f_calc_pi(model_x_star[[1]])
        phi_star_hat_b <- f_estimate_phi(model = model_x_star[[1]], dispersion = current_disp)
        
        pi_star_hat_b_list[[b]] <- pi_star_hat_b
        phi_star_hat_b_list[[b]] <- phi_star_hat_b
        
        y_star_hat_b_list[[b]] <- f_calc_y_hat(pi_star_hat_b, current_m)
        se_pred_star_b_list[[b]] <- f_calc_se_pred_star_internal(phi_star_hat_b, pi_star_hat_b, current_m, current_N_hist)
        
        # Generate ONE future bootstrap sample (if needed by methods)
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
              message("Main bootstrap y: An error in data generation was caught. Retrying...")
              count_error_mbty <- count_error_mbty + 1
              
              if (count_error_mbty > 10) {
                message("Main bootstrap y: Ten errors occurred and were caught. Retrying with phi*0.99")
                current_phi_m <- current_phi_m*0.99
              }
            })
          }
          
          y_star_b_list[[b]] <- y_star_b
        }
      }
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
  
  
  # use Percentile Bootstrap (Gemini Method) to compute prediction intervals
  if (l_methods$percentile_bt_gemini) {
    
    # This part remains the same: calculate the matrix of pivotal quantities (standardized residuals)
    # where each row is a bootstrap replication and each column is a category.
    l_df_zc <- mapply(FUN = f_calc_zc_main,
                      l_y_star = l_y_star_b_vec,
                      l_y_star_hat = l_y_star_hat_b_vec,
                      l_se = l_se_pred_star_b_vec,
                      MoreArgs = list(n_cat = ncol(m_true_pi_vec)),
                      SIMPLIFY = FALSE)
    
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
      
      if (length(valid_zc) < 2) { 
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


  # Bayesian Hierarchical MCMC Model (Beta Prior on Rho) ---
  if (l_methods$bayesian_mcmc_beta) {
    
    l_pred_int_bayesian_beta <- mapply(
      FUN = f_calc_pi_bayesian_mcmc,
      df_hist = l_x_k_vec,
      m = as.list(df_sim_settings$m),
      alpha = as.list(df_sim_settings$alpha),
      MoreArgs = list(stan_model = mod_beta),
      SIMPLIFY = FALSE
    )
    
    # Extract the individual lists of data frames
    l_pred_int_bayesian_beta_mean <- lapply(l_pred_int_bayesian_beta, `[[`, "mean_centered")
    l_pred_int_bayesian_beta_marg <- lapply(l_pred_int_bayesian_beta, `[[`, "marginal")
    l_pred_int_bayesian_beta_scs  <- lapply(l_pred_int_bayesian_beta, `[[`, "scs_rank")
    
    # --- 1. Bayesian MCMC beta (Mean Centered) ---
    l_res_bayesian_beta_mean <- lapply(seq_along(l_pred_int_bayesian_beta_mean), function(i) {
      if (is.null(l_pred_int_bayesian_beta_mean[[i]])) return(NA)
      df_counts <- l_y_vec[[i]]
      df_intervals <- l_pred_int_bayesian_beta_mean[[i]]
      category_cols <- setdiff(names(df_counts), "level")
      counts <- as.numeric(df_counts[1, category_cols])
      if (length(counts) != nrow(df_intervals)) stop(paste("Mismatch at index", i))
      all(counts >= df_intervals$lower & counts <= df_intervals$upper)
    })
    l_cwi_bayesian_beta_mean <- mapply(FUN = f_compare_with_intervals, data = l_y_vec, intervals = l_pred_int_bayesian_beta_mean, SIMPLIFY = FALSE)
    df_res_bayesian_beta_mean <- as.data.frame(do.call(rbind, l_res_bayesian_beta_mean))
    colnames(df_res_bayesian_beta_mean) <- "coverage_bayesian_mcmc_beta_mean"
    df_cwi_bayesian_beta_mean <- as.data.frame(do.call(rbind, l_cwi_bayesian_beta_mean))
    colnames(df_cwi_bayesian_beta_mean) <- paste0("bayesian_mcmc_beta_mean_", colnames(df_cwi_bayesian_beta_mean))
    df_width_stats_bayesian_beta_mean <- do.call(rbind, lapply(l_pred_int_bayesian_beta_mean, f_calc_interval_width_stats))
    colnames(df_width_stats_bayesian_beta_mean) <- paste0("bayesian_mcmc_beta_mean_", colnames(df_width_stats_bayesian_beta_mean))
    df_sim_res <- cbind(df_sim_res, df_res_bayesian_beta_mean, df_cwi_bayesian_beta_mean, df_width_stats_bayesian_beta_mean)
    
    # --- 2. Bayesian MCMC beta (Marginal) ---
    l_res_bayesian_beta_marg <- lapply(seq_along(l_pred_int_bayesian_beta_marg), function(i) {
      if (is.null(l_pred_int_bayesian_beta_marg[[i]])) return(NA)
      df_counts <- l_y_vec[[i]]
      df_intervals <- l_pred_int_bayesian_beta_marg[[i]]
      category_cols <- setdiff(names(df_counts), "level")
      counts <- as.numeric(df_counts[1, category_cols])
      if (length(counts) != nrow(df_intervals)) stop(paste("Mismatch at index", i))
      all(counts >= df_intervals$lower & counts <= df_intervals$upper)
    })
    l_cwi_bayesian_beta_marg <- mapply(FUN = f_compare_with_intervals, data = l_y_vec, intervals = l_pred_int_bayesian_beta_marg, SIMPLIFY = FALSE)
    df_res_bayesian_beta_marg <- as.data.frame(do.call(rbind, l_res_bayesian_beta_marg))
    colnames(df_res_bayesian_beta_marg) <- "coverage_bayesian_mcmc_beta_marg"
    df_cwi_bayesian_beta_marg <- as.data.frame(do.call(rbind, l_cwi_bayesian_beta_marg))
    colnames(df_cwi_bayesian_beta_marg) <- paste0("bayesian_mcmc_beta_marg_", colnames(df_cwi_bayesian_beta_marg))
    df_width_stats_bayesian_beta_marg <- do.call(rbind, lapply(l_pred_int_bayesian_beta_marg, f_calc_interval_width_stats))
    colnames(df_width_stats_bayesian_beta_marg) <- paste0("bayesian_mcmc_beta_marg_", colnames(df_width_stats_bayesian_beta_marg))
    df_sim_res <- cbind(df_sim_res, df_res_bayesian_beta_marg, df_cwi_bayesian_beta_marg, df_width_stats_bayesian_beta_marg)
    
    # --- 3. Bayesian MCMC beta (SCS Rank) ---
    l_res_bayesian_beta_scs <- lapply(seq_along(l_pred_int_bayesian_beta_scs), function(i) {
      if (is.null(l_pred_int_bayesian_beta_scs[[i]])) return(NA)
      df_counts <- l_y_vec[[i]]
      df_intervals <- l_pred_int_bayesian_beta_scs[[i]]
      category_cols <- setdiff(names(df_counts), "level")
      counts <- as.numeric(df_counts[1, category_cols])
      if (length(counts) != nrow(df_intervals)) stop(paste("Mismatch at index", i))
      all(counts >= df_intervals$lower & counts <= df_intervals$upper)
    })
    l_cwi_bayesian_beta_scs <- mapply(FUN = f_compare_with_intervals, data = l_y_vec, intervals = l_pred_int_bayesian_beta_scs, SIMPLIFY = FALSE)
    df_res_bayesian_beta_scs <- as.data.frame(do.call(rbind, l_res_bayesian_beta_scs))
    colnames(df_res_bayesian_beta_scs) <- "coverage_bayesian_mcmc_beta_scs"
    df_cwi_bayesian_beta_scs <- as.data.frame(do.call(rbind, l_cwi_bayesian_beta_scs))
    colnames(df_cwi_bayesian_beta_scs) <- paste0("bayesian_mcmc_beta_scs_", colnames(df_cwi_bayesian_beta_scs))
    df_width_stats_bayesian_beta_scs <- do.call(rbind, lapply(l_pred_int_bayesian_beta_scs, f_calc_interval_width_stats))
    colnames(df_width_stats_bayesian_beta_scs) <- paste0("bayesian_mcmc_beta_scs_", colnames(df_width_stats_bayesian_beta_scs))
    df_sim_res <- cbind(df_sim_res, df_res_bayesian_beta_scs, df_cwi_bayesian_beta_scs, df_width_stats_bayesian_beta_scs)
  }

  return(df_sim_res)
  
}

## Global Simulation Function ---------------------------------------------

f_global_sim <- function(m_true_pi_vec, df_sim_settings, l_methods) {
  
  # Determine the directory based on the number of columns in m_true_pi_vec
  ncols <- ncol(m_true_pi_vec)
  dir_path <- switch(as.character(ncols),
                     "3" = "./Results/intermediate_results/threecats/mult_sim_",
                     "5" = "./Results/intermediate_results/fivecats/mult_sim_",
                     "10" = "./Results/intermediate_results/tencats/mult_sim_",
                     stop("Invalid number of columns in m_true_pi_vec. Must be 3, 5, or 10."))
  
  l_prop <- replicate(df_sim_settings$nsim[1], m_true_pi_vec, simplify = FALSE)
  
  print(m_true_pi_vec)
  
  # Check if a parallel plan is configured
  if (is(future::plan(), "multisession") || is(future::plan(), "multicore") || is(future::plan(), "cluster")) {
    
    df_multsim_res <- NULL
    failures <- 0
    
    while (is.null(df_multsim_res)) {
      
      tryCatch({

        if (failures > 0) {
          message(paste("Retrying parallel computation... (Previous attempts failed", failures, "time(s))"))
        } else {
          message("Starting parallel computation with future_lapply...")
        }
        
        result_list <- future_lapply(l_prop, 
                                     f_predint_mult_sim, 
                                     df_sim_settings = df_sim_settings, 
                                     l_methods = l_methods, 
                                     mod_gamma = mod_gamma,
                                     mod_cauchy = mod_cauchy, 
                                     mod_beta = mod_beta,
                                     future.seed = TRUE)
        
        df_multsim_res <- do.call(rbind, result_list)
        
        message("Parallel computation completed successfully.")
        
      }, error = function(e) {
        
        failures <<- failures + 1
        
        if (failures > 1) {
          plan(sequential)
          Sys.sleep(5)
          plan(multisession)
        }
        
        message(paste("An error occurred. This was failure number", failures, "."))
        message("Error message: ", e$message) 
        message("Retrying in 10 seconds...")
        Sys.sleep(10) # Pause for 10 seconds before the next attempt
      })
    }
    
  } else {
    message("Starting sequential computation with lapply")
    df_multsim_res <-
      do.call(rbind,
              lapply(l_prop, 
                     f_predint_mult_sim, 
                     df_sim_settings, 
                     l_methods, 
                     mod_gamma,
                     mod_cauchy,
                     mod_beta
                    ))
  }
  
  print("dfdone")
  print(Sys.time())
  
  probsname <- paste(m_true_pi_vec, collapse = ".")
  df_multsim_res <- cbind(probs = probsname, df_multsim_res)
  
  file_path <- paste0(dir_path, probsname, ".rds")
  
  print(file_path)
  
  if (file.exists(file_path)) {
    boolFalse <- FALSE
    while (!boolFalse) {
      try({
        df_complete_simres <- readRDS(file_path)
        boolFalse <- TRUE
      }, silent = FALSE)
    }
    
    saveRDS(
      bind_rows(df_complete_simres, df_multsim_res),
      file_path,
      compress = FALSE
    )
  } else {
    saveRDS(df_multsim_res,
            file_path,
            compress = FALSE)
  }
}

