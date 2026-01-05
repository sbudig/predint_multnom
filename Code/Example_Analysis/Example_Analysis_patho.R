# Setup -------------------------------------------------------------------

library(tidyverse)
library(VGAM)
library(MCMCpack)
library(mvtnorm)
library(cmdstanr)
library(forcats)
library(xtable)
library(here)

# Functions ---------------------------------------------------------------

## Estimate Dispersion Parameter ------------------------------------------
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

## Function to calculate probabilities for multinomial model -------------
# calculates the probabilities for each category c
f_calc_pi <- function(mult_model) {
  # Extract the intercept coefficients
  coefs <- unname(coef(mult_model))
  
  # Compute the exponential of the intercepts
  exp_coefs <- exp(coefs)
  
  # Calculate probabilities using the softmax function
  probs <- exp_coefs / (1 + sum(exp_coefs))
  
  # Add the probability for the first category (which is 1 minus the sum of other probabilities)
  prob_first <- 1 - sum(probs)
  
  # Combine the probabilities for all categories
  return(c(prob_first, probs))
}

## Calculate the expected counts with pi and cluster size ----------------
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

## Fit a multinomial model function ----------------------------------------
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


## Internal function to calculate standard error of the historical  --------

f_calc_se_pred_star_internal <- function(phi_hat, pi_hat_vec, m, N_hist){
  
  # Variance of fut. random variable
  var_y_star <- m * phi_hat * pi_hat_vec * (1-pi_hat_vec)
  
  # Variance of fut. expectation
  var_y_star_hat <- (m^2 * phi_hat * pi_hat_vec * (1-pi_hat_vec))/N_hist
  
  # Prediction SE
  se_pred <- sqrt(var_y_star + var_y_star_hat)
  
  
  return(se_pred)
}

## Calculate the prediciton intervals -------------------------------------

f_calc_prediction_interval <- function(y_hat, q, se_pred, alternative)
{
  
  if (is.vector(q)){
    
    if (alternative == "both") {
      if (length(q) == 1) {
        lower <- pmax(0, y_hat - q * se_pred)
        upper <- y_hat + q * se_pred
        
        out <- data.frame(lower, upper)
        return(out)
      }
      
      if (length(q) == 2) {
        lower <- pmax(0, y_hat - q[1] * se_pred)
        upper <- y_hat + q[2] * se_pred
        
        out <- data.frame(lower, upper)
        return(out)
      }
      
    }
    
    if (alternative == "lower") {
      lower <- pmax(0, y_hat - q * se_pred)
      
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
      
      lower <- pmax(0, y_hat[i] - q[i,1] * se_pred[i])
      upper <- y_hat[i] + q[i,2] * se_pred[i]
      
      out <- rbind(out, data.frame(lower, upper))
      
    }
    return(out)
    
  } else {
    return(NULL)
  }
  
}

## Calculate adjusted q ---------------------------------------------------
# adjust q depending on number of n (Bonferroni adjustment)
f_adjust_q <- function(n, alpha) {
  qnorm(1 - (alpha/n)/ 2)
}

## Function to assess the coverage probability ----------------------------
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


## Bisection function -----------------------------------------------------
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
    
    # If another iteration is required,
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    if(sign(runval)==1){
      quant_min <- c}
    
    else if(sign(runval)==-1){
      quant_max <- c}
    
    
  }
  # If the max number of iterations is reached and no root has been found,
  # return message and end function.
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
  
  # Return NA if inputs are null
  if (is.null(phi_hat) || is.null(pi_hat)) {
    return(NA)
  }
  
  C <- length(pi_hat)
  
  # If there's only one category, it's not a multinomial problem.
  # Return standard normal quantile.
  if (C <= 1) {
    return(qnorm(1 - alpha / 2))
  }
  
  # Use a tryCatch block in case the correlation matrix is not positive definite
  q_mvn <- tryCatch({
    
    # Calculate the prediction covariance matrix
    Sigma_pi_hat <- diag(pi_hat) - pi_hat %*% t(pi_hat)
    Sigma_pred <- phi_hat * m * (1 + m / N_hist) * Sigma_pi_hat
    
    # Convert to a correlation matrix
    R_pred <- cov2cor(Sigma_pred)
    
    # Get the simultaneous critical value
    qmvnorm(1 - alpha, corr = R_pred, tail = "both.tails")$quantile
    
  }, error = function(e) {
    # If qmvnorm fails, return NA
    return(NA)
  })
  
  return(q_mvn)
}

f_calc_zc_help <- function(y_star, y_star_hat, se){
  zc <- (y_star-y_star_hat)/se
  return(zc)
}

f_calc_zc_main <- function(l_y_star, l_y_star_hat, l_se, n_cat) {
  
  # If any list is empty, return NULL
  if (length(l_y_star) == 0 || length(l_y_star_hat) == 0 || length(l_se) == 0) {
    return(NULL)
  }
  
  # Find indices of valid bootstrap replicates. A replicate is valid only if
  # its y_star, y_star_hat, and se all have the correct length (n_cat).
  valid_indices <- which(
    sapply(l_y_star, function(x) !is.null(x) && length(x) == n_cat) &
      sapply(l_y_star_hat, function(x) !is.null(x) && length(x) == n_cat) &
      sapply(l_se, function(x) !is.null(x) && length(x) == n_cat)
  )
  
  # If no valid replicates are found after filtering, return NULL.
  if (length(valid_indices) == 0) {
    return(NULL)
  }
  
  # Subset the lists to include only data from the valid replicates.
  l_y_star_valid <- l_y_star[valid_indices]
  l_y_star_hat_valid <- l_y_star_hat[valid_indices]
  l_se_valid <- l_se[valid_indices]
  
  # Now, perform the calculation only on the consistently-sized vectors.
  zc_list <- mapply(
    FUN = f_calc_zc_help,
    y_star = l_y_star_valid,
    y_star_hat = l_y_star_hat_valid,
    se = l_se_valid,
    SIMPLIFY = FALSE
  )
  
  # This rbind will now work because all elements in zc_list have the same length.
  do.call(rbind, zc_list)
  
}

#  Helper function to construct a simultaneous PI from posterior p --------
f_calc_pi_bayesian_standardized <- function(post_pred_samples, alpha) {
  
  # 1. Calculate moments of the posterior predictive distribution
  y_bar_bayes <- colMeans(post_pred_samples)
  y_sd_bayes  <- apply(post_pred_samples, 2, sd)
  
  # Handle zero variance if a category never occurs (prevent div by zero)
  y_sd_bayes[y_sd_bayes == 0] <- 1e-6 
  
  # 2. Calculate STANDARDIZED residuals (Pivotal Quantity) for every sample
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
    category_cols <- 1:(ncol(df_hist))
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
}

#  Helper function to construct a simultaneous PI from posterior p --------
# Helper function to construct a simultaneous PI from posterior predictive samples
f_construct_simultaneous_pi <- function(post_pred_samples, alpha) {
  
  # Step 1: Calculate the center of the predictive distribution
  y_hat_bayes <- colMeans(post_pred_samples)
  
  # Step 2: For each sample, calculate its maximum absolute deviation from the center
  # This is the "extremeness" score for each simulated vector
  abs_deviations <- abs(sweep(post_pred_samples, 2, y_hat_bayes, FUN = "-"))
  max_abs_devs <- apply(abs_deviations, 1, max)
  
  # Step 3: Find the (1-alpha) quantile of these scores
  # This is the critical value that defines the boundary of our simultaneous region
  q_crit <- quantile(max_abs_devs, probs = 1 - alpha, na.rm = TRUE)
  
  # Step 4: Construct the simultaneous prediction interval
  lower <- y_hat_bayes - q_crit
  upper <- y_hat_bayes + q_crit
  
  return(data.frame(lower, upper))
}

# Data --------------------------------------------------------------------

#set.seed(1577)
set.seed(206       )


df_pat <- data.frame(Minimal = c(15,18,6,4,7),
                     Slight = c(22,20,19,18,25),
                     Moderate = c(8,5,17,22,9),
                     Severe = c(1,0,4,2,0),
                     Massive = c(0,0,0,0,1)
                     )

fit_pat <- vglm(
  formula = cbind(Slight, Moderate, Severe, Massive, Minimal) ~ 1,
  family = multinomial,
  model = TRUE,
  data = df_pat
)

v_pi_hat <- f_calc_pi(fit_pat)
phi_hat <- f_estimate_phi(model = fit_pat, dispersion = "afroz")



xy_vec <- f_rdirmultinom(
  phi = phi_hat,
  k = 11,
  n = 46,
  v_pi = v_pi_hat
)

df_pat <- as.data.frame(xy_vec) %>%
  rename(
    Minimal = V1,
    Slight = V2,
    Moderate = V3,
    Severe = V4,
    Massive = V5
  ) %>%
  mutate(Trial = c( rep("Historical", times = 10), "Current"))
# 
# write.csv(df_pat, ".\\Code\\Example_Analysis\\df_pat.csv", row.names = FALSE)

# df_pat <- read.csv(".\\Code\\Example_Analysis\\df_pat.csv")

df_pat <- df_pat %>% rowwise() %>% mutate(Total = sum(c_across(Minimal:Massive))) %>% ungroup() %>% mutate(Study = c(1,seq(1:(nrow(df_pat)-1))))

df_pat_long <- df_pat %>% 
  pivot_longer(
    cols = c(Minimal, Slight, Moderate, Severe, Massive), 
    names_to = "Categories"
  ) %>%
  mutate(Proportion = value/Total)

# Define the order of the categories
category_order <- c("Minimal", "Slight", "Moderate", "Severe", "Massive")

# Apply the order to the 'Categories' column
df_pat_long <- df_pat_long %>%
  mutate(Categories = factor(Categories, levels = category_order),
         Trial = factor(Trial, levels = c("Historical", "Current")))

# Plots -------------------------------------------------------------------

plot_pat_trial <- df_pat_long %>%
  mutate(Trial = fct_relevel(Trial, "Historical", "Current")) %>%
  ggplot(
    aes(x = Categories)
  ) +
  geom_jitter(
    aes(y = value, color = factor(Study)), 
    height = 0, 
    width = 0.05,
    size = 4
  ) +
  facet_grid(~Trial) +
  theme_bw(base_size = 18) +
  labs(color = "Study", y = "Count", x = "Category")

plot_pat_trial

# ggsave("C:\\Users\\Budig\\Google_Drive\\Uni\\Phd\\03_Predint_multinomial\\Code\\Code_and_Data\\Figures\\plot_pat_trial.png",
#        plot = plot_pat_trial,
#        width = 12, height = 8,
#        dpi = 900)


# ggplot(df_pat_long, aes(x = Categories, y = value)) +
#   geom_col(fill = "steelblue") +
#   facet_wrap(~Study) + # Facets will also be in the correct order
#   labs(
#     title = "Patient Data by Study and Category",
#     x = "Injury Category",
#     y = "Count"
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Prediction Intervals ----------------------------------------------------

dat_pat_hist <- filter(df_pat, Trial == "Historical")
dat_pat_curr <- filter(df_pat, Trial == "Current")

fit_pat <- vglm(
  formula = cbind(Slight, Moderate, Severe, Massive, Minimal) ~ 1,
  family = multinomial,
  model = TRUE,
  data = dat_pat_hist
)

current_disp <- "afroz"

pi_hat_vec <- f_calc_pi(fit_pat)
phi_hat <- f_estimate_phi(model = fit_pat, dispersion = current_disp)

# the total number of observations across all historical studies is calculated.
# $N_{hist}=\sum_{k=1}^{K}n_k$
N_hist <- sum(rowSums(dat_pat_hist[c(1,2,3,4)]))

# number of historical studies
k_hist <- nrow(dat_pat_hist)

# historical cluster size (takes the maximum of all cluster sizes)
n_hist <- max(dat_pat_hist$Total)

m <- dat_pat_curr$Total

# Expected count for each category
# $\hat{\bm{y}}=(\bm{y}_1,...,\bm{y}_C$
y_hat_vec <- f_calc_y_hat(pi_hat_vec, m)

# Variance of future observation
var_y_vec <- f_calc_var_y(phi_hat, m, pi_hat_vec)

# Variance of the estimated future observation
var_y_hat_vec <- f_calc_var_y_hat(phi_hat, m, pi_hat_vec, N_hist)

# Prediction standard error for each category for each setting
se_pred_vec <- f_calc_se_pred(var_y_vec, var_y_hat_vec)

# Number of bootstrap iterations
B <- 1000

f_parboot <- function(pi_vec, phi, k, n, B, m, N_hist){
  
  current_phi <- phi
  current_k <- k
  current_n <- n
  current_pi <- pi_vec
  current_m <- m
  current_N_hist <- N_hist
  
  # Pre-allocate lists to store results from each bootstrap iteration
  pi_star_hat_b_list <- vector("list", B)
  phi_star_hat_b_list <- vector("list", B)
  y_star_hat_b_list <- vector("list", B)
  se_pred_star_b_list <- vector("list", B)
  y_star_b_list <- vector("list", B) # For bisection methods
  
  for (b in 1:B) {
    
    # Initialize a variable to hold the result for this iteration
    x_star_b <- NULL
    count_error_mbt <- 0
    # This loop will continue until x_star_b is successfully assigned a value
    while (is.null(x_star_b)) {
      
      tryCatch({
        # 1. Generate ONE historical bootstrap sample
        x_star_b <- as.data.frame(f_rdirmultinom(
          phi = current_phi,
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
      
      # 4. Generate ONE future bootstrap sample 
      y_star_b <- NULL
      count_error_mbty <- 0
      
      # This loop will continue until y_star_b is successfully assigned a value
      while (is.null(y_star_b)) {
        tryCatch({
          y_star_b <- as.data.frame(f_rdirmultinom(
            phi = current_phi,
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

bootstrap_results <- f_parboot(pi_hat_vec, phi_hat, k_hist, n_hist, B, m, N_hist)

# Unpack the results into the list structures your existing code expects
pi_star_hat_b_vec <- bootstrap_results[["pi_star_hat_b_vec"]]
phi_star_hat_b <- bootstrap_results[["phi_star_hat_b"]]
y_star_hat_b_vec <- bootstrap_results[["y_star_hat_b_vec"]]
se_pred_star_b_vec <-  bootstrap_results[["se_pred_star_b_vec"]]
y_star_b_vec <- bootstrap_results[["y_star_b_vec"]]

# set alpha level
alpha <- 0.05

# set alternative
alternative <- "both"

## Pointwise Normal Approximation -----------------------------------------
# use Pointwise Normal Approximation to compute prediction intervals

# calculate vector of quantiles depending on respective alphas
v_q_norm <- qnorm((1 - alpha / 2))

# calculate pointwise prediction intervals
pred_int_pointwise <- f_calc_prediction_interval(
  y_hat = y_hat_vec,
  q = v_q_norm,
  se_pred = se_pred_vec,
  alternative = alternative
)

pred_int_pointwise


# Multivariate Normal Approximation ---------------------------------------

q_mvn <- f_calc_q_mvn(
  phi_hat = phi_hat,
  pi_hat = pi_hat_vec,
  m = m,
  N_hist = N_hist,
  alpha = alpha)

# calculate prediction intervals with MVN-adjusted q
pred_int_mvn <- f_calc_prediction_interval(
  y_hat = y_hat_vec,
  q = q_mvn,
  se_pred = se_pred_vec,
  alternative = alternative
)

pred_int_mvn

## Bonferroni Adjustment --------------------------------------------------
# use Bonferroni Adjustment to compute prediction intervals

# list with elements of number of categories
n_cat <- length(y_hat_vec)

q_norm <- f_adjust_q(n_cat, alpha)

# calculate pointwise prediction intervals
pred_int_bonf <- f_calc_prediction_interval(
  y_hat = y_hat_vec,
  q = q_norm,
  se_pred = se_pred_vec,
  alternative = alternative
)

pred_int_bonf


## Symmetrical Bisection --------------------------------------------------
# use Symmetrical Bisection to compute prediction intervals

# calibrate q with bisection
q_calib_sym <- f_bisection(
  y_star_hat = y_star_hat_b_vec,
  se_pred = se_pred_star_b_vec,
  y_star = y_star_b_vec,
  alternative = alternative,
  quant_min = 0.01,
  quant_max = 10,
  n_bisec = 100,
  tol = 1e-3,
  alpha = alpha,
)

# calculate pointwise prediction intervals
pred_int_bisection_sym <- f_calc_prediction_interval(
  y_hat = y_hat_vec,
  q = q_calib_sym,
  se_pred = se_pred_vec,
  alternative = alternative
)

pred_int_bisection_sym


# Asymmetrical Bisection --------------------------------------------------

# calibrate q with bisection for lower border
alpha_calib_asym_lower <- f_bisection(
  y_star_hat = y_star_hat_b_vec,
  se_pred = se_pred_star_b_vec,
  y_star = y_star_b_vec,
  quant_min = 0.01,
  quant_max = 10,
  n_bisec = 100,
  tol = 1e-3,
  alpha = alpha / 2,
  alternative = "lower"
)

# calibrate q with bisection for upper boarder
alpha_calib_asym_upper <- f_bisection(
  y_star_hat = y_star_hat_b_vec,
  se_pred = se_pred_star_b_vec,
  y_star = y_star_b_vec,
  quant_min = 0.01,
  quant_max = 10,
  n_bisec = 100,
  tol = 1e-3,
  alpha = alpha / 2,
  alternative = "upper"
)

q_calib_asym <- c(alpha_calib_asym_lower, alpha_calib_asym_upper)

# calculate prediction intervals with calibrated q
pred_int_bisection_asym <- f_calc_prediction_interval(
  y_hat = y_hat_vec,
  q = q_calib_asym,
  se_pred = se_pred_vec,
  alternative = alternative
)

pred_int_bisection_asym


# Bonferroni Bisection ----------------------------------------------------

# list with elements of number of categories

alpha_calib_bonf_bisec <- f_single_bisec(
  n_cat = length(y_hat_vec),
  y_star_hat = y_star_hat_b_vec,
  se_pred = se_pred_star_b_vec,
  y_star = y_star_b_vec,
  quant_min = 0.01,
  quant_max = 10,
  n_bisec = 100,
  tol = 1e-3,
  alpha = alpha,
  alternative = "both"
)

# calculate prediction intervals with bonferroni adjusted q
pred_int_bonf_bisec <- f_calc_prediction_interval(
  y_hat = y_hat_vec,
  q = alpha_calib_bonf_bisec,
  se_pred = se_pred_vec,
  alternative = alternative
)

pred_int_bonf_bisec

# SCS Rank ----------------------------------------------------------------

df_zc <- f_calc_zc_main(
  l_y_star = y_star_b_vec,
  l_y_star_hat = y_star_hat_b_vec,
  l_se = se_pred_star_b_vec,
  n_cat = length(y_hat_vec)
)

scs_zc <- abs(SCSrank(as.matrix(df_zc))$conf.int)

# calculate prediction intervals with rectangular simultaneous confidence set
pred_int_scsrank <- f_calc_prediction_interval(
  y_hat = y_hat_vec,
  q = scs_zc,
  se_pred = se_pred_vec,
  alternative = alternative
)

pred_int_scsrank

# Studentized Bootstrap ---------------------------------------------------

df_zc <- f_calc_zc_main(
  l_y_star = y_star_b_vec,
  l_y_star_hat = y_star_hat_b_vec,
  l_se = se_pred_star_b_vec,
  n_cat = length(y_hat_vec)
)

# For each bootstrap replication (row), find the maximum absolute standardized residual.
max_abs_residuals <- apply(abs(df_zc), 1, max)

# The critical value 'q' is the (1-alpha) quantile of this distribution of maximums.
q_val <- quantile(max_abs_residuals, probs = 1 - alpha, na.rm = TRUE)


# Calculate prediction intervals with the new single critical value q
pred_int_stud <- f_calc_prediction_interval(
  y_hat = y_hat_vec,
  q = q_val,
  se_pred = se_pred_vec,
  alternative = alternative
)

pred_int_stud

# Bayesian Hierarchical Model ---------------------------------------------

## Beta Prior -------------------------------------------------------------
mod_beta <- cmdstan_model(here("Code", "dirichlet_multinomial_beta_rho.stan"), compile = TRUE)

pred_int_bayesian_beta <- f_calc_pi_bayesian_mcmc(
  df_hist = dat_pat_hist[, 1:5],
  m = m,
  alpha = alpha,
  stan_model = mod_beta)

pred_int_bayesian_beta_mean <- pred_int_bayesian_beta$mean_centered
pred_int_bayesian_beta_marg <- pred_int_bayesian_beta$marginal
pred_int_bayesian_beta_scs <- pred_int_bayesian_beta$scs_rank

## Half Cauchy Prior ------------------------------------------------------
mod_cauchy <- cmdstan_model(here("Code", "dirichlet_multinomial_cauchy.stan"), 
                            compile = TRUE)

pred_int_bayesian_cauchy <- f_calc_pi_bayesian_mcmc(
  df_hist = dat_pat_hist[, 1:5],
  m = m,
  alpha = alpha,
  stan_model = mod_cauchy)

pred_int_bayesian_cauchy_mean <- pred_int_bayesian_cauchy$mean_centered
pred_int_bayesian_cauchy_marg <- pred_int_bayesian_cauchy$marginal
pred_int_bayesian_cauchy_scs <- pred_int_bayesian_cauchy$scs_rank

# Result Table ------------------------------------------------------------

# 3. DEFINE CONSTANTS & INPUTS
# -------------------------------------------------------------------
# Define the order for categories and trials
category_order <- c("Minimal", "Slight", "Moderate", "Severe", "Massive")
trial_order    <- c("Historical", "Current")

# Observed/Predicted data for plotting
y_c_vec <- unlist(as.vector(dat_pat_curr[1:5]))

# Methods to include in the main plot
chosen_methods <- c(
  "Sym. Calibration", "Asym. Calibration", "Marg. Calibration",
  "Rank-Based SCS", "Max. Abs. Res.", "B: Cauchy (SCS)"
)

# chosen_methods <- c(
#   "Pointwise", "Bonferroni", "Sym. Calibration", "Asym. Calibration", "Marg. Calibration", 
#   "Rank-Based SCS", "Max. Abs. Res.", "B: Cauchy (SCS)"
# )

# Plotting aesthetics
dodge_width <- 0.6


# 4. MASTER DATA PREP: Consolidate all intervals
# -------------------------------------------------------------------
# Create a single named list of all interval data frames
all_intervals <- list(
  `Pointwise` = pred_int_pointwise,
  `MVN` = pred_int_mvn,
  `Bonferroni` = pred_int_bonf,
  `Sym. Calibration` = pred_int_bisection_sym,
  `Asym. Calibration` = pred_int_bisection_asym,
  `Marg. Calibration` = pred_int_bonf_bisec,
  `Rank-Based SCS` = pred_int_scsrank,
  `Max. Abs. Res.` = pred_int_stud,
  `B: Beta (mean)` = pred_int_bayesian_beta_mean,
  `B: Beta (marg)` = pred_int_bayesian_beta_marg,
  `B: Beta (SCS)` = pred_int_bayesian_beta_scs,
  `B: Cauchy (mean)` = pred_int_bayesian_cauchy_mean,
  `B: Cauchy (marg)` = pred_int_bayesian_cauchy_marg,
  `B: Cauchy (SCS)` = pred_int_bayesian_cauchy_scs

)

# 5. OUTPUT 1: LaTeX Table (for Supplementary Material)
# -------------------------------------------------------------------
# This section creates the wide table with all 10 methods
formatted_list <- lapply(all_intervals, function(df) {
  df <- as.data.frame(df)
  lower_bound <- pmax(df$lower, 0)
  upper_bound <- df$upper
  sprintf("[%.2f, %.2f]", lower_bound, upper_bound)
})

comparison_table <- as.data.frame(do.call(cbind, formatted_list))

# Get row names, or create them
cat_names <- rownames(all_intervals[[1]])
if (is.null(cat_names)) {
  cat_names <- paste("Category", 1:nrow(comparison_table))
}

# Add category names and predicted means
comparison_table <- data.frame(
  Category = cat_names,
  y = y_c_vec, # Using 'y_hat_vec' from assumptions
  comparison_table
)

# Generate xtable object
xtable_obj <- xtable(comparison_table,
                     caption = "Comparison of Simultaneous Prediction Intervals by Method",
                     label = "tab:pi_comparison",
                     align = c("c", "l", "r", rep("c", 14)) # l(eft), r(ight), c(enter)
)

# Print LaTeX code to console
print(xtable_obj,
      include.rownames = FALSE,
      booktabs = TRUE,
      caption.placement = "top",
      comment = FALSE
)


# 6. OUTPUT 2: Main Manuscript Plot (plot_combined_clean)
# -------------------------------------------------------------------

# 6.1. Create 'plot_data': The master long-format data for all intervals
plot_data_list <- lapply(names(all_intervals), function(method_name) {
  df <- as.data.frame(all_intervals[[method_name]])
  data.frame(
    Method = method_name,
    Category = factor(category_order, levels = category_order),
    Lower = pmax(df$lower, 0), # Apply 0-bound
    Upper = df$upper,
    y_hat_c = y_hat_vec,
    y_c = y_c_vec
  )
})
plot_data <- do.call(rbind, plot_data_list)

# 6.2. Prepare 'df_hist': Data for the "Historical" facet
df_hist <- df_pat_long %>%
  filter(Trial == "Historical") %>%
  mutate(
    Categories = factor(Categories, levels = category_order),
    Trial = factor(Trial, levels = trial_order)
  )

# 6.3. Prepare 'df_intervals': Data for the "Current" facet (Intervals)
df_intervals <- plot_data %>%
  filter(Method %in% chosen_methods) %>%
  mutate(
    Categories = Category, # Match column name
    Trial = factor("Current", levels = trial_order),
    Method = factor(Method, levels = chosen_methods)
  )

# 6.4. Prepare 'df_current_points': Data for "Current" facet (Points)
df_current_points <- df_intervals %>%
  group_by(Categories, Trial) %>%
  summarize(
    y_c = first(y_c),
    y_hat_c = first(y_hat_c),
    .groups = 'drop' # Drop grouping
  )

# 6.5. Generate the combined plot
plot_combined_clean <- ggplot() +
  
  # Layer 1: Historical Jitter
  geom_jitter(
    data = df_hist,
    aes(x = Categories, y = value, shape = factor(Study)),
    height = 0,
    width = 0.05,
    size = 3,
    color = "gray40"
  ) +
  
  # Layer 2: Prediction Intervals (Dodged)
  geom_errorbar(
    data = df_intervals,
    aes(
      x = Categories,
      ymin = Lower,
      ymax = Upper,
      color = Method
    ),
    width = 0.4,
    linewidth = 1.1,
    position = position_dodge(width = dodge_width)
  ) +
  
  # Layer 3: Observed Current Count (Central)
  geom_point(
    data = df_current_points,
    aes(x = Categories, y = y_c),
    shape = 4, # 'X'
    color = "red",
    size = 4,
    stroke = 1.5
  ) +
  
  # Layer 4: Predicted Mean (Central)
  geom_point(
    data = df_current_points,
    aes(x = Categories, y = y_hat_c),
    shape = 21, # Circle
    fill = "white",
    color = "black",
    size = 3,
    stroke = 1
  ) +
  
  # Layer 5: Facets, Theme, and Labels
facet_grid(
    ~Trial,
    scales = "free_x",
    space = "free_x",
    labeller = as_labeller(c("Historical" = "Historical Control Data", 
                             "Current" = "Current Trial"))
  ) +
  theme_bw(base_size = 10) +
  labs(
    y = "Count",
    x = "Category",
    shape = "Historical Study",
    color = "Prediction Interval Method"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  scale_shape_manual(values = 1:10) +
  guides(shape = "none")


print(plot_combined_clean)

ggsave(here("Figures", "plot_pat_trial.png"),
       plot = plot_combined_clean,
       width = 7,       # Full text width
       height = 5,      # Max text height
       units = "in",
       dpi = 900)     

# Alternative 1 -----------------------------------------------------------

# Create the plot
pi_plot <- ggplot(plot_data, aes(y = Method)) +
  
  # Add the prediction interval (as a line range)
  geom_linerange(aes(xmin = Lower, xmax = Upper), 
                 linewidth = 1.2, color = "gray50") +
  
  # Add the predicted mean (y_hat) as a point
  geom_point(aes(x = y_hat_c), 
             shape = 21, size = 2.5, fill = "white", color = "black") +
  
  # Add the observed value (y_c) as a red cross
  geom_point(aes(x = y_c), 
             shape = 4, size = 3, color = "red", stroke = 1.5) +
  
  # --- This is the most important part ---
  #   Create a separate plot for each Category
  facet_wrap(~ Category, scales = "free_x") +
  
  # --- Theming and Labels ---
  labs(
    title = "Comparison of 95% Simultaneous Prediction Intervals by Method and Category",
    x = "Observed/Predicted Count",
    y = "Method",
    caption = "Intervals [Lower, Upper]; Circle = Predicted (y_hat); Red Cross = Observed (y_c)"
  ) +
  theme_bw() + # A clean theme
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 12), # Facet titles
    legend.position = "none"
  )

# Display the plot
print(pi_plot)


# Alternative 2 -----------------------------------------------------------


# The corrected plot code
# plot_pat_with_intervals <- df_pat_long %>%
#   mutate(Trial = fct_relevel(Trial, "Historical", "Current")) %>%
#   ggplot(
#     aes(x = Categories) # Only map 'x' globally, since both layers use it
#   ) +
#   
#   # Layer 1: Prediction intervals (uses 'interval_df')
#   # This layer does NOT need or use 'value'
#   geom_linerange(
#     data = interval_df,
#     aes(
#       ymin = lower,
#       ymax = upper,
#       x = Categories
#     ),
#     alpha = 0.5,
#     color = "red",
#     linewidth = 1
#   ) +
#   
#   # Layer 2: Data points (uses the main 'df_pat_long')
#   # This layer gets its own 'y' and 'color' aesthetics
#   geom_jitter(
#     aes(y = value, color = factor(Study)), 
#     height = 0, 
#     width = 0.05,
#     size = 4
#   ) +
#   
#   # Apply faceting and theme as before
#   facet_grid(~Trial) +
#   theme_bw(base_size = 18) +
#   labs(color = "Study", y = "Count", x = "Category")
# 
# # Now this will print without error
# print(plot_pat_with_intervals)
