#library(VIM)
#data(sleep, package = "VIM")
#data(tao, package = "VIM")
#data(chorizonDL, package = "VIM")
#data(testdata, package = "VIM")

####################
### EM Algorithm ###
####################

initialize <- function(data) {
  x <- data[,1]
  y <- na.omit(data[,2])
  n <- length(x)
  
  # Initalize
  mu <- c(mean(x), mean(y))
  #cov <- matrix(c(var(x), 0, 0, var(y)), nrow = 2, ncol = 2) # Unbiased
  s_xx <- (1/n)*(sum(x^2) - sum(x)^2/n)
  s_yy <- (1/n)*(sum(y^2) - sum(y)^2/n)
  cov <- matrix(c(s_xx, 0, 0, s_yy), nrow = 2, ncol = 2)
  
  return(list(mu = mu, cov = cov))
}

emstep <- function(data, indicator, params) {
  mu <- params$mu
  cov <- params$cov
  x <- data[,1]
  y <- data[,2]
  n <- length(x)
  
  ## E step, update data
  
  # Regression estimation
  beta_1 <- cov[1,2]/cov[1,1]
  beta_0 <- mu[2] - beta_1*mu[1]
  sig_y_x <- cov[2,2] - beta_1^2*cov[1,1]
  
  # Imputation
  y[indicator] <- beta_0 + beta_1*x[indicator]
  
  # Estimate sufficient statistics
  sum_y <- sum(y)
  sum_xy <- sum(x*y)
  sum_yy <- sum(c(y[indicator]^2 + sig_y_x, y[!indicator]^2))
  
  ## M step, update parameters
  
  mu[2] <- mean(y)
  #cov[1,2] <- cov[2,1] <-  (1/(n-1))*(sum_xy - (sum_y*sum(x))/n) # Unbiased
  #cov[2,2] <- (1/(n-1))*(sum_yy - sum_y^2/n) # Unbiased
  
  cov[1,2] <- cov[2,1] <-  (1/n)*(sum_xy - (sum_y*sum(x))/n)
  cov[2,2] <- (1/n)*(sum_yy - sum_y^2/n)
  
  # Update data
  data[,2] <- y
  
  # Return data and params
  return(list(data = data, params = list(mu = mu, cov = cov)))
}

objective <- function(data, params) {
  # Calculate 2D log likelihood
  const <- -log(2*pi)
  det <- -(1/2)*log(det(params$cov))
  maha <- -(1/2)*mahalanobis(data, params$mu, params$cov)
  return(sum(const + det + maha))
}

em_algorithm <- function(data, tolerance = 0.00001) {
  
  params <- initialize(data)
  indicator <- is.na(y)
  loglike <- NULL
  vec_params <- list(params)
  
  repeat {
    # EM step
    result <- emstep(data, indicator, params)
    data <- result$data # Data for next step
    params <- result$params # Parameter for next step
    vec_params <- append(vec_params, list(params))
    
    # Calculate objective (likelihood)
    loglike <- append(loglike, objective(data, params))
    
    # Check if convergence criterion is hit
    last <- length(loglike)
    if(last > 2 && abs(loglike[last] - loglike[last-1]) < tolerance) {
      break;
    }
  }
  
  return(list(loglike = loglike, params = vec_params))
}

# Data

# Enders example
#x <- c(78,84,84,85,87,91,92,94,94,96,99,105,105,106,108,112,113,115,118,134)
#y <- c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,7,10,11,15,10,10,12,14,16,12)

# Paper example 3
x <- c(8,6,11,22,14,17,18,24,19,23,26,40,4,4,5,6,8,10)
y <- c(59,58,56,53,50,45,43,42,39,38,30,27,NA,NA,NA,NA,NA,NA)

data <- matrix(c(x,y), ncol = 2)

# Run algorithm
result <- em_algorithm(data, tolerance = 0.000001)

# Plot convergence of log likelihood
plot(result$loglike, type = "l", xlab = "# of EM steps", ylab = "Log-Likelihood")
grid()

# Plot behavior of mu2
mu2 <- sapply(result$params, function(param){param$mu[2]})
plot(mu2, type = "l", xlab = "# of EM steps", ylab = "mu_2")
grid()

#####################
### SEM Algorithm ###
#####################

sem_algorithm <- function(data, params) {
  # !! Do it just for mu2 !! TODO: Extend for general case
  indicator <- is.na(data[,2])
  r_mu <- NULL
  
  # MLE estimate of EM
  mle <- tail(params, 1)[[1]] # Get last parameter element
  
  for(i in 1:(length(params)-1)) { # TODO: replace by repeat and convergence criterion
    # Get current state of theta from em computation
    theta_t <- params[[i]]
    
    # Define theta_i as the final mle and replace i-th value with current state
    theta_t_i <- mle
    theta_t_i$mu[2] <- theta_t$mu[2]
    
    # Do EM step using theta_i
    res <- emstep(data, indicator, theta_t_i)
    theta_t_1_i <- res$params
    data <- res$data
    
    # Calculate ratio
    r_mu <- append(r_mu, (theta_t_1_i$mu[2] - mle$mu[2])/(params[[i]]$mu[2] - mle$mu[2]))
  }
  
  return(r_mu)
}

r_mu2 <- sem_algorithm(data, result$params)




















