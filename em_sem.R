#library(VIM)
#data(sleep, package = "VIM")
#data(tao, package = "VIM")
#data(chorizonDL, package = "VIM")
#data(testdata, package = "VIM")

library(mixtools)

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
  
  # Return data and params
  return(params = list(mu = mu, cov = cov))
}

objective <- function(data, indicator, params) {
  # Calculate 2D log likelihood where x and y is observed
  const <- -log(2*pi)
  det <- -(1/2)*log(det(params$cov))
  maha <- -(1/2)*mahalanobis(data[!indicator,], params$mu, params$cov)
  log1 <- sum(const + det + maha)

  # Calculate 1D log likelihood where just x is observed
  const <- -(1/2)*log(2*pi)
  det <- -(1/2)*log(params$cov[1,1])
  maha <- -(1/2)*(data[indicator,][,1]-params$mu[1])^2/params$cov[1,1]
  log2 <- sum(const + det + maha)

  return(log1+log2)
}

em_algorithm <- function(data, tolerance) {
  
  params <- initialize(data)
  indicator <- is.na(y)
  vecLoglike <- NULL
  vecParams <- list(params)
  
  repeat {
    # EM step
    params <- emstep(data, indicator, params)
    vecParams <- append(vecParams, list(params))
    
    # Calculate objective (likelihood)
    vecLoglike <- append(vecLoglike, objective(data, indicator, params))
    
    # Check if convergence criterion is hit
    last <- length(vecLoglike)
    if(last >= 2 && abs(vecLoglike[last] - vecLoglike[last-1]) < tolerance) {
      break
    }
  }
  
  return(list(loglike = vecLoglike, params = vecParams))
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
result <- em_algorithm(data, tolerance = 10^-8)

# Plot convergence of log likelihood
plot(result$loglike, type = "l", xlab = "# of EM steps", ylab = "Log-Likelihood")
grid()

# Plot behavior of mu2
mu2 <- sapply(result$params, function(param){param$mu[2]})
plot(mu2, type = "l", xlab = "# of EM steps", ylab = "mu_2")
grid()

# contour plot
par(mfrow=c(1,4))
for(i in c(1,2,3,length(result$params))) {
  plot(data, pch = 20, xlim = c(-10,45), ylim = c(15,85),
       xlab = "x", ylab = "y", main = paste("EM iteration", i))
  grid()
  for(j in 1:9) ellipse(mu=result$params[[i]]$mu, sigma=result$params[[i]]$cov,
                        alpha = (j*0.1), npoints = 250, col=rgb(0,0,0,alpha=0.4))
}
par(mfrow=c(1,1))

#####################
### SEM Algorithm ###
#####################

calculateRatio <- function(data, params, i, j, tolerance) {
  indicator <- is.na(data[,2])
  vecR <- NULL
  
  # Get last parameter element (MLE)
  mle <- unlist(tail(params, 1)[[1]])
  
  t <- 1
  repeat {
    # Get current state of theta from em computation
    theta_t <- unlist(params[[t]])
    
    # Define theta_i as the final mle and replace i-th value with current state
    theta_t_i <- mle
    theta_t_i[[i]] <- theta_t[[i]]
    
    # Do EM step using theta_i, first convet parameters in list and then do step
    names(theta_t_i) <- NULL
    paramList <- list(mu = theta_t_i[1:2], cov = matrix(theta_t_i[3:6], nrow = 2, ncol = 2))
    theta_t_1_i <- unlist(emstep(data, indicator, paramList))
    
    # Calculate ratio
    #print(paste(theta_t_1_i[[j]],mle[[j]],theta_t[[i]],mle[[i]]))
    vecR <- append(vecR, (theta_t_1_i[[j]] - mle[[j]])/(theta_t[[i]] - mle[[i]]))
    
    # Check if convergence criterion is hit
    last <- length(vecR)
    if(last > 2 && abs(vecR[last] - vecR[last-1]) < tolerance) {
      break
    }
    t <- t+1
    if(t > length(result$params)) {
      break
    }
  }
  
  # Just return last rate
  return(vecR[length(vecR)])
}

calculateDM <- function(data, params, estimates, tolerance) {
  # Number of parameters to calculate variance for
  d <- length(estimates)
  
  # Define empty DM matrix
  DM <- matrix(nrow = d, ncol = d)
  
  # Calculate any r_ij and store in DM
  for(i in 1:d) {
    for(j in 1:d) {
      print(paste(i,j))
      DM[i,j] <- calculateRatio(data, params, estimates[i], estimates[j], tolerance)
    }
  }
  
  # Return whole DM
  return(DM)
}

sem_algorithm <- function(data, params, estimates, tolerance) {
  # Get DM matrix
  DM <- calculateDM(data, params, estimates, tolerance)
  
  # Calculate I_oc
  # ...
}

# 1: mu1, 2: mu2, 3: s_xx, 4/5: s_xy, 6: s_yy
DM <- calculateDM(data, params_trans, c(2,6), 10^-4)
#r_mu2 <- sem_algorithm(data, result$params, c(2,4,6), 10^-4)















