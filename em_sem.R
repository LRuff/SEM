#library(VIM)
#data(sleep, package = "VIM")
#data(tao, package = "VIM")
#data(chorizonDL, package = "VIM")
#data(testdata, package = "VIM")

library(MASS)
library(mixtools)
library(plyr)

####################
### EM Algorithm ###
####################

initialize <- function(data, type = "canonical") {
  x <- data[,1]
  y <- na.omit(data[,2])
  n <- length(x)
  
  # Initalize
  mu <- c(mean(x), mean(y))
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
  
  cov[1,2] <- cov[2,1] <-  (1/n)*(sum_xy - (sum_y*sum(x))/n)
  cov[2,2] <- (1/n)*(sum_yy - sum_y^2/n)
  
  # Return data and params
  return(params = list(mu = mu, cov = cov))
}

objective <- function(data, indicator, params) {
  # Calculate 2D log likelihood where x and y is observed
  if(sum(!indicator) > 0) {
    const <- -log(2*pi)
    det <- -(1/2)*log(det(params$cov))
    maha <- -(1/2)*mahalanobis(data[!indicator,], params$mu, params$cov)
    log1 <- sum(const + det + maha)
  } else {
    log1 <- 0
  }

  # Calculate 1D log likelihood where just x is observed
  if(sum(indicator) > 0) {
    const <- -(1/2)*log(2*pi)
    det <- -(1/2)*log(params$cov[1,1])
    maha <- -(1/2)*(data[indicator,][,1]-params$mu[1])^2/params$cov[1,1]
    log2 <- sum(const + det + maha)
  } else {
    log2 <- 0
  }

  return(log1+log2)
}

em_algorithm <- function(data, tolerance) {
  
  params <- initialize(data)
  indicator <- is.na(data[,2])
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
result <- em_algorithm(data, tolerance = 10^-12)

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
    
    # Do EM step using theta_i: convert parameters to list, then do step
    paramList <- list(mu = theta_t_i[1:2],
                      cov = matrix(theta_t_i[3:6], nrow = 2, ncol = 2))
    theta_t_1_i <- unlist(emstep(data, indicator, paramList))
    
    # Calculate ratio
    #print(paste(theta_t_1_i[[j]],mle[[j]],theta_t[[i]],mle[[i]]))
    vecR <- append(vecR, (theta_t_1_i[[j]] - mle[[j]])/(theta_t[[i]] - mle[[i]]))
    
    # Increase iteration
    t <- t+1
    
    # Check if convergence criterion is hit or we're running out of original estimations
    last <- length(vecR)
    if((last >= 2 && abs(vecR[last] - vecR[last-1]) < tolerance)) {
      break
    }
    if(t >= length(params)) {
      warning("SEM did not converge for one component.")
      break
    }
  }
  
  #print(vecR)
  
  # Just return last rate
  return(vecR[length(vecR)])
}

calculateDM <- function(data, params, tolerance) {
  # Parameters to estimate in DM*
  # 1: mu1, 2: mu2, 3: s_xx, 4/5: s_xy, 6: s_yy
  estimates <- c(2,4,6)
  
  # Number of parameters to calculate variance for
  d <- length(estimates)
  
  # Define empty DM matrix
  DM <- matrix(nrow = d, ncol = d)
  
  # Calculate any r_ij and store in DM
  for(i in 1:d) {
    for(j in 1:d) {
      DM[i,j] <- calculateRatio(data, params, estimates[i], estimates[j], tolerance)
      #print(paste(i,j,":",DM[i,j]))
    }
  }
  
  # Return whole DM
  return(DM)
}

sem_algorithm <- function(data, params, tolerance) {
  n <- length(data[,1])
  
  # Get DM* matrix
  DM <- calculateDM(data, params, tolerance)
  
  # Get covariance matrix of MLE estimate (last step from em algorithm)
  cov <- tail(params, 1)[[1]]$cov
  
  # Calculate G3 of I_oc
  G3_11 <- cov[2,2]
  G3_12 <- 0
  G3_13 <- 0
  G3_22 <- det(cov)^2*((cov[1,1]*cov[2,2])^2 - cov[1,2]^4)/
    (-cov[1,2]^6 + 3*cov[1,1]*cov[2,2]*cov[1,2]^4
     - 3*(cov[1,1]*cov[2,2]*cov[1,2])^2 + (cov[1,1]*cov[2,2])^3)
  G3_23 <- 2*cov[2,2]*cov[1,2]
  G3_33 <- 2*cov[2,2]^2
  
  G3 <- (1/n)*matrix(c(G3_11, G3_12, G3_13, G3_12, G3_22, G3_23,
                       G3_13, G3_23, G3_33), nrow = 3, ncol = 3)
  
  # Compute Delta V*
  DV <- G3%*%DM%*%solve(diag(3)-DM)
  
  return(G3 + DV)
}

#DM <- calculateDM(data, result$params, 10^-4)
V <- sem_algorithm(data, result$params, 10^-4)

##################
### Simulation ###
##################

simulation <- function(data, N, nSim, tolerance) {
  results <- vector("list", length = nSim)
  loglikes <- vector("list", length = nSim)
  mles <- matrix(3*nSim, nrow = nSim, ncol = 3)
  
  for(iSim in 1:nSim) {
    results[iSim] <- list(em_algorithm(data[((iSim-1)*N+1):(iSim*N),], tolerance^2))
    loglikes[iSim] <- list(results[iSim][[1]]$loglike)
    #print(str(results[iSim][[1]]$loglike))
    mles[iSim,] <- unlist(tail(results[[iSim]]$params, 1)[[1]])[c(2,4,6)]
  }
  
  # Calculate MC SD values
  #mcSd <- apply(mles, 2, function(col) {sd(col)})
  mcSd <- apply(mles, 2, function(col) {
    return(sqrt((1/nSim)*(sum(col^2) - sum(col)^2/nSim)))
  })
  names(mcSd) <- c("sd(mu2)","sd(cov12)","sd(cov22)")
  
  # Calculate SEM
  semVars <- vector("list", length = nSim)
  semVars <- matrix(3*nSim, nrow = nSim, ncol = 3)
  for(iSim in 1:nSim) {
    sem <- sem_algorithm(data[((iSim-1)*N+1):(iSim*N),], results[[iSim]]$params, tolerance)
    #print(sem)
    semVars[iSim,] <- suppressWarnings(sqrt(diag(sem)))
  }
  ## FIXME: Hack! diag(sem) should always be positive!!
  if(sum(is.nan(semVars)) > 0) {
    warning("Sem V matrix probably not positive semidefinite. Probably reached saddle point. Increase EM iterations.")
  }
  #semSd <- apply(semVars, 2, function(col) {mean(col[!is.nan(col)])})
  semSd <- apply(semVars, 2, function(col) {mean(col)})
  names(semSd) <- c("sd(mu2)","sd(cov12)","sd(cov22)")
  
  # Calculate average observed information values
  # ...
  
  return(list(loglikes = loglikes, mles = mles, sd = list(mc = mcSd, sem = semSd)))
}

nSim <- 1000 # 1000
N <- 250 # 250

data <- data_mcar <- data_mar <-
  mvrnorm(N*nSim, c(100,12), matrix(c(169,19.5,19.5,9), nrow = 2, ncol = 2))

# Generate MCAR data
na <- rbinom(N*nSim,1,0.5)
data_mcar[,2][as.logical(na)] <- NA

results_mcar <- simulation(data_mcar, N, nSim, 10^-4)

meanMu2 <- mean(results_mcar$mles[,1])
mcSdMu2 <- results_mcar$sd$mc[[1]]
plot(results_mcar$mles[,1], rep(1,nSim), pch = 16, col = rgb(0,0,0,0.1), cex = 0.8,
     ylim = c(-0.5,1.5), ylab = "", xlab = "", yaxt = 'n',
     main = "Standard deviation for mu_2")
axis(2, at = c(1,0), labels = c("MC", "SEM"))
arrows(meanMu2-mcSdMu2,1,meanMu2+mcSdMu2,1, col = "blue", code=3, length=0.1, angle = 90, lwd = 2)
points(meanMu2, 1, pch = 16, col = "red", cex = 1.2)

semSdMu2 <- results_mcar$sd$sem[[1]]
points(results_mcar$mles[,1], rep(0,nSim), pch = 16, col = rgb(0,0,0,0.1), cex = 0.8)
arrows(meanMu2-semSdMu2,0,meanMu2+semSdMu2,0, col = "blue", code=3, length=0.1, angle = 90, lwd = 2)
points(meanMu2, 0, pch = 16, col = "red", cex = 1.2)



# Plot MCAR result
logs <- ldply(results_mcar$loglikes, rbind)
plot(1, xlab = "Iteration", ylab = "Observed loglikelihood",
     xlim = c(1,dim(logs)[2]),
     ylim = c(min(logs, na.rm = TRUE), max(logs, na.rm = TRUE)))
for(i in 1:length(results_mcar$loglikes)) {
  lines(1:length(logs[i,]), logs[i,], col = rgb(0,0,1,0.1))
}
lines(apply(logs, 2, mean, na.rm = TRUE), col = "red")

# Generate MAR data
data_mar[,2][data_mar[,1] < median(data_mar[,1])] <- NA

results_mar <- simulation(data_mar, N, nSim, 10^-4)








