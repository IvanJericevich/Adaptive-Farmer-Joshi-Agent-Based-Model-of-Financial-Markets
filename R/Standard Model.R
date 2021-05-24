# Title: Honours Project - Calibrating the Farmer-Joshi Agent-Based Model of Financial Markets
# Authors: Ivan Jericevich & Murray McKechnie
# Supervisor: A/Prof. Tim Gebbie

# 1. Load packages, data, and clean data
# 2. Get empirical moments and calculate weight matrix for the objective function
# 3. Load functions necessary for GA calibrations
# 4. Conduct GA calibrations, and calculate relevant confidence intervals
# 5. Load functions necessary for NMTA algorithm calibrations
# 6. Conduct NMTA algorithm calibrations, and calculate relevant confidence intervals
# 7. Plot results from simulations based on best parameters compared to empirical data
# 8. Other code - testing stability of objective function

### 1. Load packages, data, and clean data
####################################################################################################
rm(list = ls()) # Clear environment
if(!is.null(dev.list())) dev.off() # Clear plots
list.of.packages = c("openxlsx", "Rcpp", "timeSeries", "ggplot2", "timeDate", "pracma", "GA", "doParallel", "plotly", "dplyr", "forecast", "tseries", "fracdiff", "e1071", "extremefit") # This automatically installs and loads packages not already installed on the users computer
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])] # openxlsx = read.xlsx, timeSeries = timeSeries, timeDate = kurtosis, pracma = hurstexp, GA = ga, doParallel = implementation of ga in parallel, forecast = auto.arima
if(length(new.packages) > 0) {install.packages(new.packages)}
lapply(list.of.packages, require, character.only = T)

# Load and clean data
df = read.xlsx("BRICs-EOD-Compact-2005-2015.xlsx", sheet = 2, detectDates = TRUE)
df = df[, -c(length(colnames(df)), length(colnames(df)) - 1)] # Remove last two columns
colnames(df) = df[1, ]
df = df[-1, ]
names(df)[1] = "Date"
newColNames = strsplit(colnames(df), "J") # Clean column names
for (m in 3:length(newColNames)) {
  names(df)[m] = newColNames[[m]][1]
}
names(df)[2] = "JTOPI"
df = df[-c(nrow(df), nrow(df) - 1),] # Remove last two rows
df[, 1] = as.Date(as.numeric(df[, 1]), origin = "1899-12-30") # Convert excel numeric date format to dates
df[, 2:ncol(df)] = sapply(df[, 2:ncol(df)], as.numeric) # Convert character data to numeric
tsBRICS = timeSeries(df[, 2:ncol(df)], df[, 1]) # Convert data into time series object
setFinCenter(tsBRICS) = "Johannesburg"
price = na.omit(log(window(tsBRICS, start = "2005-01-01", end = "2016-04-29")$AGL))[-c(321)] # Log price of individual stock from 2013 to end
date = as.Date(time(tsBRICS)[which(!is.na(tsBRICS$AGL))])[-c(1:201, 321)]
returns = diff(price[-c(1:201)]) # Exclude data used to initialise model
####################################################################################################

### 2. Get empirical moments and calculate weight matrix for the objective function
####################################################################################################
# Calculate Empirical Moments
Moments = function(returns1, returns2) {
  mean_e = mean(returns1)
  stdev_e = sd(returns1)
  kurtosis_e = kurtosis(returns1)
  hurst_e = hurstexp(returns1, display = FALSE)$Hs
  ks_e = ks.test(x = returns1, y = returns2)$statistic # Kolmogorov-Smirnov statistic is zero since we are dealing with actual data. Represents difference between cummulative emperical distribution functions
  gph_e = fdGPH(abs(returns1))$d # GPH estimator representing long-range dependence
  adf_e = adf.test(returns1, k = 0)$statistic # ADF statistic representing random walk property of returns
  garch_e = sum(garch(returns1, trace = F)$coef[2:3]) # GARCH paramaters representing short-range dependence
  hill_e = mean(hill(returns1)$hill[(0.05 * length(returns1)):(0.1 * length(returns1))]) # Hill estimator
  return(c(mean_e, stdev_e, kurtosis_e, ks_e, hurst_e, gph_e, adf_e, garch_e, hill_e))
}
actual_moments <- Moments(returns, returns)

# Moving Block Bootsrap with overlapping blocks to calculate weight matrix (same as for adaptive model)

b = 100 # Size of window
k = 9 # Number of moments
bootstrap_moments = matrix(NA, nrow = 10000, ncol = k)

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

bootstrap_moments <- foreach(i = 1:10000, .combine = "rbind", 
                             .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit")) %dopar% {
  set.seed(i)
  index = sample(x = 1:(length(returns) - b + 1), 
                 size = ceiling(length(returns)/b), replace = TRUE) # Get block indices for bootstrap sample
  boot = c()
  for(j in 1:length(index)) { # Join the blocks together until we create a time series of at least length(p_returns)
    boot = c(boot, returns[(index[j]):(index[j]  + b - 1)])
  }
  boot = boot[1:length(returns)] # Truncate last block if necessary
  c(mean(boot), 
    sd(boot), 
    kurtosis(boot), 
    ks.test(boot, returns)$statistic, 
    hurstexp(boot, display = FALSE)$Hs,  
    fdGPH(abs(boot))$d, 
    adf.test(boot, k = 0)$statistic, 
    sum(garch(boot, trace = F)$coef[2:3]), 
    mean(hill(boot)$hill[(0.05 * length(boot)):(0.1 * length(boot))]))
}
stopCluster(cl)

W = solve(cov(bootstrap_moments)) # Inverse of covariance matrix from distribution of bootstrapped moments

####################################################################################################

### 3. Load functions necessary for GA calibrations
####################################################################################################
# Market-maker price setting function
MarketImpact = function(omega, lambda, p_t, zeta) {
  return(p_t + (1/lambda)*omega + zeta)
}

# Farmer-Joshi simulation function
Simulate = function(N, lambda, a, d_max, d_min, mu_eta, sigma_eta, sigma_zeta, 
                    T_max, tau_min, v_max, v_min, T_min = 0, tau_max = 0, I = 50) {
  # Initialization
  tau = matrix(NA, nrow = 2*N, ncol = I) # Each column refers to a different Monte-Carlo replication and each element is an agent
  p_sim = matrix(NA, nrow = length(price), ncol = I) # Matrix of I vectors of simulated prices
  d = sapply(X = rep(x = N, times = I), FUN = sample, x = d_min:d_max, replace = TRUE) # Matrix of I vectors of lags for chartists
  Thresh = sapply(X = rep(x = 2*N, times = I), FUN = runif, min = T_min, max = T_max) # Matrix of I vectors of entry thresholds (T) for agents
  for (i in 1:I) {
    tau_constraint = cbind(Thresh[, i], -tau_min) # Used to ensure -tau < Thresh
    tau[, i] = runif(2*N, -colMins(t(tau_constraint)), tau_max) # Exit thresholds (tau) for each agent for ith replication
  }
  c_value = a*(Thresh[1:N, ] - tau[1:N, ]) # Matrix of I vectors of c's for fundamentalists
  c_trend = a*(Thresh[(N + 1):(2*N), ] - tau[(N + 1):(2*N), ]) # Matrix of I vectors of c's for chartists
  p_sim[1:(200 + 2), ] = price[1:(200 + 2)] # The simulated and actual prices are the same for this interval so that the lagged price does not fall behind our window of data
  v_previous = sapply(X = rep(x = N, times = I), FUN = runif, min = v_min + p_sim[201], max = v_max + p_sim[201]) # The initial value for each value investor is exogenous
  
  # Inital mispricings
  m <- matrix(price[201], nrow = N, ncol = I, byrow = F) - v_previous # Mispricing that each fundamentlist believes there to be
  mc <-matrix(price[201 - d], nrow = N, ncol = I, byrow = F) - matrix(price[201], nrow = N, ncol = I, byrow = F)
  
  # Initial positions assuming an entry threshold is met
  x_previous_trend = c_trend*sign(matrix(price[201], nrow = N, ncol = I, byrow = TRUE) - matrix(price[201 - d], nrow = N, ncol = I, byrow = FALSE)) # The initial positions of chartists
  x_previous_value = c_value*sign(v_previous - matrix(price[201], nrow = N, ncol = I, byrow = TRUE)) # The initial positions of fundamentalists
  
  # Initial positions are only taken if outside of entry thresholds
  x_previous_value[which(m < Thresh[1:N, ] & m > -Thresh[1:N, ], arr.ind = TRUE)] = 0 # Initial position is 0
  x_previous_trend[which(mc < Thresh[1:N, ] & mc > -Thresh[1:N, ], arr.ind = TRUE)] = 0 # Initial position is 0
  
  for (t in (200 + 3):(length(price))) {
    # Chartists
    x_new_trend = c_trend*sign(matrix(p_sim[t - 1, ], nrow = N, ncol = I, byrow = TRUE) - 
                                 matrix(p_sim[as.vector(t - 1 - d + length(price) * matrix(0:(I - 1), nrow = N, ncol = I, byrow = T))], nrow = N, ncol = I, byrow = FALSE)) # The positions chartists want to now move to
    # Fundamentalists
    v_new = v_previous + sapply(X = rep(N, times = I), FUN = rnorm, mean = mu_eta, sd = sigma_eta) # The fundamentalists' new value judgements
    m = matrix(p_sim[t - 1, ], nrow = N, ncol = I, byrow = F) - v_new # Mispricing that each fundamentlist believes there to be
    mc <- matrix(p_sim[as.vector(t - 1 - d + length(price) * matrix(0:(I - 1), nrow = N, ncol = I, byrow = T))], 
                 nrow = N, ncol = I, byrow = F) - matrix(p_sim[t - 1, ], nrow = N, ncol = I, byrow = T)
    x_new_value = c_value*sign(v_new - matrix(p_sim[t - 1, ], nrow = N, ncol = I, byrow = TRUE)) # The position they wish to now take (assuming thresholds are met for fundamentalists)
    
    # State dependent value strategy
    # Fundamentalist was in a long position; will only trade if below tau
    x_new_value[which(x_previous_value > 0 & m < -tau[1:N, ], arr.ind = TRUE)] =
      x_previous_value[which(x_previous_value > 0 & m < -tau[1:N, ], arr.ind = TRUE)] # No trade takes place ()
    x_new_value[which(x_previous_value > 0 & m < Thresh[1:N, ] & m > -tau[1:N, ], arr.ind = TRUE)] = 0 # Exits long but does not enter short position as m doesn't meet entry threshold
    # Fundamentalist was in a short position; will only trade if above -tau
    x_new_value[which(x_previous_value < 0 & m > tau[1:N, ], arr.ind = TRUE)] = 
      x_previous_value[which(x_previous_value < 0 & m > tau[1:N, ], arr.ind = TRUE)] # No trade takes place
    x_new_value[which(x_previous_value < 0 & (m > -Thresh[1:N, ]) & (m < tau[1:N, ]), arr.ind = TRUE)] = 0 # Exits short but does not enter long position as m doesn't meet entry threshold
    # Fundamentalist was taking no position; will only trade if above Thresh or below -Thresh
    x_new_value[which(x_previous_value == 0 & m < Thresh[1:N, ] & m > -Thresh[1:N, ], arr.ind = TRUE)] = 0 # No trade takes place
    
    # Chartist was in a long position; will only trade if below tau
    x_new_trend[which(x_previous_trend > 0 & mc < -tau[(N+1):(2*N), ], arr.ind = TRUE)] =
      x_previous_trend[which(x_previous_trend > 0 & mc < -tau[(N+1):(2*N), ], arr.ind = TRUE)] # No trade takes place ()
    x_new_trend[which(x_previous_trend > 0 & mc < Thresh[(N+1):(2*N), ] & mc > -tau[(N+1):(2*N), ], arr.ind = TRUE)] = 0 # Exits long but does not enter short position as m doesn't meet entry threshold
    # Fundamentalist was in a short position; will only trade if above -tau
    x_new_trend[which(x_previous_trend < 0 & mc > tau[(N+1):(2*N), ], arr.ind = TRUE)] = 
      x_previous_trend[which(x_previous_trend < 0 & mc > tau[(N+1):(2*N), ], arr.ind = TRUE)] # No trade takes place
    x_new_trend[which(x_previous_trend < 0 & (mc > -Thresh[(N+1):(2*N), ]) & (mc < tau[(N+1):(2*N), ]), arr.ind = TRUE)] = 0 # Exits short but does not enter long position as m doesn't meet entry threshold
    # Fundamentalist was taking no position; will only trade if above Thresh or below -Thresh
    x_new_trend[which(x_previous_trend == 0 & mc < Thresh[(N+1):(2*N), ] & mc > -Thresh[(N+1):(2*N), ], arr.ind = TRUE)] = 0 # No trade takes place
    
    # Fundamentalists and chartists
    omega = colSums(x_new_trend - x_previous_trend) + colSums(x_new_value - x_previous_value) # The net order of all I sets of agents for the day
    v_previous <- v_new # Prepare for the next loop by updating v_previous
    x_previous_value <- x_new_value # Prepare for the next loop
    x_previous_trend <- x_new_trend # Prepare for the next loop
    p_sim[t, ] = MarketImpact(lambda = lambda, p_t = p_sim[t - 1, ], omega = omega, zeta = rnorm(I, 0, sigma_zeta)) # MarketImpact evaluated at a vector of parameters
  }
  return(p_sim) # I simulated price paths; First 202 days are real data
}

# Objective function
ObjectiveFunction = function(parameters) {
  I = 50 # Number of Monte Carlo replications
  set.seed(seed = seed + 100)
  p_sim = Simulate(N = ceiling(parameters[1]), # Ceiling and floor functions act as integer constraints
                   lambda = parameters[2], 
                   a = parameters[3], 
                   d_max = ceiling(parameters[5]) + floor(parameters[4]), # d_min + distance to d_max
                   d_min = ceiling(parameters[5]), 
                   mu_eta = parameters[6], 
                   sigma_eta = parameters[7], 
                   sigma_zeta = parameters[8], 
                   T_max = parameters[9] + parameters[13], 
                   tau_min = parameters[10], 
                   v_max = parameters[11] + parameters[12], # v_min + distance to v_max
                   v_min = parameters[12], 
                   T_min = parameters[13], 
                   tau_max = parameters[14] + parameters[10], 
                   I = 50)
  sim_returns <- diff(p_sim[-c(1:201), ])
  mean_sim = colMeans(sim_returns)
  stdev_sim = colSds(sim_returns)
  kurtosis_sim = apply(X = sim_returns, MARGIN = 2, FUN = kurtosis)
  hurst_sim = numeric(I)
  ks_sim = numeric(I)
  garch_sim = numeric(I)
  gph_sim = numeric(I)
  adf_sim = numeric(I)
  hill_sim <- numeric(I)
  for (i in 1:I)
  {
    if(max(abs(na.omit(sim_returns[, i]))) > 10 | sum(is.na(sim_returns[, i])) > 0) {
      return(NaN) # Prevents errors and quickly notices crazy behaviour
    } else {
      hurst_sim[i] = hurstexp(sim_returns[, i], display = FALSE)$Hs
      ks_sim[i] = ks.test(x = sim_returns[, i], y = returns)$statistic
      gph_sim[i] <- fdGPH(abs(sim_returns[, i]))$d
      adf_sim[i] <- adf.test(sim_returns[, i], k = 0)$statistic
      garch_sim[i] <- sum(garch(sim_returns[, i], trace = F)$coef[2:3])
      hill_sim[i] <- mean(hill(sim_returns[, i])$hill[(0.05 * length(sim_returns[, i])):(0.1 * length(sim_returns[, i]))])
    }
  }
  G = rbind(mean_sim, stdev_sim, kurtosis_sim, ks_sim, hurst_sim, 
            gph_sim, adf_sim, garch_sim, hill_sim) - 
    matrix(actual_moments, 
           ncol = I, nrow = k, byrow = FALSE)
  G = apply(G, 1, trimmean, percent = 70)
  return(-t(G)%*%W%*%G) # Returns a negative value as GA maximises
}

# Converts GA object to Simulate parameters
Decode1 = function(ga_object) { # Need to convert the optimised float parameters that are actually integers into the correct integers using ceiling and floor functions
  parameters = c(ceiling(summary(ga_object)$solution[1]), 
                 summary(ga_object)$solution[2], 
                 summary(ga_object)$solution[3], 
                 floor(summary(ga_object)$solution[4]) + ceiling(summary(ga_object)$solution[5]), 
                 ceiling(summary(ga_object)$solution[5]), 
                 summary(ga_object)$solution[6], 
                 summary(ga_object)$solution[7], 
                 summary(ga_object)$solution[8], 
                 summary(ga_object)$solution[9] + summary(ga_object)$solution[13], 
                 summary(ga_object)$solution[10], 
                 summary(ga_object)$solution[11] + summary(ga_object)$solution[12], 
                 summary(ga_object)$solution[12], 
                 summary(ga_object)$solution[13], 
                 summary(ga_object)$solution[14] + summary(ga_object)$solution[10])
  return(parameters)
}

# Converts ObjectiveFunction parameters to Simulate parameters
Decode2 = function(vectorparams) { # Need to convert the optimised float parameters that are actually integers into the correct integers using ceiling and floor functions
  parameters = c(ceiling(vectorparams[1]), 
                 vectorparams[2], 
                 vectorparams[3], 
                 floor(vectorparams[4]) + ceiling(vectorparams[5]), 
                 ceiling(vectorparams[5]), 
                 vectorparams[6], 
                 vectorparams[7], 
                 vectorparams[8], 
                 vectorparams[9] + vectorparams[13], 
                 vectorparams[10], 
                 vectorparams[11] + vectorparams[12], 
                 vectorparams[12], 
                 vectorparams[13], 
                 vectorparams[14] + vectorparams[10])
  return(parameters)
}

# Converts Simulate parameters to ObjectiveFunction parameters
Encode = function(vectorparams) { # Converts Simulate parameters to objective function parameters
  parameters = c(ceiling(vectorparams[1]), 
                 vectorparams[2], 
                 vectorparams[3], 
                 floor(vectorparams[4]) - ceiling(vectorparams[5]), 
                 ceiling(vectorparams[5]), 
                 vectorparams[6], 
                 vectorparams[7], 
                 vectorparams[8], 
                 vectorparams[9] - vectorparams[13], 
                 vectorparams[10], 
                 vectorparams[11] - vectorparams[12], 
                 vectorparams[12], 
                 vectorparams[13], 
                 vectorparams[14] - vectorparams[10])
  return(parameters)
}
####################################################################################################

### 4. Conduct GA calibrations, and calculate relevant confidence intervals
####################################################################################################

# GA Parameter CIs
# Run the GA many times for each calibration but 
# Only use a small number of iterations and set suggestedsol to the best solution from the previous one
parameters_CI = matrix(NA, nrow = 10, ncol = 14)
seed = 123
for(i in 1:10) {
  GA = ga(type = "real-valued", fitness = ObjectiveFunction, parallel = detectCores(), # parameters are: N, lambda, a, d_max-dmin, d_min, mu_eta, sigma_eta, sigma_zeta, T_max, tau_min, v_max-vmin, v_min
          lower = c(40, 0, 0, 0, 0.01, -0.01, 0, 0, 0, -1, 0, -0.5, 0, 0), seed = seed, popSize = 60, 
          upper = c(240, 15, 1.5, 100, 100, 0.01, 0.05, 0.05, 1, 0, 1, 0, 1, 1), maxiter = 25)
  
  seed = seed + 1
  
  GA = ga(type = "real-valued", fitness = ObjectiveFunction, parallel = detectCores(), # parameters are: N, lambda, a, d_max-dmin, d_min, mu_eta, sigma_eta, sigma_zeta, T_max, tau_min, v_max-vmin, v_min
          lower = c(40, 0, 0, 0, 0.01, -0.01, 0, 0, 0, -1, 0, -0.5, 0, 0), seed = seed,  popSize = 60, 
          upper = c(240, 15, 1.5, 100, 100, 0.01, 0.05, 0.05, 1, 0, 1, 0, 1, 1), maxiter = 15, 
          suggestions = summary(GA)$solution)
  
  seed = seed + 1
  
  GA = ga(type = "real-valued", fitness = ObjectiveFunction, parallel = detectCores(), # parameters are: N, lambda, a, d_max-dmin, d_min, mu_eta, sigma_eta, sigma_zeta, T_max, tau_min, v_max-vmin, v_min
          lower = c(40, 0, 0, 0, 0.01, -0.01, 0, 0, 0, -1, 0, -0.5, 0, 0), seed = seed,  popSize = 60, 
          upper = c(240, 15, 1.5, 100, 100, 0.01, 0.05, 0.05, 1, 0, 1, 0, 1, 1), maxiter = 15, 
          suggestions = summary(GA)$solution)
  
  seed = seed + 1
  
  parameters_GA[i, ] <- Decode1(GA)
  
  write.table(cbind(c("N", "lambda", "a", "d_max", "d_min", "mu_eta", "sigma_eta", "sigma_zeta", "T_max", "tau_min", "v_max", "v_min", "T_min", "tau_max"), parameters_CI), file = "Parameters_CI.txt", append = FALSE, quote = FALSE, sep = ",", dec = ".", row.names = FALSE)
}

seed = 153 # Keep updating seed

for(i in 1:10) { # Keep optimising each calibration run until fitness exceeds -40
  while(fitness[i] < -40) {
    GA = ga(type = "real-valued", fitness = ObjectiveFunction, parallel = detectCores(), # parameters are: N, lambda, a, d_max-dmin, d_min, mu_eta, sigma_eta, sigma_zeta, T_max, tau_min, v_max-vmin, v_min
            lower = c(40, 0, 0, 0, 0.01, -0.01, 0, 0, 0, -1, 0, -0.5, 0, 0), seed = seed,  popSize = 60, 
            upper = c(240, 15, 1.5, 100, 100, 0.01, 0.05, 0.05, 1, 0, 1, 0, 1, 1), maxiter = 15, 
            suggestions = Encode(parameters_GA[i, ]))
    
    parameters_GA[i, ] <- Decode1(GA)
    fitness_GA[i] <- summary(GA)$fitness
    seed = seed + 1
  }
  write.table(cbind(c("N", "lambda", "a", "d_max", "d_min", "mu_eta", "sigma_eta", "sigma_zeta", "T_max", "tau_min", "v_max", "v_min", "T_min", "tau_max"), parameters_CI), file = "Parameters_CI.txt", append = FALSE, quote = FALSE, sep = ",", dec = ".", row.names = FALSE)
}

par_GA <- parameters_GA[1:10, -1]
names(par_CI) <- c("N", "lambda", "a", "d_max", "d_min", "mu_eta", "sigma_eta", "sigma_zeta", "T_max", "tau_min", "v_max", "v_min", "T_min", "tau_max")

x_bar <- colMeans(par_GA)
samp_sd <- colSds(par_GA)

par_GA_lower_CI <- x_bar + qt(0.025, (dim(par_GA)[1] - 1)) * samp_sd/sqrt(dim(par_GA)[1])
par_GA_upper_CI <- x_bar - qt(0.025, (dim(par_GA)[1] - 1)) * samp_sd/sqrt(dim(par_GA)[1])

# GA Moment CIs + GA Price Path CIs
set.seed(1)
simulation_price = do.call(Simulate, as.list(c(par_CI[7, ], 1000)))[-c(1:201), ]
sim_returns = diff(simulation_price)

mean_sim = colMeans(sim_returns)
stdev_sim = colSds(sim_returns)
kurtosis_sim = apply(X = sim_returns, MARGIN = 2, FUN = kurtosis)
hurst_sim = numeric(1000)
ks_sim = numeric(1000)
garch_sim = numeric(1000)
gph_sim = numeric(1000)
adf_sim = numeric(1000)
hill_sim <- numeric(1000)
for (i in 1:1000) {
  hurst_sim[i] = hurstexp(sim_returns[, i], display = FALSE)$Hs
  ks_sim[i] = ks.test(x = sim_returns[, i], y = returns)$statistic
  gph_sim[i] <- fdGPH(abs(sim_returns[, i]))$d
  adf_sim[i] <- adf.test(sim_returns[, i], k = 0)$statistic
  garch_sim[i] <- sum(garch(sim_returns[, i], trace = F)$coef[2:3])
  hill_sim[i] <- mean(hill(sim_returns[, i])$hill[(0.05 * length(sim_returns[, i])):(0.1 * length(sim_returns[, i]))])
}

mean_sim = sort(mean_sim)
stdev_sim = sort(stdev_sim)
kurtosis_sim = sort(kurtosis_sim)
hurst_sim = sort(hurst_sim)
ks_sim = sort(ks_sim)
garch_sim = sort(garch_sim)
gph_sim = sort(gph_sim)
adf_sim = sort(adf_sim)
hill_sim <- sort(hill_sim)

lower_bounds_GA <- c(mean_sim[26], stdev_sim[26], kurtosis_sim[26], ks_sim[26], hurst_sim[26], 
                  garch_sim[26], gph_sim[26], adf_sim[26], hill_sim[26])

upper_bounds_GA <- c(mean_sim[975], stdev_sim[975], kurtosis_sim[975], ks_sim[975], hurst_sim[975], 
                  garch_sim[975], gph_sim[975], adf_sim[975], hill_sim[975])

sorted_prices_GA <- t(apply(simulation_price, 1, sort)) # Sort simulated prices for use in price path CIs in section 7 below
####################################################################################################

### 5. Load functions necessary for NMTA calibrations
####################################################################################################
# Threshold accepting portion of algorithm
threshold_accepting = function(vertices, fitnesses) {
  tau <- c(3, 1.5, 0) # Each round, we reduce how much worse the best fitness can be than the previous solution. These numbers will almost certainly need tweaking.
  n <- length(fitnesses) - 1 # Number of parameters
  
  # Now, we try shifting a bunch of parameters one at a time at random, keeping whatever improves the fitness (or at least makes it worse by less than tau)
  for(i in 1:3) { # Each round we reduce tau and increase the number of replications
    for(j in 1:7) { # Within each round we shift up to 5 parameters
      direction <- sample(size = 1, x = 1:n) # Determine which parameter to shift
      magnitude <- runif(n = 1, min = -0.5, max = 0.5) * mean(vertices[direction, ]) # Determine how much to shift it by
      # Get new candidate solutions
      temp_vertices <- vertices
      temp_vertices[direction, ] <- temp_vertices[direction, ] + magnitude
      
      # Restrict shifts to possible values only
      temp_vertices[1, ] <- pmax(1, temp_vertices[1, ])
      temp_vertices[2, ] <- pmax(0.001, temp_vertices[2, ])
      temp_vertices[3, ] <- pmax(0, temp_vertices[3, ])
      temp_vertices[4, ] <- pmax(0, temp_vertices[4, ])
      temp_vertices[4, ] <- pmin(100, temp_vertices[4, ])
      temp_vertices[5, ] <- pmax(1, temp_vertices[5, ])
      temp_vertices[5, ] <- pmin(90, temp_vertices[5, ])
      temp_vertices[7, ] <- pmax(0, temp_vertices[7, ])
      temp_vertices[8, ] <- pmax(0.00000001, temp_vertices[8, ])
      temp_vertices[9, ] <- pmax(0, temp_vertices[9, ])
      temp_vertices[10, ] <- pmin(0, temp_vertices[10, ])
      temp_vertices[11, ] <- pmax(0, temp_vertices[11, ])
      temp_vertices[14, ] <- pmax(0, temp_vertices[14, ])
      # Get their fitnesses
      temp_fitnesses <- foreach(k = 1:(n+1), .combine = "c", .export = ls(globalenv()), 
                                .packages = list.of.packages) %dorng% {
                                  ObjectiveFunction(c(temp_vertices[, k], 15 + 5 * i))
                                }
      # Replace the old solutions with the new ones if the best solution from the new set is good enough
      if(max(temp_fitnesses) > max(fitnesses) - tau[i]) {
        vertices <- temp_vertices
        fitnesses <- temp_fitnesses
      }
    }
  }
  return(rbind(vertices, fitnesses))
}

# Full NMTA algorithm - requires a maximum number of iterations and a set number of replications
Calibrate_NMTA = function(maxiter, I) { # Based off of Platt. 'I' is the number of replications
  
  # Get starting vertices and fitnesses
  initials <- foreach(k = 1:15, .combine = "cbind", .export = ls(globalenv()), 
                      .packages = list.of.packages) %dorng% {
                        vertex <- c(runif(1, 40, 240), runif(1, 0, 15), runif(1, 0, 1.5), 
                                    runif(1, 0, 100), runif(1, 0.01, 90), runif(1, -0.01, 0.01), 
                                    runif(1, 0, 0.05), runif(1, 0, 0.05), runif(1, 0, 1), 
                                    runif(1, -1, 0), runif(1, 0, 1), runif(1, -0.5, 0), 
                                    runif(1, 0, 1), runif(1, 0, 1))
                        fitness <- ObjectiveFunction(c(vertex, I))
                        while(fitness == -100000) { # Ensure that simulation does not result in crazy behaviour
                          vertex <- c(runif(1, 40, 240), runif(1, 0, 15), runif(1, 0, 1.5), 
                                      runif(1, 0, 100), runif(1, 0.01, 90), runif(1, -0.01, 0.01), 
                                      runif(1, 0, 0.05), runif(1, 0, 0.05), runif(1, 0, 1), 
                                      runif(1, -1, 0), runif(1, 0, 1), runif(1, -0.5, 0), 
                                      runif(1, 0, 1), runif(1, 0, 1))
                          fitness <- ObjectiveFunction(c(vertex, I))
                        }
                        c(fitness, vertex)
                      }
  
  vertices <- initials[-1, ]
  fitnesses <- initials[1, ]
  
  # Sort from best to worst
  vertices <- vertices[, order(fitnesses, decreasing = T)]
  fitnesses <- sort(fitnesses, decreasing = T)
  
  # Keep track of best ever result
  best_ever_iteration <- 0
  best_ever_fitness <- fitnesses[1]
  best_ever_vertex <- vertices[, 1]
  
  progress <- c(0, fitnesses[1], vertices[, 1])
  names(progress) <- c("Iteration", "fitness", "N", "lambda", "a", "damx - dmin", "dmin", "mu_eta", "sigma_eta", "sigma_zeta", "T_max - T_min", "tau_min", "v_max - vmin", "v_min", "Tmin", "tau_max - tau_min")
  print(progress)
  
  # Set the calibration parameters. These all might need tweaking.
  alpha <- 0.15 # Probability of choosing threshold accepting section
  sigmaNM <- 0.5 # Shrinkage parameter
  rhoNM <- 0.5 # Contraction parameter
  alphaNM <- 1 # Reflection parameter
  gammaNM <- 2 # Expansion parameter
  
  for (r in 1:maxiter) {
    if (runif(1, 0, 1) < alpha) { # Threshold accepting
      new_parameters <- threshold_accepting(vertices, fitnesses)
      vertices <- new_parameters[1:(14), ]
      fitnesses <- new_parameters[15, ]
    }
    else {
      # Reflect the worst point to the other side of the average of each of the other parameters
      centroid <- rowMeans(vertices[, -15])
      reflected_point <- centroid + alphaNM * (centroid - vertices[, 15])
      
      # Restrict shifts to possible values only
      reflected_point[1] <- max(1, reflected_point[1])
      reflected_point[2] <- max(0.001, reflected_point[2])
      reflected_point[3] <- max(0, reflected_point[3])
      reflected_point[4] <- max(0, reflected_point[4])
      reflected_point[4] <- min(100, reflected_point[4])
      reflected_point[5] <- max(1, reflected_point[5])
      reflected_point[5] <- min(90, reflected_point[5])
      reflected_point[7] <- max(0, reflected_point[7])
      reflected_point[8] <- max(0.00000001, reflected_point[8])
      reflected_point[9] <- max(0, reflected_point[9])
      reflected_point[10] <- min(0, reflected_point[10])
      reflected_point[11] <- max(0, reflected_point[11])
      reflected_point[14] <- max(0, reflected_point[14])
      
      ref_fitness <- ObjectiveFunction(c(reflected_point, I))
      # Check if new solution is in between best and worst fitnesses
      if(ref_fitness < fitnesses[1] & ref_fitness > fitnesses[14]) {
        vertices[, 15] <- reflected_point
        fitnesses[15] <- ref_fitness
      }
      # Check if new parameters are new best solution; if so try improve further by moving more in that direction
      else if(ref_fitness > fitnesses[1]) {
        # Expand
        expanded_point <- centroid + gammaNM * (reflected_point - centroid)
        
        # Restrict shifts to possible values only
        expanded_point[1] <- max(1, expanded_point[1])
        expanded_point[2] <- max(0.001, expanded_point[2])
        expanded_point[3] <- max(0, expanded_point[3])
        expanded_point[4] <- max(0, expanded_point[4])
        expanded_point[4] <- min(100, expanded_point[4])
        expanded_point[5] <- max(1, expanded_point[5])
        expanded_point[5] <- min(90, expanded_point[5])
        expanded_point[7] <- max(0, expanded_point[7])
        expanded_point[8] <- max(0.00000001, expanded_point[8])
        expanded_point[9] <- max(0, expanded_point[9])
        expanded_point[10] <- min(0, expanded_point[10])
        expanded_point[11] <- max(0, expanded_point[11])
        expanded_point[14] <- max(0, expanded_point[14])
        
        exp_fitness <- ObjectiveFunction(c(expanded_point, I))
        # Check if the solution is even better now
        if(exp_fitness > ref_fitness) {
          vertices[, 15] <- expanded_point
          fitnesses[15] <- exp_fitness
        } else { # If not then just use the original improvement
          vertices[, 15] <- reflected_point
          fitnesses[15] <- ref_fitness
        }
      }
      # Check if new parameters are still the worst. If so, try contracting
      else if(ref_fitness < fitnesses[14]) {
        contracted_point <- centroid + rhoNM * (vertices[, 15] - centroid)
        con_fitness <- ObjectiveFunction(c(contracted_point, I))
        # Check if contracting managed to actually get a better result than the original worst parameter
        if(con_fitness > fitnesses[15]) {
          vertices[, 15] <- contracted_point
          fitnesses[15] <- con_fitness
        } else {
          # If none of the above helped, shrink everything towards best solution
          vertices <- vertices[, 1] + sigmaNM * (vertices - vertices[, 1])
          fitnesses <- c(fitnesses[1], foreach(k = 2:15, .combine = "c", .export = ls(globalenv()), 
                                               .packages = list.of.packages) %dorng% {
                                                 ObjectiveFunction(c(vertices[, k], I))
                                               })
        }
      }
    }
    # Sort from best to worst and try again
    vertices <- vertices[, order(fitnesses, decreasing = T)]
    fitnesses <- sort(fitnesses, decreasing = T)
    
    if(fitnesses[1] > best_ever_fitness) { # Keep track of best ever result
      best_ever_iteration <- r
      best_ever_fitness <- fitnesses[1]
      best_ever_vertex <- vertices[, 1]
    }
    
    progress <- c(r, fitnesses[1], vertices[, 1])
    names(progress) <- c("Iteration", "fitness", "N", "lambda", "a", "damx - dmin", "dmin", "mu_eta", "sigma_eta", "sigma_zeta", "T_max - T_min", "tau_min", "v_max - vmin", "v_min", "T_min", "tau_max - tau_min")
    print(progress)
  }
  return(c(best_ever_iteration, best_ever_fitness, best_ever_vertex))
}

# Slightly altered ObjectiveFunction
ObjectiveFunction = function(parameters) {
  I = parameters[15] # Number of Monte Carlo replications
  p_sim = Simulate(N = ceiling(parameters[1]), # Ceiling and floor functions act as integer constraints
                   lambda = parameters[2], 
                   a = parameters[3], 
                   d_max = ceiling(parameters[5]) + floor(parameters[4]), # d_min + distance to d_max
                   d_min = ceiling(parameters[5]), 
                   mu_eta = parameters[6], 
                   sigma_eta = parameters[7], 
                   sigma_zeta = parameters[8], 
                   T_max = parameters[9] + parameters[13], 
                   tau_min = parameters[10], 
                   v_max = parameters[11] + parameters[12], # v_min + distance to v_max
                   v_min = parameters[12], 
                   T_min = parameters[13], 
                   tau_max = parameters[14] + parameters[10], 
                   I = parameters[15])
  sim_returns <- diff(p_sim[-c(1:201), ])
  mean_sim = colMeans(sim_returns)
  stdev_sim = colSds(sim_returns)
  kurtosis_sim = apply(X = sim_returns, MARGIN = 2, FUN = kurtosis)
  hurst_sim = numeric(I)
  ks_sim = numeric(I)
  garch_sim = numeric(I)
  gph_sim = numeric(I)
  adf_sim = numeric(I)
  hill_sim <- numeric(I)
  for (i in 1:I)
  {
    if(max(abs(na.omit(sim_returns[, i]))) > 10 | sum(is.na(sim_returns[, i])) > 0) {
      return(NaN) # Prevents errors and quickly notices crazy behaviour
    } else {
      hurst_sim[i] = hurstexp(sim_returns[, i], display = FALSE)$Hs
      ks_sim[i] = ks.test(x = sim_returns[, i], y = returns)$statistic
      gph_sim[i] <- fdGPH(abs(sim_returns[, i]))$d
      adf_sim[i] <- adf.test(sim_returns[, i], k = 0)$statistic
      garch_sim[i] <- sum(garch(sim_returns[, i], trace = F)$coef[2:3])
      hill_sim[i] <- mean(hill(sim_returns[, i])$hill[(0.05 * length(sim_returns[, i])):(0.1 * length(sim_returns[, i]))])
    }
  }
  G = rbind(mean_sim, stdev_sim, kurtosis_sim, ks_sim, hurst_sim, 
            gph_sim, adf_sim, garch_sim, hill_sim) - 
    matrix(actual_moments, 
           ncol = I, nrow = k, byrow = FALSE)
  G = apply(G, 1, trimmean, percent = 70)
  return(-t(G)%*%W%*%G)
}

# Simplified Hurst Exponent function - taken directly from hurstexp function
rssimple <- function(x, d = 50) {
  
  stopifnot(is.numeric(x), is.numeric(d))
  d <- max(2, floor(d[1]))
  N <- length(x)
  if (N%%2 != 0) {
    x <- c(x, (x[N - 1] + x[N])/2)
    N <- N + 1
  }
  
  divisors <- function(n, n0 = 2) {
    n0n <- n0:floor(n/2)
    dvs <- n0n[n%%n0n == 0]
    return(dvs)
  }
  
  N <- length(x)
  dmin <- d
  N0 <- min(floor(0.99 * N), N - 1)
  N1 <- N0
  dv <- divisors(N1, dmin)
  
  for (i in (N0 + 1):N) {
    dw <- divisors(i, dmin)
    if (length(dw) > length(dv)) {
      N1 <- i
      dv <- dw
    }
  }
  
  OptN <- N1
  d <- dv
  x <- x[1:OptN]
  n <- length(x)
  y <- x - mean(x)
  s <- cumsum(y)
  rs <- (max(s) - min(s))/sd(x)
  
  log(rs)/log(n)
}
####################################################################################################

### 6. Conduct NMTA calibrations, and calculate relevant confidence intervals
####################################################################################################
## NMTA Parameter CIs
# Conduct 15 calibrations
cl <- makeCluster(detectCores())
registerDoParallel(cl)
set.seed(1)
cal1full <- Calibrate_NMTA(65, 50)
cal1 <- Decode2(cal1full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
set.seed(2)
cal2full <- Calibrate_NMTA(65, 50)
cal2 <- Decode2(cal2full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores()) # replace
registerDoParallel(cl)
set.seed(3)
cal3full <- Calibrate_NMTA(65, 50)
cal3 <- Decode2(cal3full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
set.seed(4)
cal4full <- Calibrate_NMTA(65, 50)
cal4 <- Decode2(cal4full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
set.seed(5)
cal5full <- Calibrate_NMTA(65, 50)
cal5 <- Decode2(cal5full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores()) # replace
registerDoParallel(cl)
set.seed(6)
cal6full <- Calibrate_NMTA(65, 50)
cal6 <- Decode2(cal6full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
set.seed(7)
cal7full <- Calibrate_NMTA(65, 50)
cal7 <- Decode2(cal7full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores()) # replace
registerDoParallel(cl)
set.seed(8)
cal8full <- Calibrate_NMTA(65, 50)
cal8 <- Decode2(cal8full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores()) #replace
registerDoParallel(cl)
set.seed(9)
cal9full <- Calibrate_NMTA(65, 50)
cal9 <- Decode2(cal9full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
set.seed(10)
cal10full <- Calibrate_NMTA(65, 50)
cal10 <- Decode2(cal10full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
set.seed(11)
cal11full <- Calibrate_NMTA(65, 50)
cal11 <- Decode2(cal11full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
set.seed(12)
cal12full <- Calibrate_NMTA(65, 50)
cal12 <- Decode2(cal12full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
set.seed(13)
cal13full <- Calibrate_NMTA(65, 50)
cal13 <- Decode2(cal13full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
set.seed(14)
cal14full <- Calibrate_NMTA(65, 50)
cal14 <- Decode2(cal14full[3:16])
stopCluster(cl)

cl <- makeCluster(detectCores())
registerDoParallel(cl)
set.seed(15)
cal15full <- Calibrate_NMTA(65, 50)
cal15 <- Decode2(cal15full[3:16])
stopCluster(cl)

# Retain best 10 calibrations
parameters_NMTA <- rbind(cal1, cal2, cal4, cal5, cal7, cal10, cal11, cal12, cal13, cal15)[, 3:16]

# Calculate CIs
x_barNMTA <- colMeans(parameters_NMTA)
samp_sdNMTA <- colSds(parameters_NMTA)

lower_CI <- x_barNMTA + qt(0.025, (dim(parameters_NMTA)[1] - 1)) * samp_sdNMTA/sqrt(dim(parameters_NMTA)[1])
upper_CI <- x_barNMTA - qt(0.025, (dim(parameters_NMTA)[1] - 1)) * samp_sdNMTA/sqrt(dim(parameters_NMTA)[1])


## NMTA Moment CIs + NMTA Price Path CIs
set.seed(1)
simulation_priceNM = do.call(Simulate, as.list(c(cal5[3:16], 1000)))[-c(1:201), ]
sim_returnsNM = diff(simulation_priceNM)

mean_simNM = colMeans(sim_returnsNM)
stdev_simNM = colSds(sim_returnsNM)
kurtosis_simNM = apply(X = sim_returnsNM, MARGIN = 2, FUN = kurtosis)
hurst_simNM = numeric(1000)
ks_simNM = numeric(1000)
garch_simNM = numeric(1000)
gph_simNM = numeric(1000)
adf_simNM = numeric(1000)
hill_simNM = numeric(1000)
for (i in 1:1000) {
  hurst_simNM[i] = hurstexp(sim_returnsNM[, i], display = FALSE)$Hs
  ks_simNM[i] = ks.test(x = sim_returnsNM[, i], y = returns)$statistic
  gph_simNM[i] <- fdGPH(abs(sim_returnsNM[, i]))$d
  adf_simNM[i] <- adf.test(sim_returnsNM[, i], k = 0)$statistic
  garch_simNM[i] <- sum(garch(sim_returnsNM[, i], trace = F)$coef[2:3])
  hill_simNM[i] <- mean(hill(sim_returnsNM[, i])$hill[(0.05 * length(sim_returnsNM[, i])):(0.1 * length(sim_returnsNM[, i]))])
}

mean_simNM = sort(mean_simNM)
stdev_simNM = sort(stdev_simNM)
kurtosis_simNM = sort(kurtosis_simNM)
hurst_simNM = sort(hurst_simNM)
ks_simNM = sort(ks_simNM)
garch_simNM = sort(garch_simNM)
gph_simNM = sort(gph_simNM)
adf_simNM = sort(adf_simNM)
hill_simNM <- sort(hill_simNM)

lower_boundsNM <- c(mean_simNM[26], stdev_simNM[26], kurtosis_simNM[26], ks_simNM[26], hurst_simNM[26], 
                  garch_simNM[26], gph_simNM[26], adf_simNM[26], hill_simNM[26])

upper_boundsNM <- c(mean_simNM[975], stdev_simNM[975], kurtosis_simNM[975], ks_simNM[975], hurst_simNM[975], 
                  garch_simNM[975], gph_simNM[975], adf_simNM[975], hill_simNM[975])

sorted_prices_NMTA <- t(apply(simulation_priceNM, 1, sort)) # Sort simulated prices for use in price path CIs in section 7 below
####################################################################################################

### 7. Plot results from simulations based on best parameters compared to empirical data
####################################################################################################
## Get GA and NMTA simulations for graphs
set.seed(1) # Simulate using GA optimised parameters
simulation_price_GA = do.call(Simulate, as.list(c(par_CI[7, ], 1)))[-c(1:201)]
simulation_returns_GA = diff(simulation_price_GA)
Moments(simulation_returns, returns) # Point estimates of moments for GA calibration 

set.seed(1) # Simulate using NMTA optimised parameters
simulation_price_NMTA = do.call(Simulate, as.list(c(Decode2(cal5[-(1:2)]), 1)))[-c(1:201)]
simulation_returns_NMTA = diff(simulation_price_NMTA)
Moments(simulation_returns, returns) # Point estimates of moments for NMTA algorithm calibration

## Plot Graphs
# Log Price
pdf("Price paths.pdf", width = 10, height = 5)
data.frame(simulation_price = c(price[-c(1:201)], simulation_price_GA, simulation_price_NMTA), type = rep(c("Actual", "GA Simulation", "NMTA Simulation"), each = length(simulation_price_GA)), date = rep(date, times = 3)) %>%
  ggplot(data = ., mapping = aes(x = date, y = c(price[-c(1:201)], simulation_price_GA, simulation_price_NMTA), colour = rep(c("Actual", "GA Simulation", "NMTA Simulation"), each = length(simulation_price_GA)))) +
  geom_line() +
  scale_colour_manual(name = "Legend", values = c(Actual = "black", `GA Simulation` = "blue", `NMTA Simulation` = "red")) +
  geom_ribbon(data = data.frame(date = rep(date, times = 3), up = rep(sorted_prices_GA[,975], 3), low = rep(sorted_prices_GA[,26], 3)), 
              aes(ymin = low, ymax = up), fill="blue", alpha="0.1", colour = NA) + 
  geom_ribbon(data = data.frame(date = rep(date, times = 3), up = rep(sorted_prices_NMTA[,975], 3), low = rep(sorted_prices_NMTA[,26], 3)), 
              aes(ymin = low, ymax = up), fill="red", alpha="0.1", colour = NA) + 
  labs(x = "Date", y = "Log Price")
dev.off()

# Log Returns
pdf("Return paths.pdf", width = 10, height = 5)
data.frame(returns = c(returns, simulation_returns_GA, simulation_returns_NMTA), date = rep(date[-1], times = 3), type = rep(c("Actual", "GA Simulation", "NMTA Simulation"), each = length(simulation_returns_GA))) %>%
  ggplot(data = ., mapping = aes(x = date, y = returns, color = type)) +
  geom_line() +
  scale_colour_manual(name = "Legend", values = c(Actual = "black", `GA Simulation` = "blue", `NMTA Simulation` = "red")) +
  facet_wrap(~ type, nrow = 3) +
  labs(x = "Date", y = "Log Returns")
dev.off()

# Normal probability plots for log returns
pdf("Normal probability plots.pdf", width = 10, height = 5)
data.frame(returns = c(scale(returns), scale(simulation_returns_GA), scale(simulation_returns_NMTA)), type = rep(c("Actual", "GA Simulation", "NMTA Simulation"), each = length(simulation_returns_GA))) %>%
  ggplot(data = ., mapping = aes(sample = returns)) +
  stat_qq_point() +
  stat_qq_line(col = "black") +
  stat_qq_band(fill = rep(c("black", "blue", "red"), each = length(returns)), alpha = 0.4) + # alpha = transparency
  facet_wrap(~ type, ncol = 3) + # Returns should be log normal
  labs(x = "Theoretical", y = "Sample")
dev.off()

# Autocorrelation of log returns
pdf("Autocorrelation of log returns.pdf", width = 10, height = 5)
cbind(rbind(with(acf(scale(returns), type = "covariance", lag.max = 50, plot = FALSE), data.frame(lag, acf)), with(acf(scale(simulation_returns_GA), type = "covariance", lag.max = 50, plot = FALSE), data.frame(lag, acf)), with(acf(scale(simulation_returns_NMTA), type = "covariance", lag.max = 50, plot = FALSE), data.frame(lag, acf))), type = rep(c("Actual", "GA Simulation", "NMTA Simulation"), each = 51)) %>%
  ggplot(data = ., mapping = aes(x = lag, y = acf)) +
  geom_point(mapping = aes(colour = type)) +
  geom_segment(mapping = aes(x = lag, y = 0, xend = lag, yend = acf)) +
  geom_hline(yintercept = c(-1, 1)*qnorm((1 + 0.95)/2)/sqrt(length(returns)), lty = 2, col = "blue") +
  scale_colour_manual(name = "Legend", values = c(Actual = "black", `GA Simulation` = "blue", `NMTA Simulation` = "red")) +
  facet_wrap(~ type, ncol = 3) +
  labs(x = "Lag", y = "Autocorrelation")
dev.off()

# Autocorrelation of absolute log returns
pdf("Autocorrelation of absolute log returns.pdf", width = 10, height = 5)
cbind(rbind(with(acf(scale(abs(returns)), type = "covariance", lag.max = 50, plot = FALSE), data.frame(lag, acf)), with(acf(scale(abs(simulation_returns_GA)), type = "covariance", lag.max = 50, plot = FALSE), data.frame(lag, acf)), with(acf(scale(abs(simulation_returns_NMTA)), type = "covariance", lag.max = 50, plot = FALSE), data.frame(lag, acf))), type = rep(c("Actual", "GA Simulation", "NMTA Simulation"), each = 51)) %>%
  ggplot(data = ., mapping = aes(x = lag, y = acf)) +
  geom_point(mapping = aes(colour = type)) +
  geom_segment(mapping = aes(x = lag, y = 0, xend = lag, yend = acf)) +
  geom_hline(yintercept = c(-1, 1)*qnorm((1 + 0.95)/2)/sqrt(length(returns)), lty = 2, color = "blue") +
  scale_colour_manual(name = "Legend", values = c(Actual = "black", `GA Simulation` = "blue", `NMTA Simulation` = "red")) +
  facet_wrap(~ type, ncol = 3) +
  labs(x = "Lag", y = "Autocorrelation")
dev.off()
####################################################################################################

### 8. Other code - testing stability of the objective function
####################################################################################################
cl <- makeCluster(detectCores())
registerDoParallel(cl)

Objective_Variance <- foreach(i = 1:10, .combine = "c", 
                              .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
                                trim = 10 * (i - 1)
                                values = c()
                                for(j in 1:10) {
                                  seed = j
                                  values = c(values, ObjectiveFunction(summary(GA)$solution))
                                }
                                var(values)
                              }

stopCluster(cl)

### 9. Save workspace
save.image("Jericevich and McKechnie Honours Project - Standard Model.RData")
