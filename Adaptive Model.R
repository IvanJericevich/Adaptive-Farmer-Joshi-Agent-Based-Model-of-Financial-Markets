# Title: Calibrating an Adabtive Farmer-Joshi Agent-Based Model of Financial Markets
# Authors: Ivan Jericevich & Murray McKechnie
# Supervisor: A/Prof. Tim Gebbie
# 1. Set working directory and load packages
# 2. Supplementary functions
# 3. Load and prepare data sets
# 4. Get empirical moments
# 5. Apply moving block bootstrap with overlapping blocks to obtain distribution of actual moments
# 6. Simulate price path for given set of parameters
#     - Set up dataframe keeping track of all agents' strategies
#     - Obtain omega for chartists using state dependent trend strategy
#     - Obtain omega for fundamentalists using state dependent value strategy
#     - Update positions, profits and strategies
#     - Allow for agents to adapt their strategy according to the profits made from each strategy
# 7. Construct objective function using simulated price paths
#     - Calculate simulated moments using Monte Carlo method
#     - Return deviations from actual moments weighted by the uncertainty in each moment (W)
# 8. Implement genetic algorithm
#     - Apply GA to price paths
#     - Decode parameters given by GA so they can be entered into simulate function
#     - Get confidence intervals for price paths, parameters and moments
# 9. Implement Nelder-Mead with threshold accepting algorithm
#     - Apply NMTA to price paths
#     - Get confidence intervals for price paths, parameters and moments
# 10. Store, visualise results and test replication of stylised facts about financial markets
#     - Log price paths
#     - Normal probability plots
#     - Autocorrelation of log returns plots
#     - Autocorrelation of absolute log returns plots
#     - Objective surface plots
# 11. Save workspace
####################################################################################################

### 1. Preliminaries
cat("\014") # Clear console
rm(list = ls()) # Clear environment
if(!is.null(dev.list())) dev.off() # Clear plots
setwd("~/Honours Project")
load("Workspace1.RData") # Load data from previous workspace
list.of.packages = c("openxlsx", "timeSeries", "ggplot2", "timeDate", "pracma", "GA", "doParallel", "plotly", "dplyr", "forecast", "tseries", "fracdiff", "e1071", "extremefit", "latex2exp", "qqplotr", "doRNG", "plot3D", "pracma", "TTR") # This automatically installs and loads packages not already installed on the users computer
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])] # openxlsx = read.xlsx, timeSeries = timeSeries, timeDate = kurtosis, pracma = hurstexp, GA = ga, doParallel = foreach, forecast = auto.arima, tseries = garch, fracdiff = fdGPH, e1071, extremefit = hill, latex2exp = TeX, qqplotr = stat_qq_point, doRNG = seed in ga, plot3D = persp3D, pracma = trimmean, doSNOW = registerDoSNOW, TTR = SMA
if(length(new.packages) > 0) {install.packages(new.packages)}
lapply(list.of.packages, require, character.only = T)
####################################################################################################

### 2. Supplementary Functions
Decode1 = function(ga_object) { # Need to convert the optimised float parameters that are actually integers into the correct integers using ceiling and floor functions
  parameters = c(ceiling(summary(ga_object)$solution[1]), # N
                 summary(ga_object)$solution[2], # lambda
                 summary(ga_object)$solution[3], # a
                 floor(summary(ga_object)$solution[4]) + ceiling(summary(ga_object)$solution[5]), # d_max
                 ceiling(summary(ga_object)$solution[5]), # d_min
                 summary(ga_object)$solution[6], # mu_eta
                 summary(ga_object)$solution[7], # sigma_eta
                 summary(ga_object)$solution[8], # sigma_zeta
                 summary(ga_object)$solution[9] + summary(ga_object)$solution[15], # T_max
                 summary(ga_object)$solution[10], # tau_min
                 summary(ga_object)$solution[11] + summary(ga_object)$solution[12], # v_max
                 summary(ga_object)$solution[12], # v_min
                 summary(ga_object)$solution[13], # Gamma
                 ceiling(summary(ga_object)$solution[14]), # H
                 summary(ga_object)$solution[15], # T_min
                 summary(ga_object)$solution[16] + summary(ga_object)$solution[10]) # tau_max
  return(parameters)
}
Decode2 = function(vectorparams) { # Need to convert the optimised float parameters that are actually integers into the correct integers using ceiling and floor functions
  parameters = c(ceiling(vectorparams[1]), 
                 vectorparams[2], 
                 vectorparams[3], 
                 floor(vectorparams[4]) + ceiling(vectorparams[5]), 
                 ceiling(vectorparams[5]), 
                 vectorparams[6], 
                 vectorparams[7], 
                 vectorparams[8], 
                 vectorparams[9] + vectorparams[15], 
                 vectorparams[10], 
                 vectorparams[11] + vectorparams[12], 
                 vectorparams[12], 
                 vectorparams[13], 
                 ceiling(vectorparams[14]), 
                 vectorparams[15], 
                 vectorparams[16] + vectorparams[10])
  return(parameters)
}
####################################################################################################

### 3. Data Preparation
## BRICS data
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
tsBRICS = timeSeries(df[, 2:ncol(df)], df[, 1]) # Convert data into time series object 467
setFinCenter(tsBRICS) = "Johannesburg"
price = na.omit(log(window(tsBRICS, start = "2005-01-01", end = "2016-04-29")$AGL))[-321] # Log price of individual stock from 2013 to end
date = as.Date(time(tsBRICS)[which(!is.na(tsBRICS$AGL))])[-c(1:201, 321)]
returns = diff(price[-c(1:201)]) # Exclude data used to initialise model
####################################################################################################

### 4. Emperical moments of returns
Moments = function(returns1, returns2) {
  mean_e = mean(returns1)
  stdev_e = sd(returns1)
  kurtosis_e = kurtosis(returns1)
  ks_e = ks.test(x = returns1, y = returns2)$statistic # Kolmogorov-Smirnov statistic is zero since we are dealing with actual data. Represents difference between cummulative emperical distribution functions
  hurst_e = hurstexp(returns1, display = FALSE)$Hs
  gph_e = fdGPH(abs(returns1))$d # GPH estimator representing long-range dependence
  adf_e = adf.test(returns1, k = 0)$statistic # ADF statistic representing random walk property of returns
  garch_e = sum(garch(returns1, trace = F)$coef[2:3]) # GARCH paramaters representing short-range dependence
  hill_e = mean(hill(returns1)$hill[(0.05 * length(returns1)):(0.1 * length(returns1))]) # Hill estimator
  return(c(mean_e, stdev_e, kurtosis_e, ks_e, hurst_e, gph_e, adf_e, garch_e, hill_e))
}
####################################################################################################

### 5. Moving Block Bootstrap with overlapping blocks
b = 100 # Size of window
k = 9 # Number of moments
bootstrap_moments = matrix(NA, nrow = 10000, ncol = k)
cl = makeCluster(detectCores())
registerDoParallel(cl)
bootstrap_moments = foreach(i = 1:10000, .combine = "rbind", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit")) %dopar% {
  set.seed(i)
  index = sample(x = 1:(length(returns) - b + 1), size = ceiling(length(returns)/b), replace = TRUE) # Get block indices for bootstrap sample
  boot = c()
  for(j in 1:length(index)) { # Join the blocks together until we create a time series of at least length(returns)
    boot = c(boot, returns[(index[j]):(index[j]  + b - 1)])
  }
  boot = boot[1:length(returns)] # Truncate last block if necessary
  c(mean(boot), sd(boot), kurtosis(boot), ks.test(boot, returns)$statistic, hurstexp(boot, display = FALSE)$Hs, fdGPH(abs(boot))$d, adf.test(boot, k = 0)$statistic, sum(garch(boot, trace = F)$coef[2:3]), mean(hill(boot)$hill[(0.05 * length(boot)):(0.1 * length(boot))]))
}
stopCluster(cl)
W = solve(cov(bootstrap_moments)) # Inverse of covariance matrix from distribution of bootstrapped moments
write.table(cbind(c("mean", "stdev", "kurtosis", "ks", "hurst", "gph", "adf", "garch", "hill"), round(W, 3)), file = "W.txt", append = FALSE, quote = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Moment", "mean", "stdev", "kurtosis", "hurst", "ks", "gph", "adf", "garch", "hill"))
####################################################################################################

### 6. Simulation
MarketImpact = function(omega, lambda, p_t, zeta) {
  return(p_t + (1/lambda)*omega + zeta)
}
Profit = function(p_cur, p_prev, x_prev) {
  x_prev*(p_cur - p_prev) # Difference between trade value today and the value of the order placed yesterday which comes into effect today 
}
Simulate = function(price = price, N, lambda, a, d_max, d_min, mu_eta, sigma_eta, sigma_zeta, T_max, tau_min, v_max, v_min, Gamma = 10^(-9), H = 2, T_min = 0, tau_max = 0) {
  ## Initialization
  p_sim = numeric(length(price)) # Vector of simulated prices
  p_sim[1:(200 + 2)] = price[1:(200 + 2)] # The simulated and actual prices are the same for this interval so that the lagged price does not fall behind our window of data
  # Tracking
  N_c = numeric(length(price)) # Keep track of number of active chartists at each time point
  N_f = numeric(length(price)) # Keep track of number of active fundamentalists at each time point
  pi_c = numeric(length(price)) # Keep track of profit for chartists at each time point
  pi_f = numeric(length(price)) # Keep track of profit for fundamentalists at each time point
  # Data frame keeping track of all agents' strategies
  agents = data.frame(d = sample(size = N, x = d_min:d_max, replace = TRUE), # Vector of lags for chartists
                      Thresh = runif(N, T_min, T_max), # Vector of entry thresholds
                      tau = NA, c_value = NA, c_trend = NA, 
                      v_previous = runif(N, min = v_min + p_sim[201], max = v_max + p_sim[201]), # The initial value for each investor is exogenous
                      v_new = NA, x_previous_value = NA, x_previous_trend = NA, x_previous_actual = numeric(N), 
                      x_new_value = NA, x_new_trend = NA,  x_new_actual = numeric(N), 
                      m_value = NA, m_trend = NA, m_actual = numeric(N), 
                      profit_value = rep(0, times = N), profit_trend = rep(0, times = N), # No profit initially
                      full_profit_value = rep(0, times = N), full_profit_trend = rep(0, times = N), 
                      phi = rep(0.5, times = N), # phi = probability of being a chartist
                      chartist = as.logical(rbinom(N, size = 1, prob = 0.5)), # Randomly initialize strategies for each agent based on binomial distribution with probability 0.5
                      prev_day_chartist = NA, 
                      fundamentalist = NA, prev_day_fundamentalist = NA)
  agents$tau = runif(N, -colMins(t(cbind(agents$Thresh, -tau_min))), tau_max) # Used to ensure -tau < Thresh. Vector of exit thresholds
  agents$c_value = a*(agents$Thresh - agents$tau) # Vector of c's for all agents' value investing strategies
  agents$c_trend = agents$c_value # The ith agent has the same c for both his strategies (remains constant throughout)
  # Get mispricing for both fundamentalists and chartists
  agents$m_value = p_sim[201] - agents$v_previous # Mispricing that each fundamentlist believes there to be
  agents$m_trend = p_sim[201 - agents$d] - p_sim[201] # Mispricing that each chartist believes there to be
  # Set their initial positions as if the thresholds have been met
  agents$x_previous_trend = agents$c_trend*sign(p_sim[201] - p_sim[201 - agents$d]) # The initial positions for all agents' trend following strategies
  agents$x_previous_value = agents$c_value*sign(agents$v_previous - p_sim[201]) # The initial positions for all agents' value investing strategies
  # Check that thresholds have been met. Agents only take the positions above if they are outside one of the entry thresholds
  agents$x_previous_value[which(agents$m_value < agents$Thresh & agents$m_value > -agents$Thresh)] = 0 # No trade takes place
  agents$x_previous_trend[which(agents$m_trend < agents$Thresh & agents$m_trend > -agents$Thresh)] = 0 # No trade takes place
  agents$fundamentalist = !agents$chartist # If they're not currently a chartist then they're a fundamentalist
  agents$prev_day_chartist = agents$chartist
  agents$prev_day_fundamentalist = agents$fundamentalist
  past_profits_trend <- matrix(0, nrow = N, ncol = H)
  past_profits_value <- matrix(0, nrow = N, ncol = H)
  agents$x_previous_actual[which(agents$chartist)] = agents$x_previous_trend[which(agents$chartist)]
  agents$x_previous_actual[which(agents$fundamentalist)] = agents$x_previous_value[which(agents$fundamentalist)]
  # In what follows, both strategies for each agent is updated as normal but only their specific strategy for that time is used to calculate omega
  for (t in (200 + 3):(length(price))) {
    ## Chartists
    agents$x_new_trend = agents$c_trend*sign(p_sim[t - 1] - p_sim[t - 1 - agents$d]) # The positions chartists want to now move to. Update all agents' trend following strategies
    agents$m_trend = p_sim[t - 1 - agents$d] - p_sim[t - 1]
    ## Fundamentalists
    agents$v_new = agents$v_previous + rnorm(1, mean = mu_eta, sd = sigma_eta) # The fundamentalists' new value judgements. Update all agents' value investing strategies
    agents$m_value = p_sim[t - 1] - agents$v_new # Mispricing that each fundamentlist believes there to be
    agents$x_new_value = agents$c_value*sign(agents$v_new - p_sim[t - 1]) # The position they wish to now take (assuming thresholds are met for fundamentalists)
    ## State dependent value strategy
    # Fundamentalist was in a long position; will only trade if below tau
    agents$x_new_value[which(agents$x_previous_value > 0 & agents$m_value < -agents$tau)] =
      agents$x_previous_value[which(agents$x_previous_value > 0 & agents$m_value < -agents$tau)] # No trade takes place
    agents$x_new_value[which(agents$x_previous_value > 0 & agents$m_value < agents$Thresh & agents$m_value > -agents$tau)] = 0 # Exits long but does not enter short position as m doesn't meet entry threshold
    # Fundamentalist was in a short position; will only trade if above -tau
    agents$x_new_value[which(agents$x_previous_value < 0 & agents$m_value > agents$tau)] = 
      agents$x_previous_value[which(agents$x_previous_value < 0 & agents$m_value > agents$tau)] # No trade takes place
    agents$x_new_value[which(agents$x_previous_value < 0 & agents$m_value > -agents$Thresh & agents$m_value < agents$tau)] = 0 # Exits short but does not enter long position as m doesn't meet entry threshold
    # Fundamentalist was taking no position; will only trade if above Thresh or below -Thresh
    agents$x_new_value[which(agents$x_previous_value == 0 & agents$m_value < agents$Thresh & agents$m_value > -agents$Thresh)] = 0 # No trade takes place
    ## State dependent trend strategy
    # Chartist was in a long position; will only trade if below tau
    agents$x_new_trend[which(agents$x_previous_trend > 0 & agents$m_trend < -agents$tau)] =
      agents$x_previous_trend[which(agents$x_previous_trend > 0 & agents$m_trend < -agents$tau)] # No trade takes place
    agents$x_new_trend[which(agents$x_previous_trend > 0 & agents$m_trend < agents$Thresh & agents$m_trend > -agents$tau)] = 0 # Exits long but does not enter short position as m doesn't meet entry threshold
    # Chartist was in a short position; will only trade if above -tau
    agents$x_new_trend[which(agents$x_previous_trend < 0 & agents$m_trend > agents$tau)] = 
      agents$x_previous_trend[which(agents$x_previous_trend < 0 & agents$m_trend > agents$tau)] # No trade takes place
    agents$x_new_trend[which(agents$x_previous_trend < 0 & agents$m_trend > -agents$Thresh & agents$m_trend < agents$tau)] = 0 # Exits short but does not enter long position as m doesn't meet entry threshold
    # Chartist was taking no position; will only trade if above Thresh or below -Thresh
    agents$x_new_trend[which(agents$x_previous_trend == 0 & agents$m_trend < agents$Thresh & agents$m_trend > -agents$Thresh)] = 0 # No trade takes place
    ## Actual strategy
    agents$x_new_actual[which(agents$chartist)] = agents$x_new_trend[which(agents$chartist)]
    agents$m_actual[which(agents$chartist)] = agents$m_trend[which(agents$chartist)]
    agents$x_new_actual[which(agents$fundamentalist)] = agents$x_new_value[which(agents$fundamentalist)]
    agents$m_actual[which(agents$fundamentalist)] = agents$m_value[which(agents$fundamentalist)]
    # Agent was in a long position; will only trade if below tau
    agents$x_new_actual[which(agents$x_previous_actual > 0 & agents$m_actual < -agents$tau)] =
      agents$x_previous_actual[which(agents$x_previous_actual > 0 & agents$m_actual < -agents$tau)] # No trade takes place
    agents$x_new_actual[which(agents$x_previous_actual > 0 & agents$m_actual < agents$Thresh & agents$m_actual > -agents$tau)] = 0 # Exits long but does not enter short position as m doesn't meet entry threshold
    # Fundamentalist was in a short position; will only trade if above -tau
    agents$x_new_actual[which(agents$x_previous_actual < 0 & agents$m_actual > agents$tau)] = 
      agents$x_previous_actual[which(agents$x_previous_actual < 0 & agents$m_actual > agents$tau)] # No trade takes place
    agents$x_new_actual[which(agents$x_previous_actual < 0 & agents$m_actual > -agents$Thresh & agents$m_actual < agents$tau)] = 0 # Exits short but does not enter long position as m doesn't meet entry threshold
    # Fundamentalist was taking no position; will only trade if above Thresh or below -Thresh
    agents$x_new_actual[which(agents$x_previous_actual == 0 & agents$m_actual < agents$Thresh & agents$m_actual > -agents$Thresh)] = 0 # No trade takes place
    ## Fundamentalists and chartists
    omega = sum(agents$x_new_actual - agents$x_previous_actual) # The net order of all agents for the day. Select only those agents' specific strategies
    ## Market maker
    p_sim[t] = MarketImpact(lambda = lambda, p_t = p_sim[t - 1], omega = omega, zeta = rnorm(1, 0, sigma_zeta))
    ## Update adaptive probability, profit of each agent's strategy and which strategy each agent is adopting
    agents$profit_value = Profit(p_cur = p_sim[t], p_prev = p_sim[t - 1], x_prev = agents$x_previous_value)
    agents$profit_trend = Profit(p_cur = p_sim[t], p_prev = p_sim[t - 1], x_prev = agents$x_previous_trend)
    past_profits_trend <- cbind(agents$profit_trend, past_profits_trend[, 1:(H-1)])
    past_profits_value <- cbind(agents$profit_value, past_profits_value[, 1:(H-1)])
    agents$full_profit_trend <- rowSums(past_profits_trend)
    agents$full_profit_value <- rowSums(past_profits_value)
    agents$prev_day_chartist = agents$chartist
    agents$prev_day_fundamentalist = agents$fundamentalist
    ## Tracking
    N_c[t-1] = sum(agents$chartist)
    N_f[t-1] = sum(agents$fundamentalist)
    pi_c[t] = sum(agents$full_profit_trend)
    pi_f[t] = sum(agents$full_profit_value)
    agents$phi = (exp(agents$full_profit_trend/Gamma))/(exp(agents$full_profit_trend/Gamma) + exp(agents$full_profit_value/Gamma)) # Gamma is the intensity of switching parameter
    agents$chartist = as.logical(sapply(X = agents$phi, FUN = rbinom, n = 1, size = 1)) # Choose strategies based on each agent's probability of switching
    agents$fundamentalist = !agents$chartist
    ## Update positions
    agents$v_previous = agents$v_new # Prepare for the next loop by updating v_previous for all agents
    agents$x_previous_value = agents$x_new_value # Prepare for the next loop by updating x_previous_value for all agents
    agents$x_previous_trend = agents$x_new_trend # Prepare for the next loop by updating x_previous_trend for all agents
    agents$x_previous_actual = agents$x_new_actual # Prepare for the next loop by updating x_previous_actual for all agents
  }
  return(data.frame(p_sim = p_sim, N_c = N_c, N_f = N_f, pi_c = pi_c, pi_f = pi_f)) # Single simulated price path
}
####################################################################################################

### 7. Objective function
ObjectiveFunction = function(parameters) {
  I = 30 # Number of replications
  G = matrix(NA, nrow = k, ncol = I) # Initialise G
  for (i in 1:I) {
    set.seed(i) # Ensure there is a different seed for each replication
    p_sim = Simulate(price = price, N = ceiling(parameters[1]), lambda = parameters[2], a = parameters[3], d_max = ceiling(parameters[5]) + floor(parameters[4]), d_min = ceiling(parameters[5]), mu_eta = parameters[6], sigma_eta = parameters[7], sigma_zeta = parameters[8], T_max = parameters[9] + parameters[15], tau_min = parameters[10], v_max = parameters[11] + parameters[12], v_min = parameters[12], Gamma = parameters[13], H = ceiling(parameters[14]), T_min = parameters[15], tau_max = parameters[10] + parameters[16])$p_sim # Ceiling and floor functions act as integer constraints
    returns_sim = diff(p_sim[-c(1:201)])
    # Controlling errors from poor initialisations
    if(max(abs(na.omit(returns_sim))) > 10 | sum(is.na(returns_sim)) > 0) {
      mean_sim = NaN
      stdev_sim = NaN
      kurtosis_sim = NaN
      ks_sim = NaN
      hurst_sim = NaN
      gph_sim = NaN
      adf_sim = NaN
      garch_sim = NaN
      hill_sim = NaN
      return(NaN)
    } else {
      moments = Moments(returns_sim, returns)
      mean_sim = moments[1]
      stdev_sim = moments[2]
      kurtosis_sim = moments[3]
      ks_sim = moments[4]
      hurst_sim = moments[5]
      gph_sim = moments[6]
      adf_sim = moments[7]
      garch_sim = moments[8]
      hill_sim = moments[9]
    }
    G[, i] = c(mean_sim, stdev_sim, kurtosis_sim, ks_sim, hurst_sim, gph_sim, adf_sim, garch_sim, hill_sim) - Moments(returns, returns)
  }
  G = matrix(apply(G, 1, mean), nrow = k, ncol = 1)
  return(-t(G)%*%W%*%G) # Fitness
}
####################################################################################################

### 8. Genetic algorithm
## 95% confidence interval for optimal parameters
bootstrap_parameters_GA = rbind(matrix(NA, nrow = 5, ncol = 17), bootstrap_parameters_GA) # Extra column for the objective function value
for(i in 1:5) {
  GA_CI = ga(type = "real-valued", fitness = ObjectiveFunction, parallel = detectCores(), # Parameters are: N, lambda, a, d_max-dmin, d_min, mu_eta, sigma_eta, sigma_zeta, T_max-T_min, tau_min, v_max-vmin, v_min, Gamma, H, T_min, tau_max-tau_min
             lower = c(40, 0, 0, 0, 0.01, -0.01, 0, 0, 0, -1, 0, -0.5, 0, 0.01, 0, 0),
             upper = c(240, 15, 1, 100, 100, 0.01, 0.05, 0.05, 1, 0, 1, 0, 1, 100, 1, 1), maxiter = 30, seed = 127 + i)
  bootstrap_parameters_GA[i, ] = c(Decode1(GA_CI), summary(GA_CI)$fitness)
}
## Optimal parameters
args_GA = list(price, bootstrap_parameters_GA[11, 1], bootstrap_parameters_GA[11, 2], bootstrap_parameters_GA[11, 3], bootstrap_parameters_GA[11, 4], bootstrap_parameters_GA[11, 5], bootstrap_parameters_GA[11, 6], bootstrap_parameters_GA[11, 7], bootstrap_parameters_GA[11, 8], bootstrap_parameters_GA[11, 9], bootstrap_parameters_GA[11, 10], bootstrap_parameters_GA[11, 11], bootstrap_parameters_GA[11, 12], bootstrap_parameters_GA[11, 13], bootstrap_parameters_GA[11, 14], bootstrap_parameters_GA[11, 15], bootstrap_parameters_GA[11, 16])
## Price path confidence intervals
CI_prices_GA = foreach(i = 1:1000, .combine = "cbind", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
  set.seed(i)
  do.call(Simulate, args = args_GA)$p_sim[-c(1:201)]
}
sorted_prices_GA = t(apply(CI_prices_GA, 1, sort))
## 95% confidence intervals for moments of simulated data
bootstrap_moments_GA = matrix(NA, nrow = 1000, ncol = 9)
cl = makeCluster(detectCores())
registerDoParallel(cl)
bootstrap_moments_GA = foreach(i = 1:1000, .combine = "rbind", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
  set.seed(i)
  Moments(diff(do.call(Simulate, args = args_GA)$p_sim[-c(1:201)]), returns)
}
stopCluster(cl)
set.seed(12)
simulation_GA = do.call(Simulate, args = args_GA)
simulation_price_GA = simulation_GA$p_sim[-c(1:201)] # Simulate using optimised parameters
simulation_returns_GA = diff(simulation_price_GA)
N_c_GA = simulation_GA$N_c[-c(1:201, length(price))]
N_f_GA = simulation_GA$N_f[-c(1:201, length(price))]
pi_c_GA = simulation_GA$pi_c[-c(1:202)]
pi_f_GA = simulation_GA$pi_f[-c(1:202)]
####################################################################################################

### 9. Nelder-Mead with threshold accepting
rssimple = function(x, d = 50) {
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

Moments = function(returns1, returns2) {
  mean_e = mean(returns1)
  stdev_e = sd(returns1)
  kurtosis_e = kurtosis(returns1)
  ks_e = ks.test(x = returns1, y = returns2)$statistic # Kolmogorov-Smirnov statistic is zero since we are dealing with actual data. Represents difference between cummulative emperical distribution functions
  hurst_e = rssimple(returns1)
  gph_e = fdGPH(abs(returns1))$d # GPH estimator representing long-range dependence
  adf_e = adf.test(returns1, k = 0)$statistic # ADF statistic representing random walk property of returns
  garch_e = sum(garch(returns1, trace = F)$coef[2:3]) # GARCH paramaters representing short-range dependence
  hill_e = mean(hill(returns1)$hill[(0.05 * length(returns1)):(0.1 * length(returns1))]) # Hill estimator
  return(c(mean_e, stdev_e, kurtosis_e, ks_e, hurst_e, gph_e, adf_e, garch_e, hill_e))
}

ObjectiveFunction = function(parameters) {
  I = parameters[17] # Number of Monte Carlo replications
  G = matrix(NA, nrow = 9, ncol = I) # Initialise G
  for (i in 1:I) {
    p_sim = Simulate(price = price, N = ceiling(parameters[1]), lambda = parameters[2], a = parameters[3], d_max = ceiling(parameters[5]) + floor(parameters[4]), d_min = ceiling(parameters[5]), mu_eta = parameters[6], sigma_eta = parameters[7], sigma_zeta = parameters[8], T_max = parameters[9] + parameters[15], tau_min = parameters[10], v_max = parameters[11] + parameters[12], v_min = parameters[12], Gamma = parameters[13], H = ceiling(parameters[14]), T_min = parameters[15], tau_max = parameters[10] + parameters[16])$p_sim # Ceiling and floor functions act as integer constraints
    returns_sim = diff(p_sim[-c(1:201)])
    if(max(abs(na.omit(returns_sim))) > 10 | sum(is.na(returns_sim)) > 0) {
      return(-100000)
    } else {
      moments = Moments(returns_sim, returns)
      mean_sim = moments[1]
      stdev_sim = moments[2]
      kurtosis_sim = moments[3]
      ks_sim = moments[4]
      hurst_sim = moments[5]
      gph_sim = moments[6]
      adf_sim = moments[7]
      garch_sim = moments[8]
      hill_sim = moments[9]
    }
    G[, i] = c(mean_sim, stdev_sim, kurtosis_sim, ks_sim, hurst_sim, gph_sim, adf_sim, garch_sim, hill_sim) - Moments(returns, returns)
  }
  G = apply(G, 1, mean)
  return(-t(G)%*%W%*%G)
}

ObjectiveFunction2 = function(parameters) {
  I = parameters[17] # Number of Monte Carlo replications
  G <- foreach(i = 1:I, .combine = "cbind", .export = ls(globalenv()), .packages = list.of.packages) %dorng% {
                 p_sim = Simulate(price = price, N = ceiling(parameters[1]), lambda = parameters[2], a = parameters[3], d_max = ceiling(parameters[5]) + floor(parameters[4]), d_min = ceiling(parameters[5]), mu_eta = parameters[6], sigma_eta = parameters[7], sigma_zeta = parameters[8], T_max = parameters[9] + parameters[15], tau_min = parameters[10], v_max = parameters[11] + parameters[12], v_min = parameters[12], Gamma = parameters[13], H = ceiling(parameters[14]), T_min = parameters[15], tau_max = parameters[10] + parameters[16])$p_sim # Ceiling and floor functions act as integer constraints
                 returns_sim = diff(p_sim[-c(1:201)])
                 if(max(abs(na.omit(returns_sim))) > 10 | sum(is.na(returns_sim)) > 0) {
                   return(rep(NA, 9))
                 } else {
                   moments = Moments(returns_sim, returns)
                   mean_sim = moments[1]
                   stdev_sim = moments[2]
                   kurtosis_sim = moments[3]
                   ks_sim = moments[4]
                   hurst_sim = moments[5]
                   gph_sim = moments[6]
                   adf_sim = moments[7]
                   garch_sim = moments[8]
                   hill_sim = moments[9]
                 }
                 c(mean_sim, stdev_sim, kurtosis_sim, ks_sim, hurst_sim, gph_sim, adf_sim, garch_sim, hill_sim) - Moments(returns, returns)
               }
  if (sum(is.na(G)) > 0) {
    return(-100000)
  }
  G = apply(G, 1, mean)
  return(-t(G)%*%W%*%G)
}

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
      temp_vertices[1, ] <- pmax(30, temp_vertices[1, ])
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
      temp_vertices[13, ] <- pmax(0.00000001, temp_vertices[13, ])
      temp_vertices[14, ] <- pmax(1, temp_vertices[14, ])
      temp_vertices[14, ] <- pmin(200, temp_vertices[14, ])
      temp_vertices[16, ] <- pmax(0, temp_vertices[16, ])
      # Get their fitnesses
      temp_fitnesses <- foreach(k = 1:(n+1), .combine = "c", .export = ls(globalenv()), .packages = list.of.packages) %dorng% {
        ObjectiveFunction(c(temp_vertices[, k], 8 + 4 * i))
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

Calibrate_NMTA = function(maxiter, I) { # Based off of Platt. 'I' is the number of replications
  # Get starting vertices and fitnesses
  initials <- foreach(k = 1:17, .combine = "cbind", .export = ls(globalenv()), .packages = list.of.packages) %dorng% {
    vertex <- c(runif(1, 40, 240), runif(1, 0, 15), runif(1, 0, 1.5), runif(1, 0, 100), runif(1, 0.01, 90), runif(1, -0.01, 0.01), runif(1, 0, 0.05), runif(1, 0, 0.05), runif(1, 0, 1), runif(1, -1, 0), runif(1, 0, 1), runif(1, -0.5, 0), runif(1, 0, 1), runif(1, 0, 100), runif(1, 0, 1), runif(1, 0, 1))
    fitness <- ObjectiveFunction(c(vertex, I))
    while(fitness == -100000) { # Ensure that simulation does not result in crazy behaviour
      vertex <- c(runif(1, 40, 240), runif(1, 0, 15), runif(1, 0, 1.5), runif(1, 0, 100), runif(1, 0.01, 90), runif(1, -0.01, 0.01), runif(1, 0, 0.05), runif(1, 0, 0.05), runif(1, 0, 1), runif(1, -1, 0), runif(1, 0, 1), runif(1, -0.5, 0), runif(1, 0, 1), runif(1, 0, 100), runif(1, 0, 1), runif(1, 0, 1))
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
  names(progress) <- c("Iteration", "fitness", "N", "lambda", "a", "d_max - d_min", "d_min", "mu_eta", "sigma_eta", "sigma_zeta", "T_max - T_min", "tau_min", "v_max - v_min", "v_min", "Gamma", "H", "T_min", "tau_max - tau_min")
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
      vertices <- new_parameters[1:(16), ]
      fitnesses <- new_parameters[17, ]
    }
    else {
      # Reflect the worst point to the other side of the average of each of the other parameters
      centroid <- rowMeans(vertices[, -17])
      reflected_point <- centroid + alphaNM * (centroid - vertices[, 17])
      # Restrict shifts to possible values only
      reflected_point[1] <- max(30, reflected_point[1])
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
      reflected_point[13] <- max(0.00000001, reflected_point[13])
      reflected_point[14] <- max(1, reflected_point[14])
      reflected_point[14] <- min(200, reflected_point[14])
      reflected_point[16] <- max(0, reflected_point[16])
      ref_fitness <- ObjectiveFunction2(c(reflected_point, I))
      # Check if new solution is in between best and worst fitnesses
      if(ref_fitness < fitnesses[1] & ref_fitness > fitnesses[16]) {
        vertices[, 17] <- reflected_point
        fitnesses[17] <- ref_fitness
      }
      # Check if new parameters are new best solution; if so try improve further by moving more in that direction
      else if(ref_fitness > fitnesses[1]) {
        # Expand
        expanded_point <- centroid + gammaNM * (reflected_point - centroid)
        # Restrict shifts to possible values only
        expanded_point[1] <- max(30, expanded_point[1])
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
        expanded_point[13] <- max(0.00000001, expanded_point[13])
        expanded_point[14] <- max(1, expanded_point[14])
        expanded_point[14] <- min(200, expanded_point[14])
        expanded_point[16] <- max(0, expanded_point[16])
        exp_fitness <- ObjectiveFunction2(c(expanded_point, I))
        # Check if the solution is even better now
        if(exp_fitness > ref_fitness) {
          vertices[, 17] <- expanded_point
          fitnesses[17] <- exp_fitness
        } else { # If not then just use the original improvement
          vertices[, 17] <- reflected_point
          fitnesses[17] <- ref_fitness
        }
      }
      # Check if new parameters are still the worst. If so, try contracting
      else if(ref_fitness < fitnesses[16]) {
        contracted_point <- centroid + rhoNM * (vertices[, 17] - centroid)
        con_fitness <- ObjectiveFunction2(c(contracted_point, I))
        # Check if contracting managed to actually get a better result than the original worst parameter
        if(con_fitness > fitnesses[17]) {
          vertices[, 17] <- contracted_point
          fitnesses[17] <- con_fitness
        } else {
          # If none of the above helped, shrink everything towards best solution
          vertices <- vertices[, 1] + sigmaNM * (vertices - vertices[, 1])
          fitnesses <- c(fitnesses[1], foreach(k = 2:17, .combine = "c", .export = ls(globalenv()), .packages = list.of.packages) %dorng% {
            ObjectiveFunction(c(vertices[, k], I))
          }
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
    names(progress) <- c("Iteration", "fitness", "N", "lambda", "a", "d_max - d_min", "d_min", "mu_eta", "sigma_eta", "sigma_zeta", "T_max - T_min", "tau_min", "v_max - v_min", "v_min", "Gamma", "H", "T_min", "tau_max - tau_min")
    print(progress)
  }
  return(c(best_ever_iteration, best_ever_fitness, best_ever_vertex))
}

## 95% confidence interval for optimal parameters
bootstrap_parameters_NMTA = rbind(matrix(NA, nrow = 9, ncol = 17), bootstrap_parameters_NMTA)
for(i in 1:9) {
  cl = makeCluster(detectCores())
  registerDoParallel(cl)
  set.seed(i + 2)
  NMTA_CI = Calibrate_NMTA(50, 30)
  stopCluster(cl)
  bootstrap_parameters_NMTA[i,] = c(ceiling(NMTA_CI[3]), NMTA_CI[4], NMTA_CI[5], ceiling(NMTA_CI[7]) + floor(NMTA_CI[6]), ceiling(NMTA_CI[7]), NMTA_CI[8], NMTA_CI[9], NMTA_CI[10], NMTA_CI[11] + NMTA_CI[17], NMTA_CI[12], NMTA_CI[13] + NMTA_CI[14], NMTA_CI[14], NMTA_CI[15], ceiling(NMTA_CI[16]), NMTA_CI[17], NMTA_CI[12] + NMTA_CI[18], NMTA_CI[2])
}
write.table(cbind(c("N", "lambda", "a", "d_max", "d_min", "mu_eta", "sigma_eta", "sigma_zeta", "T_max", "tau_min", "v_max", "v_min", "Gamma", "H", "T_min", "tau_max", "Fitness"), round(bootstrap_parameters_GA[11, ], 5), round(colMeans(bootstrap_parameters_GA) - qt(0.975, df = nrow(bootstrap_parameters_GA) - 1)*colStdevs(bootstrap_parameters_GA)/sqrt(nrow(bootstrap_parameters_GA)), 5), round(colMeans(bootstrap_parameters_GA) + qt(0.975, df = nrow(bootstrap_parameters_GA) - 1)*colStdevs(bootstrap_parameters_GA)/sqrt(nrow(bootstrap_parameters_GA)), 5), round(bootstrap_parameters_NMTA[11, ], 5), round(colMeans(bootstrap_parameters_NMTA) - qt(0.975, df = nrow(bootstrap_parameters_NMTA) - 1)*colStdevs(bootstrap_parameters_NMTA)/sqrt(nrow(bootstrap_parameters_NMTA)), 5), round(colMeans(bootstrap_parameters_NMTA) + qt(0.975, df = nrow(bootstrap_parameters_NMTA) - 1)*colStdevs(bootstrap_parameters_NMTA)/sqrt(nrow(bootstrap_parameters_NMTA)), 5)), file = "Parameters.txt", append = FALSE, quote = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Parameter", "GA", "Lower_GA", "Upper_GA", "NMTA", "Lower_NMTA", "Upper_NMTA"))
## Optimal parameters
args_NMTA = list(price, bootstrap_parameters_NMTA[11, 1], bootstrap_parameters_NMTA[11, 2], bootstrap_parameters_NMTA[11, 3], bootstrap_parameters_NMTA[11, 4], bootstrap_parameters_NMTA[11, 5], bootstrap_parameters_NMTA[11, 6], bootstrap_parameters_NMTA[11, 7], bootstrap_parameters_NMTA[11, 8], bootstrap_parameters_NMTA[11, 9], bootstrap_parameters_NMTA[11, 10], bootstrap_parameters_NMTA[11, 11], bootstrap_parameters_NMTA[11, 12], bootstrap_parameters_NMTA[11, 13], bootstrap_parameters_NMTA[11, 14], bootstrap_parameters_NMTA[11, 15], bootstrap_parameters_NMTA[11, 16])
## Price path confidence intervals
CI_prices_NMTA = foreach(i = 1:1000, .combine = "cbind", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
  set.seed(i)
  do.call(Simulate, args = args_NMTA)$p_sim[-c(1:201)]
}
sorted_prices_NMTA = t(apply(CI_prices_NMTA, 1, sort))
## 95% confidence intervals for moments of simulated data
bootstrap_moments_NMTA = matrix(NA, nrow = 1000, ncol = 9)
cl = makeCluster(detectCores())
registerDoParallel(cl)
bootstrap_moments_NMTA = foreach(i = 1:1000, .combine = "rbind", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
  set.seed(i)
  Moments(diff(do.call(Simulate, args = args_NMTA)$p_sim[-c(1:201)]), returns)
}
stopCluster(cl)
set.seed(12)
simulation_NMTA = do.call(Simulate, args = args_NMTA)
simulation_price_NMTA = simulation_NMTA$p_sim[-c(1:201)] # Simulate using optimised parameters
simulation_returns_NMTA = diff(simulation_price_NMTA)
N_c_NMTA = simulation_NMTA$N_c[-c(1:201, length(price))]
N_f_NMTA = simulation_NMTA$N_f[-c(1:201, length(price))]
pi_c_NMTA = simulation_NMTA$pi_c[-c(1:202)]
pi_f_NMTA = simulation_NMTA$pi_f[-c(1:202)]
write.table(cbind(c("mean", "stdev", "kurtosis", "ks", "hurst", "gph", "adf", "garch", "hill"), round(Moments(simulation_returns_GA, returns), 5), round(apply(bootstrap_moments_GA, 2, sort)[25,], 5), round(apply(bootstrap_moments_GA, 2, sort)[975,], 5), round(Moments(simulation_returns_NMTA, returns), 5), round(apply(bootstrap_moments_NMTA, 2, sort)[25,], 5), round(apply(bootstrap_moments_NMTA, 2, sort)[975,], 5), round(Moments(returns, returns), 5)), file = "Moments.txt", append = FALSE, quote = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = c("Moment", "GA", "Lower_GA", "Upper_GA", "NMTA", "Lower_NMTA", "Upper_NMTA", "Emperical"))
####################################################################################################

### 10. Results and replication of stylized facts of financial markets
## Paths
# Log Price
pdf("Price paths.pdf", width = 10, height = 5)
data.frame(simulation_price = c(price[-c(1:201)], simulation_price_GA, simulation_price_NMTA), type = rep(c("Actual", "GA Simulation", "NMTA Simulation"), each = length(simulation_price_GA)), date = rep(date, times = 3)) %>%
ggplot(data = ., mapping = aes(x = date, y = simulation_price, colour = type)) +
  geom_line() +
  scale_colour_manual(name = "Legend", values = c(Actual = "black", `GA Simulation` = "blue", `NMTA Simulation` = "red")) +
  labs(x = "Date", y = "Log Price")
dev.off()
# Log Returns
pdf("Return paths.pdf", width = 10, height = 5)
data.frame(returns = c(returns, simulation_returns_GA, simulation_returns_NMTA), date = rep(date[-1], times = 3), type = rep(c("Actual", "GA Simulation", "NMTA Simulation"), each = length(simulation_returns_GA))) %>%
ggplot(data = ., mapping = aes(x = date, y = returns, color = type)) +
  geom_line() +
  scale_colour_manual(values = c(Actual = "black", `GA Simulation` = "blue", `NMTA Simulation` = "red")) +
  facet_wrap(~ type, nrow = 3) +
  theme(legend.position = "none") +
  labs(x = "Date", y = "Log Returns")
dev.off()

## Normal probability plots of log returns
pdf("Normal probability plots.pdf", width = 10, height = 5)
data.frame(returns = c(scale(returns), scale(simulation_returns_GA), scale(simulation_returns_NMTA)), type = rep(c("Actual", "GA Simulation", "NMTA Simulation"), each = length(simulation_returns_GA))) %>%
ggplot(data = ., mapping = aes(sample = returns)) +
  stat_qq_point() +
  stat_qq_line(col = "black") +
  stat_qq_band(fill = rep(c("black", "blue", "red"), each = length(returns)), alpha = 0.4) + # alpha = transparency
  facet_wrap(~ type, ncol = 3) + # Returns should be log normal
  labs(x = "Theoretical", y = "Sample")
dev.off()

## Autocovariances
# Returns
pdf("Autocorrelation of log returns.pdf", width = 10, height = 5)
cbind(rbind(with(acf(scale(returns), type = "covariance", lag.max = 50, plot = FALSE), data.frame(lag, acf)), with(acf(scale(simulation_returns_GA), type = "covariance", lag.max = 50, plot = FALSE), data.frame(lag, acf)), with(acf(scale(simulation_returns_NMTA), type = "covariance", lag.max = 50, plot = FALSE), data.frame(lag, acf))), type = rep(c("Actual", "GA Simulation", "NMTA Simulation"), each = 51)) %>%
ggplot(data = ., mapping = aes(x = lag, y = acf)) +
  geom_point(mapping = aes(colour = type)) +
  geom_segment(mapping = aes(x = lag, y = 0, xend = lag, yend = acf)) +
  geom_hline(yintercept = c(-1, 1)*qnorm((1 + 0.95)/2)/sqrt(length(returns)), lty = 2, col = "blue") +
  scale_colour_manual(values = c(Actual = "black", `GA Simulation` = "blue", `NMTA Simulation` = "red")) +
  facet_wrap(~ type, ncol = 3) +
  theme(legend.position = "none") +
  labs(x = "Lag", y = "Autocorrelation")
dev.off()
# Absolute Returns
pdf("Autocorrelation of absolute log returns.pdf", width = 10, height = 5)
cbind(rbind(with(acf(scale(abs(returns)), type = "covariance", lag.max = 50, plot = FALSE), data.frame(lag, acf)), with(acf(scale(abs(simulation_returns_GA)), type = "covariance", lag.max = 50, plot = FALSE), data.frame(lag, acf)), with(acf(scale(abs(simulation_returns_NMTA)), type = "covariance", lag.max = 50, plot = FALSE), data.frame(lag, acf))), type = rep(c("Actual", "GA Simulation", "NMTA Simulation"), each = 51)) %>%
  ggplot(data = ., mapping = aes(x = lag, y = acf)) +
  geom_point(mapping = aes(colour = type)) +
  geom_segment(mapping = aes(x = lag, y = 0, xend = lag, yend = acf)) +
  geom_hline(yintercept = c(-1, 1)*qnorm((1 + 0.95)/2)/sqrt(length(returns)), lty = 2, color = "blue") +
  scale_colour_manual(values = c(Actual = "black", `GA Simulation` = "blue", `NMTA Simulation` = "red")) +
  facet_wrap(~ type, ncol = 3) +
  theme(legend.position = "none") +
  labs(x = "Lag", y = "Autocorrelation")
dev.off()

## Agents
# Number
pdf("Number of chartists-fundamentalists.pdf", width = 10, height = 5) # TTR
ggplot(data = data.frame(N = c(N_c_GA, N_f_GA, N_c_NMTA, N_f_NMTA), N_smooth = c(SMA(N_c_GA, n = 40), SMA(N_f_GA, n = 40), SMA(N_c_NMTA, n = 40), SMA(N_f_NMTA, n = 40)), date = rep(date[-length(date)], times = 4), type = rep(c("GA Chartists", "GA Fundamentalists", "NMTA Chartists", "NMTA Fundamentalists"), each = length(N_c_GA))), mapping = aes(x = date, y = N, color = type)) +
  geom_line() +
  geom_line(mapping = aes(x = date, y = N_smooth), col = "white") +
  scale_colour_manual(values = c(`GA Chartists` = "blue", `GA Fundamentalists` = "blue", `NMTA Chartists` = "red", `NMTA Fundamentalists` = "red")) +
  facet_wrap(~ type, ncol = 2, nrow = 2) +
  theme(legend.position = "none") +
  labs(x = "Date", y = "Number of chartists/fundamentalists")
dev.off()
# Profit
pdf("Profit of chartists-fundamentalists.pdf", width = 10, height = 5)
ggplot(data = data.frame(pi = c(pi_c_GA, pi_f_GA, pi_c_NMTA, pi_f_NMTA), type = rep(c("GA Simulation", "NMTA Simulation"), each = 2*length(pi_c_GA)), Agent = rep(c("Chartists", "Fundamentalists", "Chartists", "Fundamentalists"), each = length(pi_c_GA)), date = rep(date[-1], times = 4)), mapping = aes(x = date, y = pi, colour = Agent)) +
  geom_line() +
  scale_colour_manual(name = "Legend", values = c(`Chartists` = "green", `Fundamentalists` = "darkgreen")) +
  facet_wrap(~ type, ncol = 2) +
  labs(x = "Date", y = "Profit")
dev.off()

## 3D plots
params = c(summary(GA)$solution) # Parameters are: N, lambda, a, d_max-d_min, d_min, mu_eta, sigma_eta, sigma_zeta, T_max-T_min, tau_min, v_max-v_min, v_min, Gamma, H, T_min, tau_max-tau_min
# a, N
a = seq(from = 0.00001, to = 0.3, length.out = 48)
N = seq(from = 70, to = 110, length.out = 48)
pb <- txtProgressBar(0, length(a), style = 3)
fitness_1 = matrix(NA, nrow = length(N), ncol = length(a))
cl = makeCluster(detectCores())
registerDoParallel(cl)
for (i in 1:length(N)) {
  params[1] = N[i]
  one_row = foreach(j = 1:length(a), .combine = "c", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
    params[3] = a[j]
    ObjectiveFunction(params)
  }
  fitness_1[i, ] = one_row
  setTxtProgressBar(pb, i)
}
stopCluster(cl)
a = matrix(a, nrow = length(a), ncol = length(a))
N = matrix(N, nrow = length(N), ncol = length(N), byrow = TRUE)
pdf("3D - a, N.pdf")
persp3D(x = a, y = N, z = fitness_1, colvar = fitness_1, colkey = FALSE, facets = TRUE, bty = "g", ticktype = "detailed", phi = 20, xlab = "a", ylab = "N", zlab = "Fitness")
dev.off()
# d_min, d_max - d_min
d_min = seq(from = 20, to = 45, length.out = 48)
d_max = seq(from = 40, to = 70, length.out = 48)
pb <- txtProgressBar(0, length(d_min), style = 3)
fitness_2 = matrix(NA, nrow = length(d_max), ncol = length(d_min))
cl = makeCluster(detectCores())
registerDoParallel(cl)
for (i in 1:length(d_max)) {
  params[4] = d_max[i]
  one_row = foreach(j = 1:length(d_min), .combine = "c", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
    params[5] = d_min[j]
    ObjectiveFunction(params)
  }
  fitness_2[i, ] = one_row
  setTxtProgressBar(pb, i)
}
stopCluster(cl)
d_min = matrix(d_min, nrow = length(d_min), ncol = length(d_min))
d_max = matrix(d_max, nrow = length(d_max), ncol = length(d_max), byrow = TRUE)
pdf("3D - d_min, d_max - d_min.pdf")
persp3D(x = d_min, y = d_max, z = fitness_2, colvar = fitness_2, colkey = FALSE, facets = TRUE, bty = "g", ticktype = "detailed", phi = 20, xlab = expression(d[min]), ylab = expression(d[max] - d[min]), zlab = "Fitness")
dev.off()
# sigma_zeta, lambda
sigma_zeta = seq(from = 0.00001, to = 0.15, length.out = 24)
lambda = seq(from = 5, to = 15, length.out = 24)
pb = txtProgressBar(0, length(sigma_zeta), style = 3)
fitness_3 = matrix(NA, nrow = length(lambda), ncol = length(sigma_zeta))
cl = makeCluster(detectCores())
registerDoParallel(cl)
for (i in 1:length(lambda)) {
  params[2] = lambda[i]
  one_row = foreach(j = 1:length(sigma_zeta), .combine = "c", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
    params[8] = sigma_zeta[j]
    ObjectiveFunction(params)
  }
  fitness_3[i, ] = one_row
  setTxtProgressBar(pb, i)
}
stopCluster(cl)
sigma_zeta = matrix(sigma_zeta, nrow = length(sigma_zeta), ncol = length(sigma_zeta))
lambda = matrix(lambda, nrow = length(lambda), ncol = length(lambda), byrow = TRUE)
pdf("3D - sigma_zeta, lambda.pdf")
persp3D(x = sigma_zeta, y = lambda, z = fitness_3, colvar = fitness_3, colkey = FALSE, facets = TRUE, bty = "g", ticktype = "detailed", phi = 20, xlab = "sigma_{zeta}", ylab = "lambda", zlab = "Fitness")
dev.off()
# sigma_eta, mu_eta
sigma_eta = seq(from = 0.0001, to = 0.15, length.out = 24)
mu_eta = seq(from = -0.15, to = 0.15, length.out = 24)
pb = txtProgressBar(0, length(sigma_eta), style = 3)
fitness_4 = matrix(NA, nrow = length(mu_eta), ncol = length(sigma_eta))
cl = makeCluster(detectCores())
registerDoParallel(cl)
for (i in 1:length(mu_eta)) {
  params[6] = mu_eta[i]
  one_row = foreach(j = 1:length(sigma_eta), .combine = "c", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
    params[7] = sigma_eta[j]
    ObjectiveFunction(params)
  }
  fitness_4[i, ] = one_row
  setTxtProgressBar(pb, i)
}
stopCluster(cl)
sigma_eta = matrix(sigma_eta, nrow = length(sigma_eta), ncol = length(sigma_eta))
mu_eta = matrix(mu_eta, nrow = length(mu_eta), ncol = length(mu_eta), byrow = TRUE)
pdf("3D - sigma_eta, mu_eta.pdf")
persp3D(x = sigma_eta, y = mu_eta, z = fitness_4, colvar = fitness_4, colkey = FALSE, facets = TRUE, bty = "g", ticktype = "detailed", phi = 20, xlab = "sigma_{eta}", ylab = "mu_{eta}", zlab = "Fitness")
dev.off()
# T_min, T_max - T_min
T_min = seq(from = 0.0001, to = 1.5, length.out = 24)
T_max = seq(from = 0.0001, to = 1.5, length.out = 24)
pb = txtProgressBar(0, length(T_min), style = 3)
fitness_5 = matrix(NA, nrow = length(T_max), ncol = length(T_min))
cl = makeCluster(detectCores())
registerDoParallel(cl)
for (i in 1:length(T_max)) {
  params[9] = T_max[i]
  one_row = foreach(j = 1:length(T_min), .combine = "c", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
    params[15] = T_min[j]
    ObjectiveFunction(params)
  }
  fitness_5[i, ] = one_row
  setTxtProgressBar(pb, i)
}
stopCluster(cl)
T_min = matrix(T_min, nrow = length(T_min), ncol = length(T_min))
T_max = matrix(T_max, nrow = length(T_max), ncol = length(T_max), byrow = TRUE)
pdf("3D - T_min, T_max - T_min.pdf")
persp3D(x = T_min, y = T_max, z = fitness_5, colvar = fitness_5, colkey = FALSE, facets = TRUE, bty = "g", ticktype = "detailed", phi = 20, xlab = "T_{min}", ylab = "T_{max} - T_{min}", zlab = "Fitness")
dev.off()
# v_min, v_max - v_min
v_min = rev(seq(from = -1.5, to = 1.5, length.out = 24))
v_max = seq(from = 0.0001, to = 2, length.out = 24)
pb = txtProgressBar(0, length(v_min), style = 3)
fitness_6 = matrix(NA, nrow = length(v_max), ncol = length(v_min))
cl = makeCluster(detectCores())
registerDoParallel(cl)
for (i in 1:length(v_max)) {
  params[11] = v_max[i]
  one_row = foreach(j = 1:length(v_min), .combine = "c", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
    params[12] = v_min[j]
    ObjectiveFunction(params)
  }
  fitness_6[i, ] = one_row
  setTxtProgressBar(pb, i)
}
stopCluster(cl)
v_min = matrix(v_min, nrow = length(v_min), ncol = length(v_min))
v_max = matrix(v_max, nrow = length(v_max), ncol = length(v_max), byrow = TRUE)
pdf("3D - v_min, v_max - v_min.pdf")
persp3D(x = v_min, y = v_max, z = fitness_6, colvar = fitness_6, colkey = FALSE, facets = TRUE, bty = "g", ticktype = "detailed", phi = 20, xlab = "v_{min}", ylab = "v_{max} - v_{min}", zlab = "Fitness")
dev.off()
# H, Gamma
H = seq(from = 30, to = 60, length.out = 24)
Gamma = seq(from = -1, to = 2, length.out = 24)
pb = txtProgressBar(0, length(H), style = 3)
fitness_7 = matrix(NA, nrow = length(Gamma), ncol = length(H))
cl = makeCluster(detectCores())
registerDoParallel(cl)
for (i in 1:length(Gamma)) {
  params[13] = Gamma[i]
  one_row = foreach(j = 1:length(H), .combine = "c", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
    params[14] = H[j]
    ObjectiveFunction(params)
  }
  fitness_7[i, ] = one_row
  setTxtProgressBar(pb, i)
}
stopCluster(cl)
H = matrix(H, nrow = length(H), ncol = length(H))
Gamma = matrix(Gamma, nrow = length(Gamma), ncol = length(Gamma), byrow = TRUE)
pdf("3D - H, Gamma.pdf")
persp3D(x = H, y = Gamma, z = fitness_7, colvar = fitness_7, colkey = FALSE, facets = TRUE, bty = "g", ticktype = "detailed", phi = 20, xlab = "H", ylab = "Gamma", zlab = "Fitness")
dev.off()
# tau_min, tau_max - tau_min
tau_min = seq(from = -1.5, to = 1.5, length.out = 24)
tau_max = seq(from = 0.0001, to = 2, length.out = 24)
pb = txtProgressBar(0, length(tau_min), style = 3)
fitness_8 = matrix(NA, nrow = length(tau_max), ncol = length(tau_min))
cl = makeCluster(detectCores())
registerDoParallel(cl)
for (i in 1:length(tau_max)) {
  params[16] = tau_max[i]
  one_row = foreach(j = 1:length(tau_min), .combine = "c", .packages = c("pracma", "timeDate", "tseries", "fracdiff", "extremefit", "timeSeries")) %dopar% {
    params[10] = tau_min[j]
    ObjectiveFunction(params)
  }
  fitness_8[i, ] = one_row
  setTxtProgressBar(pb, i)
}
stopCluster(cl)
tau_min = matrix(tau_min, nrow = length(tau_min), ncol = length(tau_min))
tau_max = matrix(tau_max, nrow = length(tau_max), ncol = length(tau_max), byrow = TRUE)
pdf("3D - tau_min, tau_max - tau_min.pdf")
persp3D(x = tau_min, y = tau_max, z = fitness_8, colvar = fitness_8, colkey = FALSE, facets = TRUE, bty = "g", ticktype = "detailed", phi = 20, xlab = "tau_{min}", ylab = "tau_{max} - tau_{min}", zlab = "Fitness")
dev.off()
####################################################################################################

### 11. Save workspace
save.image("Workspace1.RData")
####################################################################################################
