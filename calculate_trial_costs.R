#========================================================
# Trial cost calculation function
#========================================================

#' Calculate trial costs for different sample sizes
#' @param sample_sizes Vector of sample sizes to calculate costs for
#' @param n_sims Number of simulations for uncertainty quantification
#' @return Dataframe with sample_size, mean_cost, lower_ci, upper_ci

calculate_trial_costs <- function(sample_sizes, n_sims = 100000) {


  
  # Model parameters (point estimates) by van Asselt et al.(2018) https://doi.org/10.1007/s40273-017-0572-7
  intercept <- 359213.26
  submission_2014 <- -148012.07
  hta <- -38625.93
  study_duration <- 4500.73
  duration_months <- 60
  sample_size_coef <- 27.12
  study_locations <- 3278.97
  n_locations <- 1
  
  # Standard errors for each component
  se_intercept <- 209241.09
  se_submission_2014 <- 204475.64
  se_hta <- 68825.53
  se_study_duration <- 3046.40
  se_study_locations <- 3816.42
  se_sample_size_coef <- 23.53
  
  # Test parameters
  snlb_cost <- 647
  snlb_sd <- 96.682
  ly75_cost <- 500
  ly75_sd <- 74.818
  
  # CPI conversion factor (2014 to 2024)
  cpi_2014_to_2024 <- 1.310965795
  
  # Calculate point estimate for a given sample size
  calculate_point_estimate <- function(n) {
    # Fixed study costs (not dependent on n)
    study_costs_2014_base <- intercept + submission_2014 + hta + 
      (study_duration * duration_months) + study_locations
    
    # Add sample size component
    total_study_2014 <- study_costs_2014_base + (sample_size_coef * n)
    
    # Convert to 2024 euros
    total_study_2024 <- total_study_2014 * cpi_2014_to_2024
    
    # Calculate testing costs
    patients_per_arm <- n / 2
    testing_costs <- (snlb_cost * n) + (ly75_cost * patients_per_arm)
    
    # Total costs
    return(total_study_2024 + testing_costs)
  }
  
  # Function to simulate uncertainty around the point estimate
  simulate_uncertainty <- function(n, n_sims) {
    # Simulate deviations from point estimates (mean = 0)
    sim_intercept <- rnorm(n_sims, mean = 0, sd = se_intercept)
    sim_submission <- rnorm(n_sims, mean = 0, sd = se_submission_2014)
    sim_hta <- rnorm(n_sims, mean = 0, sd = se_hta)
    sim_duration <- rnorm(n_sims, mean = 0, sd = se_study_duration)
    sim_locations <- rnorm(n_sims, mean = 0, sd = se_study_locations)
    sim_sample_coef <- rnorm(n_sims, mean = 0, sd = se_sample_size_coef)
    sim_snlb <- rnorm(n_sims, mean = 0, sd = snlb_sd)
    sim_ly75 <- rnorm(n_sims, mean = 0, sd = ly75_sd)
    
    # Calculate total uncertainty
    study_uncertainty_2014 <- sim_intercept + sim_submission + sim_hta + 
      (sim_duration * duration_months) + sim_locations + (sim_sample_coef * n)
    
    study_uncertainty_2024 <- study_uncertainty_2014 * cpi_2014_to_2024
    
    patients_per_arm <- n / 2
    testing_uncertainty <- (sim_snlb * n) + (sim_ly75 * patients_per_arm)
    
    # Return uncertainty (deviations from point estimate)
    return(study_uncertainty_2024 + testing_uncertainty)
  }
  
  # Function to fit gamma distribution to uncertainty-adjusted costs
  get_gamma_ci <- function(point_estimate, uncertainty, probs = c(0.025, 0.975)) {
    # Add uncertainty to point estimate to get full cost distribution
    full_costs <- point_estimate + uncertainty
    
    # Remove any non-positive values (just in case)
    full_costs <- full_costs[full_costs > 0]
    
    # Method of moments to estimate gamma parameters
    mean_sim <- mean(full_costs)
    var_sim <- var(full_costs)
    
    # Gamma parameters
    shape <- (mean_sim^2) / var_sim
    rate <- mean_sim / var_sim
    
    # Calculate quantiles from gamma distribution
    lower <- qgamma(probs[1], shape = shape, rate = rate)
    upper <- qgamma(probs[2], shape = shape, rate = rate)
    
    return(list(lower = lower, upper = upper))
  }
  
  # Calculate results for each sample size
  results <- data.frame(
    sample_size = sample_sizes,
    mean_cost = NA_real_,
    lower_ci = NA_real_,
    upper_ci = NA_real_
  )
  
  for(i in seq_along(sample_sizes)) {
    n <- sample_sizes[i]
    
    # Get deterministic point estimate
    point_estimate <- calculate_point_estimate(n)
    
    # Simulate uncertainty
    uncertainty <- simulate_uncertainty(n, n_sims)
    
    # Get gamma-based confidence intervals
    gamma_ci <- get_gamma_ci(point_estimate, uncertainty)
    
    # Store results (using point estimate as the mean)
    results$mean_cost[i] <- point_estimate
    results$lower_ci[i] <- gamma_ci$lower
    results$upper_ci[i] <- gamma_ci$upper
  }
  
  return(results)
}