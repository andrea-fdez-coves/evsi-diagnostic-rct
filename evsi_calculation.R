#========================================================
# EVSI for CAU and Intervention using voi::evsi (IS)
#========================================================

rm(list = ls())
if(dev.cur() != 1) dev.off()
set.seed(12345) 

# Load required packages
library(mgcv)
library(ggplot2)
library(dplyr)
library(progress)
library(voi)

#--------------------------------------------------------
# 1. Load and prepare data
#--------------------------------------------------------
theta        <- as.matrix(read.csv("params_all_bc.csv"))        
cost_data    <- as.matrix(read.csv("costs_bc.csv"))
effects_data <- as.matrix(read.csv("effects_bc.csv"))
wtp <- 20000

# Net Monetary Benefit matrix (each column a decision option)
NB <- (effects_data * wtp) - cost_data
colnames(NB) <- c("S1", "S2", "S6", "CAU")

# PSA inputs (parameters that will be learned)
decision_model_psa <- as.data.frame(theta[, c(
  "CAU_sensitivity", "CAU_specificity",
  "Intervention_sensitivity", "Intervention_specificity",
  "prev"
)])

#--------------------------------------------------------
# 2. Helper functions
#--------------------------------------------------------
calculate_test_probabilities <- function(prev, se, sp) {
  p_positive <- prev * se + (1 - prev) * (1 - sp) #probability of a test giving a positive result
  p_negative <- prev * (1 - se) + (1 - prev) * sp #probability of a test giving a negative result
  return(list(p_positive = p_positive, p_negative = p_negative))
}

calculate_ppv <- function(prev, se, sp) {
  (prev * se) / (prev * se + (1 - prev) * (1 - sp))
}

calculate_npv <- function(prev, se, sp) {
  ((1 - prev) * sp) / ((1 - prev) * sp + prev * (1 - se))
}

#--------------------------------------------------------
# 3. Data generation function
#--------------------------------------------------------
datagen_fn_IS_diagnostic <- function(inputs, n = 100) {
  n_arm1 <- floor(n / 2)
  n_arm2 <- ceiling(n / 2)   # ensures total = n even if n odd
  n_rows <- nrow(inputs)
  
  result <- data.frame(
    CAU_TP = integer(n_rows),
    CAU_FP = integer(n_rows),
    CAU_FN = integer(n_rows),
    CAU_TN = integer(n_rows),
    INT_TP = integer(n_rows),
    INT_FP = integer(n_rows),
    INT_FN = integer(n_rows),
    INT_TN = integer(n_rows)
  )
  
  for (i in 1:n_rows) {
    prev <- inputs$prev[i]
    CAU_se <- inputs$CAU_sensitivity[i]
    CAU_sp <- inputs$CAU_specificity[i]
    INT_se <- inputs$Intervention_sensitivity[i]
    INT_sp <- inputs$Intervention_specificity[i]
    
    # Simulate one arm given se, sp, prev, and arm size
    sim_arm <- function(se, sp, prev, n_arm) {
      probs <- c(
        prev * se,                  # TP
        (1 - prev) * (1 - sp),      # FP
        prev * (1 - se),             # FN
        (1 - prev) * sp              # TN
      )
      as.vector(rmultinom(1, size = n_arm, prob = probs))
    }
    
    CAU_counts <- sim_arm(CAU_se, CAU_sp, prev, n_arm1)
    INT_counts <- sim_arm(INT_se, INT_sp, prev, n_arm2)
    
    result$CAU_TP[i] <- CAU_counts[1]
    result$CAU_FP[i] <- CAU_counts[2]
    result$CAU_FN[i] <- CAU_counts[3]
    result$CAU_TN[i] <- CAU_counts[4]
    
    result$INT_TP[i] <- INT_counts[1]
    result$INT_FP[i] <- INT_counts[2]
    result$INT_FN[i] <- INT_counts[3]
    result$INT_TN[i] <- INT_counts[4]
  }
  
  return(result)
}

#--------------------------------------------------------
# 4. Likelihood function
#--------------------------------------------------------
likelihood_fn_IS_diagnostic <- function(Y, inputs, n = 100) {
  n_arm1 <- floor(n / 2)
  n_arm2 <- ceiling(n / 2)
  n_psa <- nrow(inputs)
  
  # Extract observed counts
  CAU_obs <- as.numeric(Y[1, c("CAU_TP", "CAU_FP", "CAU_FN", "CAU_TN")])
  INT_obs <- as.numeric(Y[1, c("INT_TP", "INT_FP", "INT_FN", "INT_TN")])
  
  log_lik <- numeric(n_psa)
  
  for (s in 1:n_psa) {
    prev <- inputs$prev[s]
    CAU_se <- inputs$CAU_sensitivity[s]
    CAU_sp <- inputs$CAU_specificity[s]
    INT_se <- inputs$Intervention_sensitivity[s]
    INT_sp <- inputs$Intervention_specificity[s]
    
    probs_CAU <- c(prev * CAU_se,
                   (1 - prev) * (1 - CAU_sp),
                   prev * (1 - CAU_se),
                   (1 - prev) * CAU_sp)
    
    probs_INT <- c(prev * INT_se,
                   (1 - prev) * (1 - INT_sp),
                   prev * (1 - INT_se),
                   (1 - prev) * INT_sp)
    
    log_lik_CAU <- dmultinom(CAU_obs, size = n_arm1, prob = probs_CAU, log = TRUE)
    log_lik_INT <- dmultinom(INT_obs, size = n_arm2, prob = probs_INT, log = TRUE)
    
    log_lik[s] <- log_lik_CAU + log_lik_INT
  }
  
  # Scale to avoid numerical underflow
  max_log <- max(log_lik)
  exp(log_lik - max_log)
}

#--------------------------------------------------------
# 5. EVSI calculation
#--------------------------------------------------------
sample_sizes <- c(30, 100, 300,  600, 900, 1200)

gam_formula <- ~ te(CAU_sensitivity, CAU_specificity,
                    Intervention_sensitivity, Intervention_specificity,
                    prev, k = 3)

evsi_results <- data.frame(sample_size = sample_sizes, evsi = NA)

for (i in seq_along(sample_sizes)) {
  size <- sample_sizes[i]
  cat("Calculating EVSI for sample size:", size, "\n")
  
  evsi_result <- evsi(
    outputs      = NB,                             
    inputs       = decision_model_psa,             
    datagen_fn   = datagen_fn_IS_diagnostic,      
    pars         = c("CAU_sensitivity", "CAU_specificity",
                     "Intervention_sensitivity", "Intervention_specificity",
                     "prev"),                       
    likelihood   = likelihood_fn_IS_diagnostic,      
    n            = size,                             
    method       = "is",                             
    Q            = 5000,                               
    npreg_method = "gam",                              
    gam_formula  = gam_formula                         
  )
  
  evsi_results$evsi[i] <- evsi_result$evsi
}

print(evsi_results)


#--------------------------------------------------------
# 6. Trial cost calculation
#--------------------------------------------------------

# Calculate trial costs
trial_costs <- calculate_trial_costs(sample_sizes)

print(trial_costs)

#--------------------------------------------------------
# 7. ENBS calculation
#--------------------------------------------------------
population_disease <- 27458

# Merge EVSI results with trial costs and calculate ENBS
enbs_results <- evsi_results %>%
  left_join(trial_costs, by = "sample_size") %>%
  mutate(enbs = (evsi * population_disease) - mean_cost) %>%
  select(sample_size, enbs)

print(enbs_results)
