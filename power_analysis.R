############################################################
### Perform Power Analysis
### Written By: Sara Camilli, February 2021
############################################################

# Conduct power analyses on my datasets
  # 1. A priori - compute N, given alpha, power, and effect size
  # 2. Sensitivity - compute effect size, given alpha, power, and N
    # Set alpha to 0.05 and power to 0.8

library(pwr)


# Use the pwr.f2.test function for a general linear model
# u and v are the numerator and denominator degrees of freedom,
# f2 is the effect size measure, sig.level is alpha, and power is power

alpha <- 0.05
power <- 0.8

# Number of predictors, p
p <- 28 # model 1
p <- 20 # model 1 without first 2 age buckets & final buckets of other categories

# Number of observations, n
n1 <- 761 # breast cancer
  # n1 <- 32 # breast cancer, tumor-normal matched
n2 <- 7864 # pan-cancer
  # n2 <- X  # pan-cancer, tumor-normal matched
n3 <- 154 # cervical cancer, tumor-normal matched

# Degrees of freedom
u <- p - 1  
v1 <- n1 - p  
v2 <- n2 - p
v3 <- n3 - p

# Effect size, f2
range_of_possible_effect_sizes <- seq(from = 0.001, to = 0.1, by = 0.001)

############################################################
### 1. A priori Power Analysis
############################################################
#' Plots an a priori power analysis in which we are looking at the sample size (N)
#' we would need to achieve a range of possible effect sizes, given a particular
#' alpha, power, and number of model predictors
#' @param alpha the alpha (significance level) value
#' @param power the power value
#' @param effect_sizes a vector of effect sizes we might like to achieve
#' @param u the number of model predictors minus 1
get_N_vals <- function(alpha, power, effect_sizes, u) {
  N_vals <- lapply(range_of_possible_effect_sizes, function(es) {
    res <- pwr.f2.test(u = u, v = NULL, f2 = es, sig.level = alpha, power = power)
    return(res$v)
  })
  print(N_vals)      
  plot(x = effect_sizes, y = N_vals, main = "A priori Power Analysis", xlab = "Effect Size", ylab = "N")
  return(N_vals)
}

N_vals <- get_N_vals(alpha, power, range_of_possible_effect_sizes, u)

############################################################
### 2. Sensitivity power analysis
############################################################
#' Conducts a sensitivity power analysis
#' Prints the effect size we can expect to detect given a particular
#' alpha, power, sample size, and number of model predictors
#' @param alpha the alpha (significance level) value
#' @param power the power value
#' @param u the number of model predictors minus 1
#' @param v the degrees of freedom
get_effect_size <- function(alpha, power, u, v) {
  es_val_res <- pwr.f2.test(u = u, v = v, f2 = NULL, sig.level = alpha, power = power)
  print(es_val_res)      
  return(es_val_res$f2)
}

es_val1 <- get_effect_size(alpha, power, u, v1)
es_val2 <- get_effect_size(alpha, power, u, v2)
es_val3 <- get_effect_size(alpha, power, u, v3)

