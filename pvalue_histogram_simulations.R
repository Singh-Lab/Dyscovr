############################################################
### Simulate Weird P-Value Distributions
### Written By: Sara Geraghty, Sept. 2021
############################################################

library(broom)

# Draw values from common sample distributions
n <- 1000

norm_vals <- rnorm(n, 0, 1)
norm_vals2 <- rnorm(n, 7, 3)
beta_vals <- rbeta(n, 24.13, 7.67)
uniform_vals <- runif(n, min=1, max=3)
exp_vals <- rexp(n, rate = 1)
f_vals <- rf(n, df1=5, df2=2)


# Run sample linear models and plot the p-value distributions
# Linear models assume normally distributed x and y values

num_runs <- 150

p_vals <- c()
for (i in 1:num_runs) {
  norm_vals <- rnorm(n, 0, 1)
  norm_vals2 <- rnorm(n, 7, 3)
  #linearly_rel_vals <- unlist(lapply(norm_vals, function(x) x * sample(c(2,3,4), 1, replace = TRUE)))
  #beta_vals <- rbeta(n, 24.13, 7.67)
  beta_vals <- rbeta(n, 24.13, 2)
  #uniform_vals <- runif(n, min=1, max=3)
  #exp_vals <- rexp(n, rate = 1)
  #poly_rel_vals <- unlist(lapply(norm_vals, function(x) x^sample(c(2,3,4), 1, replace = TRUE)))
  #poly_rel_vals2 <- unlist(lapply(norm_vals2, function(x) x^sample(c(2,3,4), 1, replace = TRUE)))
  #collinear_vals <- unlist(lapply(norm_vals, function(x) x^sample(c(2,3,4), 1, replace = TRUE)))
  
  #lm_fit <- tidy(lm(norm_vals ~ norm_vals2))
  #lm_fit <- tidy(lm(linearly_rel_vals ~ norm_vals))
  #lm_fit <- tidy(lm(beta_vals ~ norm_vals))
  #lm_fit <- tidy(lm(norm_vals ~ beta_vals))
  #lm_fit <- tidy(lm(uniform_vals ~ norm_vals))
  #lm_fit <- tidy(lm(norm_vals ~ uniform_vals))
  #lm_fit <- tidy(lm(exp_vals ~ norm_vals))
  #lm_fit <- tidy(lm(norm_vals ~ exp_vals))
  #lm_fit <- tidy(lm(poly_rel_vals2 ~ norm_vals))
  #lm_fit <- tidy(lm(norm_vals2 ~ norm_vals + collinear_vals))
  #lm_fit <- tidy(lm(exp_vals ~ rexp(n, rate = 3)))
  lm_fit <- tidy(lm(beta_vals ~ rbeta(n, 17, 5)))
  
  p_vals <- c(p_vals, lm_fit$p.value)
  
}
hist(p_vals, main = "")  