# Required packages
library(glmnet)
library(ncvreg)
library(ggplot2)

# Function to generate simulated data
generate_simulated_data <- function(n, p) {
  X <- matrix(runif(n * p) ^ rep(1:p, each = n) / 100, n, p)
  
  omega <- 7
  th <- omega * sqrt(log(p) / n)
  beta <- ceiling(runif(p, min = th, max = th + 2)) * (1:p <= 5) -
    ceiling(runif(p, min = th - 1, max = th - 0.1)) * (10 < 1:p & 1:p <= 15)
  
  t <- 2.5 * (1:n - 0.5) / n
  f <- -t * sin(-t^2)
  
  e <- rnorm(n, sd = 0.5)
  y <- X %*% beta + f + e
  
  list(X = X, beta = beta, t = t, f = f, y = y, n = n, p = p, e = e)
}

# Adaptive Lasso function
adaptive_lasso <- function(x, y) {
  beta.init <- solve(t(x) %*% x, tol = 1e-100) %*% t(x) %*% y
  
  w <- abs(beta.init)
  x2 <- scale(x, center = FALSE, scale = 1 / w)
  
  la_eq <- cv.glmnet(x2, y, nfolds = 10, intercept = FALSE, lambda = seq(0.001, 0.05, length.out = 20))
  coefs <- coef(la_eq, s = "lambda.min")
  S1 <- coefs@i[-1]  # remove the intercept index
  S1_beta <- coefs@x[-1]
  
  list(X_S1 = x2[, S1, drop = FALSE], X_S1c = x[, -S1, drop = FALSE], S1_beta = S1_beta, S1 = S1 + 1, coefs = coefs)
}

# SCAD penalty function
scad_penalty <- function(x, y) {
  scadm <- ncvreg(x, y, family = "gaussian", penalty = "SCAD", lambda = seq(0.05, 0.001, length.out = 100))
  S1_beta <- scadm$beta[-1, 80]
  S1 <- which(S1_beta != 0)
  
  list(X_S1 = x[, S1, drop = FALSE], X_S1c = x[, -S1, drop = FALSE], S1_beta = S1_beta[S1_beta != 0], S1 = S1, coefs = S1_beta)
}

# Function to optimize threshold
optimize_threshold <- function(betaWR) {
  thr_seq <- seq(min(abs(betaWR)) * 3, max(abs(betaWR)) - 0.1, length.out = 50)
  AMPE <- sapply(thr_seq, function(thr) mean((betaWR - thr)^2))
  return(thr_seq[which.min(AMPE)])
}

# Function to compute PSE
compute_PSE <- function(obj, x, y, beta_S2, S2index, lambdar) {
  M <- diag(length(y)) - obj$X_S1 %*% solve(t(obj$X_S1) %*% obj$X_S1, tol = 1e-100) %*% t(obj$X_S1)
  Tr <- sum((t(beta_S2) %*% (t(x[, S2index, drop = FALSE]) %*% M %*% x[, S2index, drop = FALSE]) %*% beta_S2) / var(y - obj$X_S1 %*% obj$S1_beta))
  return(obj$S1_beta - (min(((length(S2index) - 2) / Tr), 1)) * (abs(obj$S1_beta - obj$S1_beta)))
}

# Weighted Ridge function
weighted_ridge <- function(obj, x, y, lambdar) {
  X_S1 <- obj$X_S1
  X_S1c <- obj$X_S1c
  betaWR <- solve(t(X_S1c) %*% X_S1c + lambdar * diag(ncol(X_S1c)), tol = 1e-100) %*% t(X_S1c) %*% y
  
  opt_thr <- optimize_threshold(betaWR)
  beta_S2 <- betaWR[abs(betaWR) > opt_thr]
  S2index <- which(abs(betaWR) > opt_thr)
  
  betaPSE <- compute_PSE(obj, x, y, beta_S2, S2index, lambdar)
  
  list(betaPSE = betaPSE, opt_thr = opt_thr, X_S1 = X_S1, X_S1c = X_S1c, betaWR = betaWR)
}

# Performance metrics
rmsebeta <- function(realbeta, estbeta) {
  sqrt(mean((realbeta - estbeta) ^ 2))
}

msefhat <- function(real_f, est_f) {
  mean((real_f - est_f) ^ 2)
}

# Main simulation function with performance measures
run_simulation <- function(sim, n, p) {
  results <- vector("list", sim)
  rmse_results <- data.frame(sim = 1:sim, alas = NA, scad = NA, PSE_alas = NA, PSE_scad = NA)
  mse_results <- data.frame(sim = 1:sim, alas = NA, scad = NA, PSE_alas = NA, PSE_scad = NA)
  
  # Initialize matrices to store fitted values across simulations
  fitted_alas_all <- matrix(0, n, sim)
  fitted_scad_all <- matrix(0, n, sim)
  fitted_PSE_alas_all <- matrix(0, n, sim)
  fitted_PSE_scad_all <- matrix(0, n, sim)
  
  for (s in 1:sim) {
    data <- generate_simulated_data(n, p)
    alas <- adaptive_lasso(data$X, data$y)
    scad <- scad_penalty(data$X, data$y)
    PSE_alas <- weighted_ridge(alas, data$X, data$y, mean(seq(0.025, 1, length.out = 30)))
    PSE_scad <- weighted_ridge(scad, data$X, data$y, mean(seq(0.025, 1, length.out = 30)))
    
    # Calculate RMSE for beta
    rmse_results[s, "alas"] <- rmsebeta(data$beta, alas$S1_beta)
    rmse_results[s, "scad"] <- rmsebeta(data$beta, scad$S1_beta)
    rmse_results[s, "PSE_alas"] <- rmsebeta(data$beta, PSE_alas$betaPSE)
    rmse_results[s, "PSE_scad"] <- rmsebeta(data$beta, PSE_scad$betaPSE)
    
    # Calculate fitted values for nonparametric component
    fitted_alas <- PSE_alas$X_S1 %*% PSE_alas$betaPSE
    fitted_scad <- PSE_scad$X_S1 %*% PSE_scad$betaPSE
    
    # Store fitted values in the matrix
    fitted_alas_all[, s] <- fitted_alas
    fitted_scad_all[, s] <- fitted_scad
    fitted_PSE_alas_all[, s] <- fitted_alas
    fitted_PSE_scad_all[, s] <- fitted_scad
    
    # Calculate MSE for nonparametric component
    mse_results[s, "alas"] <- msefhat(data$f, fitted_alas)
    mse_results[s, "scad"] <- msefhat(data$f, fitted_scad)
    mse_results[s, "PSE_alas"] <- msefhat(data$f, fitted_alas)
    mse_results[s, "PSE_scad"] <- msefhat(data$f, fitted_scad)
    
    # Store results
    results[[s]] <- list(alas = alas, scad = scad, PSE_alas = PSE_alas, PSE_scad = PSE_scad)
    
    # Output iteration info
    message(s, "th Simulation ends")
  }
  
  library(ggplot2)
  library(gridExtra)
  
  # Calculate the average fitted values across simulations
  avg_fitted_alas <- rowMeans(fitted_alas_all, na.rm = TRUE)
  avg_fitted_scad <- rowMeans(fitted_scad_all, na.rm = TRUE)
  avg_fitted_PSE_alas <- rowMeans(fitted_PSE_alas_all, na.rm = TRUE)
  avg_fitted_PSE_scad <- rowMeans(fitted_PSE_scad_all, na.rm = TRUE)
  
  # Prepare data for RMSE bar plot
  rmse_avg <- colMeans(rmse_results[, -1], na.rm = TRUE)
  rmse_df <- data.frame(Method = names(rmse_avg), RMSE = abs(scale(rmse_avg)))
  
  # Prepare data for MSE bar plot
  mse_avg <- colMeans(mse_results[, -1], na.rm = TRUE)
  mse_df <- data.frame(Method = names(mse_avg), MSE = mse_avg)
  
  # Create RMSE bar plot
  rmse_plot <- ggplot(rmse_df, aes(x = Method, y = RMSE, fill = Method)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(title = "Average RMSE of Beta Estimates", y = "Scaled RMSE", x = "") +
    scale_y_continuous(labels = scales::number_format(scale = 1, suffix = "")) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14))
  
  # Create MSE bar plot
  mse_plot <- ggplot(mse_df, aes(x = Method, y = MSE, fill = Method)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(title = "Average MSE of Nonparametric Component", y = "MSE", x = "") +
    scale_y_continuous(labels = scales::number_format(scale = 1, suffix = "")) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14))
  
  # Create averaged fitted curve plot
  fitted_plot <- ggplot() +
    geom_line(aes(x = data$t, y = data$f, color = "Real"), size = 1.2) +
    geom_line(aes(x = data$t, y = avg_fitted_alas, color = "ALASSO"), linetype = "dashed", size = 1) +
    geom_line(aes(x = data$t, y = avg_fitted_scad, color = "SCAD"), linetype = "dashed", size = 1) +
    geom_line(aes(x = data$t, y = avg_fitted_PSE_alas, color = "PSE_ALASSO"), linetype = "dashed", size = 1) +
    geom_line(aes(x = data$t, y = avg_fitted_PSE_scad, color = "PSE_SCAD"), linetype = "dashed", size = 1) +
    labs(title = "Average Fitted Curves vs Real Curve", x = "t", y = "f(t)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14)) +
    scale_color_manual(values = c("Real" = "black", "ALASSO" = "red", "SCAD" = "blue", 
                                  "PSE_ALASSO" = "green", "PSE_SCAD" = "purple")) +
    guides(color = guide_legend(title = "Method"))
  
  # Arrange plots in a grid
  grid.arrange(rmse_plot, mse_plot, fitted_plot, ncol = 2, nrow = 2)
  
  return(list(results = results, rmse_results = rmse_results, mse_results = mse_results))
}

# Example of running the simulation
sim_results <- run_simulation(50, 50, 100) #sim, n, p)
rmse_results=sim_results$rmse_results
mse_results=sim_results$mse_results
performance_df <- data.frame(
  Simulation = 1:nrow(rmse_results),
  RMSE_alas = abs(scale(rmse_results$alas)),
  RMSE_scad = abs(scale(rmse_results$scad)),
  RMSE_PSE_alas = abs(scale(rmse_results$PSE_alas)),
  RMSE_PSE_scad = abs(scale(rmse_results$PSE_scad)),
  MSE_alas = abs(scale(mse_results$alas)),
  MSE_scad = abs(scale(mse_results$scad)),
  MSE_PSE_alas = abs(scale(mse_results$PSE_alas)),
  MSE_PSE_scad = abs(scale(mse_results$PSE_scad))
)

# Melt the data frame for ggplot
library(reshape2)
melted_performance_df <- melt(performance_df, id.vars = "Simulation")

# Plotting the relationship between parametric (RMSE) and nonparametric (MSE) performance
ggplot(melted_performance_df, aes(x = Simulation)) +
  geom_line(aes(y = value, color = variable), size = 1) +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  labs(title = "Performance Relationship Between Parametric and Nonparametric Components",
       x = "Simulation",
       y = "Performance Measure") +
  theme_minimal()


library(gridExtra)

# Prepare a summary table with average MSE and RMSE values
library(gridExtra)

# Prepare a summary table with average MSE and RMSE values
summary_table <- data.frame(
  Method = c("ALASSO", "SCAD", "PSE_ALASSO", "PSE_SCAD"),
  Avg_RMSE = c(mean(abs(scale(rmse_results$alas)), na.rm = TRUE), 
               mean(abs(scale(rmse_results$scad)), na.rm = TRUE), 
               mean(abs(scale(rmse_results$PSE_alas)), na.rm = TRUE), 
               mean(abs(scale(rmse_results$PSE_scad)), na.rm = TRUE)),
  Avg_MSE = c(mean(mse_results$alas, na.rm = TRUE), 
              mean(mse_results$scad, na.rm = TRUE), 
              mean(mse_results$PSE_alas, na.rm = TRUE), 
              mean(mse_results$PSE_scad, na.rm = TRUE))
)

# Convert the summary table to a grid object with grid lines for plotting
table_plot <- tableGrob(summary_table, rows = NULL, 
                        theme = ttheme_minimal(core = list(fg_params = list(hjust = 0.5, x = 0.5),
                                                           bg_params = list(fill = "white", col = "black", lwd = 0.5))))

# Display the table
grid.arrange(table_plot)
