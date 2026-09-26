library(rstan)
library(dplyr)
library(tidyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# 1) LOAD DEPENDENCIES
# -----------------------------------------------------------------------------
prep <- readRDS("data/prepared_data.rds")
fit  <- readRDS("models/stan_fit.rds")

# Extract posteriors
post_draws <- rstan::extract(fit)
alpha_draw <- post_draws$alpha                 
beta_draw  <- post_draws$beta                  
tau_draw   <- post_draws$tau                   
L_Omega_draw <- post_draws$L_Omega             
w_draw     <- post_draws$w                     

n_draw <- length(alpha_draw)

# Helper: Sample MVN with Cholesky
rmvn_chol <- function(mu, L) {
  z <- rnorm(length(mu))
  as.vector(mu + L %*% z)
}

# -----------------------------------------------------------------------------
# 2) COUNTERFACTUAL SIMULATION
# -----------------------------------------------------------------------------
y_cf <- array(NA_real_, dim = c(n_draw, prep$T_post, prep$M))
y_last_pre <- as.vector(prep$y_pre[prep$T_pre, ])

cat("Simulating counterfactual posteriors...\n")
for (s in 1:n_draw) {
  L_Sigma <- diag(tau_draw[s, ]) %*% L_Omega_draw[s, , ]
  y_prev <- y_last_pre
  
  for (t in 1:prep$T_post) {
    mu_t <- alpha_draw[s] * y_prev + as.vector(prep$X_post[t, , ] %*% beta_draw[s, ]) + w_draw[s, ]
    y_new <- rmvn_chol(mu_t, L_Sigma)
    y_cf[s, t, ] <- y_new
    y_prev <- y_new
  }
}

# -----------------------------------------------------------------------------
# 3) DIFFERENCE-IN-DIFFERENCES LOGIC
# -----------------------------------------------------------------------------
compute_diff_in_diff <- function(y_cf, y_post, y_mu, y_sd, treated_idx, control_idx, scale = "percent") {
  
  M <- length(y_mu)
  T_post <- dim(y_cf)[2]
  n_draw <- dim(y_cf)[1]
  effect_array <- array(NA_real_, dim = c(n_draw, T_post, M))
  
  for (m in 1:M) {
    if (scale == "percent") {
      # Convert to original units
      obs_orig <- y_post[, m] * y_sd[m] + y_mu[m]     
      cf_orig_mat <- y_cf[, , m] * y_sd[m] + y_mu[m]  
      obs_mat <- matrix(obs_orig, nrow = n_draw, ncol = T_post, byrow = TRUE)
      
      zero_mask <- (abs(cf_orig_mat) < 1e-8)
      p_mat <- (obs_mat - cf_orig_mat) / cf_orig_mat
      p_mat[zero_mask] <- NA_real_
      effect_array[,,m] <- p_mat
    } else {
      # Standardized units
      effect_array[,,m] <- sweep(y_cf[,,m], 2, y_post[,m], function(cf, obs) obs - cf) 
    }
  }
  
  treated_mean <- matrix(NA_real_, nrow = n_draw, ncol = T_post)
  control_mean <- matrix(NA_real_, nrow = n_draw, ncol = T_post)
  
  for (d in 1:n_draw) {
    treated_mean[d, ] <- rowMeans(effect_array[d, , treated_idx, drop = FALSE], na.rm = TRUE)
    control_mean[d, ] <- rowMeans(effect_array[d, , control_idx, drop = FALSE], na.rm = TRUE)
  }
  
  diff_draw_time <- treated_mean - control_mean     
  diff_cum_draw  <- t(apply(diff_draw_time, 1, cumsum))  
  
  return(list(
    pointwise_diff = diff_draw_time,
    cumulative_diff = diff_cum_draw,
    effect_array = effect_array,
    control_mean = control_mean
  ))
}

treated_ids <- c(5531, 5551)
treated_idx <- match(treated_ids, prep$sensors)
control_idx <- setdiff(seq_along(prep$sensors), treated_idx)

# Compute Percent Scale
did_percent <- compute_diff_in_diff(y_cf, prep$y_post, prep$y_mu, prep$y_sd, treated_idx, control_idx, scale = "percent")

# Compute Standardized Scale (for plotting)
did_std <- compute_diff_in_diff(y_cf, prep$y_post, prep$y_mu, prep$y_sd, treated_idx, control_idx, scale = "standardized")

# Summaries
q_draw <- function(x, probs = c(0.025, 0.5, 0.975)) quantile(x, probs = probs, na.rm = TRUE)
diff_point_med <- apply(did_percent$pointwise_diff, 2, median, na.rm = TRUE)
diff_cum_q <- q_draw(did_percent$cumulative_diff[, prep$T_post])
p_cum_neg <- mean(did_percent$cumulative_diff[, prep$T_post] < 0, na.rm = TRUE)

cat("\n=========================================\n")
cat("Difference-in-Differences (Treated - Control) [%]\n")
cat("First post week median =", round(diff_point_med[1], 4), "\n")
cat("Last post week median  =", round(diff_point_med[prep$T_post], 4), "\n")
cat("Cumulative difference [%]:\n")
print(round(diff_cum_q, 4))
cat("P(Cumulative Diff < 0) = ", round(p_cum_neg, 3), "\n")
cat("=========================================\n\n")

# -----------------------------------------------------------------------------
# 4) GGPLOT VISUALIZATION (Obs - CF Effect)
# -----------------------------------------------------------------------------
get_ci <- function(mat) {
  data.frame(
    mean = colMeans(mat, na.rm = TRUE),
    lo   = apply(mat, 2, quantile, 0.025, na.rm = TRUE),
    hi   = apply(mat, 2, quantile, 0.975, na.rm = TRUE)
  )
}

df_ct <- get_ci(did_std$control_mean)
time_axis <- prep$weeks[prep$post_time_idx]
mask <- time_axis <= as.Date("2014-01-01") 

plot_list <- list()

for (st_id in treated_ids) {
  idx_st <- match(st_id, prep$sensors)
  df_st <- get_ci(did_std$effect_array[,,idx_st])
  
  df_plot <- data.frame(
    Date = time_axis[mask],
    Treated_Effect = df_st$mean[mask],
    Treated_Lower = df_st$lo[mask],
    Treated_Upper = df_st$hi[mask],
    Control_Effect = df_ct$mean[mask],
    Control_Lower = df_ct$lo[mask],
    Control_Upper = df_ct$hi[mask]
  )
  
  p <- ggplot(df_plot, aes(x = Date)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    
    # Controls
    geom_ribbon(aes(ymin = Control_Lower, ymax = Control_Upper, fill = "Control Avg"), alpha = 0.15) +
    geom_line(aes(y = Control_Effect, color = "Control Avg"), linetype = "dashed", linewidth = 1) +
    
    # Treated
    geom_ribbon(aes(ymin = Treated_Lower, ymax = Treated_Upper, fill = paste("Station", st_id)), alpha = 0.25) +
    geom_line(aes(y = Treated_Effect, color = paste("Station", st_id)), linewidth = 1) +
    
    scale_color_manual(name = "Legend", values = c("blue", "red"), breaks = c("Control Avg", paste("Station", st_id))) +
    scale_fill_manual(name = "Legend", values = c("blue", "red"), breaks = c("Control Avg", paste("Station", st_id))) +
    
    labs(title = paste("Treated Station", st_id, "vs Control Average"),
         y = "Effect (Obs - CF) [Std Units]", x = "Date") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  plot_list[[as.character(st_id)]] <- p
  print(p)
}