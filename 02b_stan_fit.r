library(rstan)

# -----------------------------------------------------------------------------
# 1) LOAD PREPARED DATA
# -----------------------------------------------------------------------------
prep <- readRDS("data/prepared_data.rds")

y_list <- lapply(1:prep$T_pre, function(t) as.vector(prep$y_pre[t, ]))
X_list <- lapply(1:prep$T_pre, function(t) prep$X_pre[t, , ])

stan_data <- list(
  T = prep$T_pre,
  M = prep$M,
  K = prep$K,
  y = y_list,
  X = X_list,
  dist_mat = prep$D,
  l = max(prep$D) / 3,
  Sigma0 = cov(prep$y_pre)
)

init_fun <- function() list(
  alpha = 0,
  beta = rep(0, prep$K),
  sigma_w = 1,
  z_w = rep(0, prep$M),
  tau = rep(1, prep$M),
  L_Omega = diag(prep$M)
)

# -----------------------------------------------------------------------------
# 2) FIT OR LOAD STAN MODEL
# -----------------------------------------------------------------------------
model_path <- "models/stan_fit.rds"
stan_file <- "models/prediction_model.stan"

if (file.exists(model_path)) {
  
  cat("Loading existing Stan model from disk...\n")
  fit <- readRDS(model_path)
  
} else {
  
  cat("Compiling and sampling Stan model. This may take a while...\n")
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())
  
  model <- stan_model(file = stan_file)
  
  fit <- sampling(
    model,
    data = stan_data,
    chains = 4,
    init = init_fun,
    iter = 2000,
    warmup = 500,
    seed = 123,
    refresh = 10,
    control = list(adapt_delta = 0.995, max_treedepth = 15)
  )
  
  # Note: Create a 'models/' directory if it doesn't exist
  saveRDS(fit, model_path)
  cat("Model training complete and saved to disk.\n")
}

print(fit, pars = c("alpha", "beta", "sigma_w", "tau"))