library(dplyr)
library(tidyr)
library(geosphere)

# -----------------------------------------------------------------------------
# 1) DATA LOADING & METADATA
# -----------------------------------------------------------------------------
col_names <- c("sensor_id", "week", "pollution", "cloud", 
               "dew_point", "humidity", "temp_max", "soil_moist", 
               "soil_temp", "pressure", "wind")

# Note: Update path to match your project structure (e.g., "data/dataset_finale_senza_NA.csv")
d <- read.csv("dataset_finale_senza_NA.csv", header = FALSE, sep = ",", skip = 1)
colnames(d) <- col_names

covariates <- c("cloud", "dew_point", "humidity", "temp_max", "soil_moist", "soil_temp", "pressure", "wind")
K <- length(covariates)

d2 <- d %>%
  mutate(week = as.Date(week)) %>%
  arrange(sensor_id, week)

sensors <- sort(unique(d2$sensor_id))
weeks   <- sort(unique(d2$week))

M <- length(sensors)
S <- length(weeks)

# Area C Intervention Date
cut_date <- as.Date("2012-01-16")

# Sensor coordinates metadata
# Best practice: move this to a "data/metadata_stazioni.csv" and read it
metadata <- data.frame(
  sensor_id = sensors,
  lat = c(45.49631883, 45.47606529, 45.53476819, 45.59538900, 45.52000944, 45.53976956, 45.39541003, 45.52342919, 
          45.50985564, 45.51782742, 45.52655611, 45.43611269, 45.54594739, 45.46334886, 45.55232831, 45.48363242, 
          45.50148956, 45.32452314, 45.54343664, 45.48027247, 45.44386047, 45.57171961, 45.54759700, 45.28196356, 
          45.46242074, 45.49958569, 45.47050100, 45.43220019, 45.54852147, 45.35428881, 45.39619767, 45.48455958),
  lon = c(9.19093444, 9.14178666, 9.23610903, 8.92203900, 9.51223153, 9.48689669, 8.91341828, 9.04460250, 
          9.08997419, 8.76656744, 8.73650211, 9.09741197, 8.77611000, 9.19532511, 9.22776572, 9.32736196, 
          8.80445883, 9.13452950, 9.08072867, 9.05686386, 9.16794489, 9.07824492, 9.16698300, 8.98857592, 
          8.88021388, 9.24732731, 9.19746075, 9.18215183, 8.84732669, 9.32924395, 9.28269897, 9.47134686)
)

# Distance matrix (Haversine in km)
coord_mat <- cbind(metadata$lon, metadata$lat)
D <- distm(coord_mat, coord_mat, fun = distHaversine) / 1000

# -----------------------------------------------------------------------------
# 2) FULL GRID & MATRICES CREATION
# -----------------------------------------------------------------------------
grid <- tidyr::expand_grid(week = weeks, sensor_id = sensors)
g <- grid %>% left_join(d2, by = c("week", "sensor_id")) %>% arrange(week, sensor_id)

pre_time_idx <- which(weeks < cut_date)
post_time_idx <- which(weeks >= cut_date)
T_pre <- length(pre_time_idx)
T_post <- length(post_time_idx)

if (T_pre < 2) stop("Not enough pre-intervention weeks to estimate AR(1).")

y_raw <- matrix(g$pollution, nrow = S, ncol = M, byrow = TRUE)
X_raw <- array(NA_real_, dim = c(S, M, K))
for (k in seq_along(covariates)) {
  X_raw[,,k] <- matrix(g[[covariates[k]]], nrow = S, ncol = M, byrow = TRUE)
}

# -----------------------------------------------------------------------------
# 3) STANDARDIZATION (PRE-PERIOD ONLY TO AVOID LEAKAGE)
# -----------------------------------------------------------------------------
y_pre_raw <- y_raw[pre_time_idx, , drop = FALSE]
y_mu <- colMeans(y_pre_raw)
y_sd <- apply(y_pre_raw, 2, sd)
y_sd[y_sd == 0] <- 1 

y_scaled <- sweep(y_raw, 2, y_mu, "-")
y_scaled <- sweep(y_scaled, 2, y_sd, "/")
y_pre <- y_scaled[pre_time_idx, , drop = FALSE]
y_post <- y_scaled[post_time_idx, , drop = FALSE]

X_scaled <- X_raw
for (k in 1:K) {
  x_pre_vec <- as.vector(X_raw[pre_time_idx, , k])
  x_mu <- mean(x_pre_vec)
  x_sd <- sd(x_pre_vec)
  if (is.na(x_sd) || x_sd == 0) x_sd <- 1
  X_scaled[,,k] <- (X_raw[,,k] - x_mu) / x_sd
}

X_pre <- X_scaled[pre_time_idx, , , drop = FALSE]
X_post <- X_scaled[post_time_idx, , , drop = FALSE]

# -----------------------------------------------------------------------------
# 4) SAVE PREPARED DATA
# -----------------------------------------------------------------------------
prepared_data <- list(
  y_pre = y_pre, y_post = y_post,
  X_pre = X_pre, X_post = X_post,
  y_mu = y_mu, y_sd = y_sd,
  weeks = weeks, sensors = sensors, cut_date = cut_date,
  pre_time_idx = pre_time_idx, post_time_idx = post_time_idx,
  T_pre = T_pre, T_post = T_post, M = M, K = K, D = D
)

# Note: Create a 'data/' directory if it doesn't exist
saveRDS(prepared_data, "data/prepared_data.rds")
cat("Data preparation complete. Saved to 'data/prepared_data.rds'\n")