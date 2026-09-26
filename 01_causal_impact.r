library(CausalImpact)
library(corrplot)
library(imputeTS)
library(ggplot2)
library(bsts)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# Evaluates the causal impact, plots the results, and prints the summaries
evaluate_impact <- function(data, pre_period, post_period) {
  impact <- CausalImpact(data, pre_period, post_period)
  
  plot(impact)
  print(summary(impact))
  
  plot(impact$model$bsts.model, "coefficients")
  print(summary(impact$model$bsts.model)$coefficients)
  
  return(impact)
}

# Generates a correlation plot for the specified variables
plot_correlation <- function(data, vars) {
  cor_matrix <- cor(data[, vars], use = "complete.obs")
  corrplot(cor_matrix, method = "number", type = "upper", tl.cex = 0.7)
}

# -----------------------------------------------------------------------------
# GLOBAL PARAMETERS
# -----------------------------------------------------------------------------
pre_period <- c(1, 108)
post_period <- c(109, 262)

# -----------------------------------------------------------------------------
# AUTOMATED MODEL 5504
# -----------------------------------------------------------------------------
evaluate_impact(data_5504[, -1], pre_period, post_period)

vars_of_interest_5504 <- c("soil_temperature_0_to_100cm_mean...C.", 
                           "temperature_2m_max...C.", 
                           "temperature_2m_mean...C.", 
                           "dew_point_2m_mean...C.",
                           "relative_humidity_2m_mean....")

plot_correlation(data_5504, vars_of_interest_5504)

# CLEANING: Removing 2m mean temp (high correlation)
data_5504_clean <- data_5504[, -1]
data_5504_clean <- data_5504_clean[, -2]

evaluate_impact(data_5504_clean, pre_period, post_period)


# -----------------------------------------------------------------------------
# AUTOMATED MODEL 5506
# -----------------------------------------------------------------------------
evaluate_impact(data_5506[, -1], pre_period, post_period)

vars_of_interest_5506 <- c("wind_speed_10m_max..km.h.", 
                           "temperature_2m_max...C.", 
                           "temperature_2m_mean...C.", 
                           "dew_point_2m_mean...C.",
                           "relative_humidity_2m_mean....",
                           "wind_speed_10m_mean..km.h.")

plot_correlation(data_5506, vars_of_interest_5506)

# CLEANING
data_5506_clean <- data_5506[, -1]
data_5506_clean <- data_5506_clean[, -2]
data_5506_clean <- data_5506_clean[, -15]

evaluate_impact(data_5506_clean, pre_period, post_period)


# -----------------------------------------------------------------------------
# AUTOMATED MODEL 5531
# -----------------------------------------------------------------------------
# na_kalman uses a univariate structural model (or ARIMA) to fill gaps based on series dynamics.
# "auto.arima" automatically finds seasonality and trend.
# data_5531$Media_Valore <- na_kalman(data_5531$Media_Valore, model = "auto.arima")

evaluate_impact(data_5531[, -1], pre_period, post_period)

vars_of_interest_5531 <- c("temperature_2m_max...C.", 
                           "relative_humidity_2m_mean....", 
                           "dew_point_2m_mean...C.", 
                           "surface_pressure_mean..hPa.",
                           "soil_temperature_0_to_100cm_mean...C.",
                           "temperature_2m_mean...C.",
                           "surface_pressure_max..hPa.")

plot_correlation(data_5531, vars_of_interest_5531)

# CLEANING
data_5531_clean <- data_5531[, -1]
data_5531_clean <- data_5531_clean[, -2]
data_5531_clean <- data_5531_clean[, -12]

evaluate_impact(data_5531_clean, pre_period, post_period)


# -----------------------------------------------------------------------------
# AUTOMATED MODEL 5542
# -----------------------------------------------------------------------------
# data_5542$Media_Valore <- na_kalman(data_5542$Media_Valore, model = "auto.arima")

evaluate_impact(data_5542[, -1], pre_period, post_period)

vars_of_interest_5542 <- c("soil_moisture_0_to_100cm_mean..m..m..", 
                           "temperature_2m_max...C.", 
                           "temperature_2m_mean...C.", 
                           "relative_humidity_2m_mean....",
                           "dew_point_2m_mean...C.",
                           "surface_pressure_mean..hPa.",
                           "surface_pressure_max..hPa.")

plot_correlation(data_5542, vars_of_interest_5542)

# CLEANING
data_5542_clean <- data_5542[, -1]
data_5542_clean <- data_5542_clean[, -2]
data_5542_clean <- data_5542_clean[, -12]

evaluate_impact(data_5542_clean, pre_period, post_period)


# -----------------------------------------------------------------------------
# AUTOMATED MODEL 5550
# -----------------------------------------------------------------------------
# data_5550$Media_Valore <- na_kalman(data_5550$Media_Valore, model = "auto.arima")

evaluate_impact(data_5550[, -1], pre_period, post_period)

vars_of_interest_5550 <- c("temperature_2m_max...C.", 
                           "temperature_2m_mean...C.", 
                           "soil_temperature_0_to_100cm_mean...C.", 
                           "dew_point_2m_mean...C.",
                           "relative_humidity_2m_mean....",
                           "soil_moisture_0_to_100cm_mean..m..m..",
                           "vapour_pressure_deficit_max..kPa.")

plot_correlation(data_5550, vars_of_interest_5550)

# CLEANING
data_5550_clean <- data_5550[, -1]
data_5550_clean <- data_5550_clean[, -2]
data_5550_clean <- data_5550_clean[, -18]

evaluate_impact(data_5550_clean, pre_period, post_period)


# -----------------------------------------------------------------------------
# AUTOMATED MODEL 5551
# -----------------------------------------------------------------------------
# data_5551$Media_Valore <- na_kalman(data_5551$Media_Valore, model = "auto.arima")

evaluate_impact(data_5551[, -1], pre_period, post_period)

vars_of_interest_5551 <- c("temperature_2m_max...C.", 
                           "temperature_2m_mean...C.", 
                           "soil_temperature_0_to_100cm_mean...C.", 
                           "dew_point_2m_mean...C.",
                           "relative_humidity_2m_mean....",
                           "soil_moisture_0_to_100cm_mean..m..m..",
                           "vapour_pressure_deficit_max..kPa.")

plot_correlation(data_5551, vars_of_interest_5551)

# CLEANING
data_5551_clean <- data_5551[, -1]
data_5551_clean <- data_5551_clean[, -2]

evaluate_impact(data_5551_clean, pre_period, post_period)


# -----------------------------------------------------------------------------
# AUTOMATED MODEL 5552
# -----------------------------------------------------------------------------
# data_5552$Media_Valore <- na_kalman(data_5552$Media_Valore, model = "auto.arima")

evaluate_impact(data_5552[, -1], pre_period, post_period)

vars_of_interest_5552 <- c("temperature_2m_max...C.", 
                           "temperature_2m_mean...C.", 
                           "soil_temperature_0_to_100cm_mean...C.", 
                           "dew_point_2m_mean...C.",
                           "relative_humidity_2m_mean....",
                           "soil_moisture_0_to_100cm_mean..m..m..",
                           "cloud_cover_mean....")

plot_correlation(data_5552, vars_of_interest_5552)

# CLEANING
data_5552_clean <- data_5552[, -1]
data_5552_clean <- data_5552_clean[, -2]
data_5552_clean <- data_5552_clean[, -18]

evaluate_impact(data_5552_clean, pre_period, post_period)


# -----------------------------------------------------------------------------
# GRAPHICAL OVERVIEW (FOREST PLOT)
# -----------------------------------------------------------------------------
results <- data.frame(
  Site = c("5504", "5550", "5551", "5506", "5531", "5552", "5542"),
  Effect = c(-22.0, -15.0, -12.0, -6.5, -4.6, 1.2, 2.7),
  LowerCI = c(-26.0, -25.0, -18.0, -12.0, -12.0, -7.3, -5.6),
  UpperCI = c(-16.0, -1.0, -4.7, -0.01, 3.4, 11.0, 12.0),
  Significant = c("Yes", "Yes", "Yes", "Yes", "No", "No", "No")
)

results$Site <- factor(results$Site, levels = results$Site[order(results$Effect)])

ggplot(results, aes(x = Effect, y = Site, color = Significant)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0.3, linewidth = 1) +
  geom_point(size = 4) +
  scale_color_manual(values = c("No" = "grey60", "Yes" = "#D55E00")) +
  labs(
    title = "Intervention Impact by Site (Forest Plot)",
    subtitle = "Relative change (%) with 95% confidence intervals",
    x = "Relative Change (%)",
    y = "Site Code",
    caption = "Note: Bars crossing 0 indicate non-statistically significant results."
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )

# -----------------------------------------------------------------------------
# CUSTOM MODEL (SITE 5504)
# -----------------------------------------------------------------------------
post_period_response <- data_5504$Media_Valore[post_period[1] : post_period[2]]
post_period_response <- na.omit(post_period_response)

data_5504$Media_Valore[post_period[1] : post_period[2]] <- NA

ss <- AddLocalLevel(list(), data_5504$Media_Valore)
bsts_model <- bsts(data_5504$Media_Valore ~ . - Settimana, ss, niter = 1000, data = data_5504)
summary(bsts_model)

impact_custom <- CausalImpact(bsts.model = bsts_model,
                              post.period.response = post_period_response)

plot(impact_custom)
summary(impact_custom)
summary(impact_custom, "report")