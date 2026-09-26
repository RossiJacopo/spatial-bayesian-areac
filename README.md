# Causal Impact of Milan's Area C on NO2 Levels

This repository contains the statistical analysis and codebase used to evaluate the impact of **Area C** (a congestion charge and traffic restriction zone introduced in Milan on January 16, 2012) on atmospheric pollution, specifically Nitrogen Dioxide ($NO_2$) levels.

The project employs two complementary causal inference methodologies to isolate the intervention's effect from meteorological confounders and spatial dynamics.

## Methodological Approach

1. **Univariate Causal Impact (BSTS):** 
   A site-specific analysis using Bayesian Structural Time Series (`CausalImpact` package) to estimate the intervention's effect on individual monitoring stations independently.
2. **Spatial Bayesian Counterfactual Modeling (Stan):** 
   A custom hierarchical Bayesian model written in Stan. It explicitly models temporal autocorrelation (AR1) and spatial dependency (via an exponential covariance function based on Haversine distances) to generate multi-step ahead counterfactual scenarios for the treated stations versus a control group (Difference-in-Differences).
