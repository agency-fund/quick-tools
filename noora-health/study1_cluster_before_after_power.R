library(tidyverse)

#==========================================
# Cluster Before-After Design Power Calculator
# Study 1: Social Norms Intervention vs Default CCP for KMC
# Design: All clusters simultaneously experience control → washout → intervention
#==========================================

cluster_before_after_power <- function(
    p1,                     # baseline proportion (0–1)
    p2,                     # treatment proportion (0–1)
    clusters,               # number of hospitals (clusters)
    births_per_week,        # avg births per hospital per week
    weeks_control      = 2, # length of control period (weeks)
    weeks_washout      = 2, # length of washout period (weeks; data excluded)
    weeks_intervention = 2, # length of intervention period (weeks)
    icc_within         = 0.02,  # within-period intra-cluster correlation
    icc_between        = 0.01,  # between-period correlation (usually < icc_within)
    alpha              = 0.05,  # two-sided significance level
    use_t              = TRUE,  # use t-crit when clusters < 30
    time_trend         = 0      # linear time trend per period (optional)
) {
  
  # Input validation
  if (!is.numeric(p1) || p1 < 0 || p1 > 1) stop("p1 must be between 0 and 1")
  if (!is.numeric(p2) || p2 < 0 || p2 > 1) stop("p2 must be between 0 and 1")
  if (clusters < 2) stop("Need at least 2 clusters")
  if (!is.numeric(icc_within) || icc_within < 0 || icc_within > 1) stop("icc_within must be between 0 and 1")
  if (!is.numeric(icc_between) || icc_between < 0 || icc_between > 1) stop("icc_between must be between 0 and 1")
  if (icc_between > icc_within) warning("icc_between typically should be <= icc_within")
  
  # Cluster-period sizes (washout dropped)
  m0 <- weeks_control * births_per_week       # control period size
  m1 <- weeks_intervention * births_per_week  # intervention period size
  
  # Variance calculation for cluster before-after design
  # Based on Donner & Klar (2000) and modifications for temporal correlation
  
  # Outcome variance under null
  sigma2 <- p1 * (1 - p1)
  
  # Design effects for each period
  DE0 <- 1 + (m0 - 1) * icc_within  # control period
  DE1 <- 1 + (m1 - 1) * icc_within  # intervention period
  
  # Variance of difference between periods
  # Accounts for within-cluster correlation between time periods
  var_diff_per_cluster <- (sigma2 * DE0 / m0) + (sigma2 * DE1 / m1) - 
                         2 * sqrt((sigma2 * DE0 / m0) * (sigma2 * DE1 / m1)) * icc_between
  
  # Standard error of treatment effect
  se_beta <- sqrt(var_diff_per_cluster / clusters)
  
  # Adjust for potential time trend
  delta_adjusted <- abs(p2 - p1) - abs(time_trend)
  
  # Critical value
  if (use_t && clusters < 30) {
    df <- clusters - 1
    crit <- qt(1 - alpha/2, df = df)
  } else {
    crit <- qnorm(1 - alpha/2)
  }
  
  # Power calculation
  power <- pnorm(delta_adjusted / se_beta - crit) + 
           pnorm(-delta_adjusted / se_beta - crit)
  
  # Return comprehensive results
  list(
    design                = "cluster_before_after",
    p1                    = p1,
    p2                    = p2,
    clusters              = clusters,
    weeks_control         = weeks_control,
    weeks_washout         = weeks_washout,
    weeks_intervention    = weeks_intervention,
    births_per_week       = births_per_week,
    m0                    = m0,
    m1                    = m1,
    icc_within            = icc_within,
    icc_between           = icc_between,
    time_trend            = time_trend,
    design_effect_control = DE0,
    design_effect_treat   = DE1,
    sigma2                = sigma2,
    var_diff_per_cluster  = var_diff_per_cluster,
    se_estimate           = se_beta,
    critical_value        = crit,
    delta_unadjusted      = abs(p2 - p1),
    delta_adjusted        = delta_adjusted,
    power                 = power
  )
}

#==========================================
# Study 1 Power Analysis
#==========================================

# Average births per week across 8 Karnataka hospitals
births_per_week_karnataka <- mean(c(259.78, 169.56, 145.44, 133.33,
                                   125.00, 114.11, 102.99, 101.11)) / 31 * 7

cat("Study 1: Social Norms Intervention for KMC Practice\n")
cat("Design: Cluster before-after with washout\n")
cat("Baseline births per week per hospital:", round(births_per_week_karnataka, 1), "\n\n")

#===========================
# Primary power calculation
#===========================

primary_result <- cluster_before_after_power(
  p1                 = 0.081,  # baseline KMC practice rate
  p2                 = 0.181,  # expected rate with intervention (+10pp)
  clusters           = 8,      # 8 Karnataka hospitals
  births_per_week    = births_per_week_karnataka,
  weeks_control      = 2,
  weeks_washout      = 2,
  weeks_intervention = 2,
  icc_within         = 0.02,   # within-period ICC
  icc_between        = 0.01,   # between-period correlation (conservative)
  alpha              = 0.05,
  use_t              = TRUE
)

cat("Primary Analysis Results:\n")
cat("Expected power:", round(primary_result$power, 3), "\n")
cat("Standard error:", round(primary_result$se_estimate, 4), "\n")
cat("Effect size:", round(primary_result$delta_unadjusted, 3), "\n\n")

#===========================
# Sensitivity analyses
#===========================

# 1. Different ICC values
icc_scenarios <- data.frame(
  icc_within = c(0.01, 0.02, 0.03, 0.05),
  icc_between = c(0.005, 0.01, 0.015, 0.025)  # typically half of within-period ICC
)

icc_results <- map2_dfr(icc_scenarios$icc_within, icc_scenarios$icc_between, function(icc_w, icc_b) {
  res <- cluster_before_after_power(
    p1 = 0.081, p2 = 0.181, clusters = 8,
    births_per_week = births_per_week_karnataka,
    weeks_control = 2, weeks_washout = 2, weeks_intervention = 2,
    icc_within = icc_w, icc_between = icc_b,
    alpha = 0.05, use_t = TRUE
  )
  data.frame(
    icc_within = icc_w,
    icc_between = icc_b,
    power = res$power,
    se = res$se_estimate
  )
})

cat("Sensitivity Analysis - ICC Values:\n")
print(icc_results, digits = 3)
cat("\n")

# 2. Different effect sizes and cluster numbers
effect_sizes <- seq(0.05, 0.20, by = 0.025)  # 5pp to 20pp
cluster_numbers <- 6:10

scenario_grid <- expand.grid(
  delta = effect_sizes,
  clusters = cluster_numbers
)

scenario_results <- map2_dfr(scenario_grid$delta, scenario_grid$clusters, function(d, k) {
  res <- cluster_before_after_power(
    p1 = 0.081, p2 = 0.081 + d, clusters = k,
    births_per_week = births_per_week_karnataka,
    weeks_control = 2, weeks_washout = 2, weeks_intervention = 2,
    icc_within = 0.02, icc_between = 0.01,
    alpha = 0.05, use_t = TRUE
  )
  data.frame(
    delta = d,
    clusters = k,
    power = res$power
  )
})

# 3. Time trend sensitivity
time_trend_scenarios <- c(0, 0.01, 0.02, -0.01, -0.02)  # ±1-2pp per period

trend_results <- map_dfr(time_trend_scenarios, function(trend) {
  res <- cluster_before_after_power(
    p1 = 0.081, p2 = 0.181, clusters = 8,
    births_per_week = births_per_week_karnataka,
    weeks_control = 2, weeks_washout = 2, weeks_intervention = 2,
    icc_within = 0.02, icc_between = 0.01,
    alpha = 0.05, use_t = TRUE, time_trend = trend
  )
  data.frame(
    time_trend = trend,
    power = res$power,
    delta_adjusted = res$delta_adjusted
  )
})

cat("Sensitivity Analysis - Time Trends:\n")
print(trend_results, digits = 3)
cat("\n")

#===========================
# Visualization
#===========================

# Power by effect size and clusters
power_plot <- ggplot(scenario_results, aes(x = clusters, y = power, color = factor(delta))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.8, linetype = "dashed", alpha = 0.7) +
  scale_color_viridis_d(
    name = "Effect size\n(Δ pp)",
    labels = paste0("+", scenario_results %>% distinct(delta) %>% pull(delta) * 100, "pp")
  ) +
  scale_x_continuous(breaks = cluster_numbers) +
  labs(
    x = "Number of Hospitals (Clusters)",
    y = "Statistical Power",
    title = "Study 1: Power for Cluster Before-After Design",
    subtitle = "Social norms intervention for KMC practice\n(Baseline: 8.1%, ICC_within: 0.02, ICC_between: 0.01)"
  ) +
  ylim(0, 1) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

print(power_plot)

#===========================
# Design recommendations
#===========================

cat("DESIGN RECOMMENDATIONS:\n")
cat("1. Washout period: 2 weeks may be insufficient for behavioral interventions\n")
cat("   Consider 4+ weeks to prevent contamination\n\n")
cat("2. Temporal trends: Include hospital-specific time trends in analysis\n")
cat("   Consider external controls if available\n\n")
cat("3. ICC assumptions: Validate with pilot data if possible\n")
cat("   Between-period correlation is critical but often unknown\n\n")
cat("4. Sample size: With 8 hospitals and 10pp effect, power =", round(primary_result$power, 2), "\n")
cat("   Consider 10+ hospitals for more robust results\n\n")
cat("5. Analysis plan: Use mixed-effects models with:\n")
cat("   - Hospital random effects\n")
cat("   - Period fixed effects\n")
cat("   - Robust standard errors\n")