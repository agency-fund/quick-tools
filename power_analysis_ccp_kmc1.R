#==========================================
# Function to compute clustered CRT power
#==========================================
cluster_power <- function(
    p1,                    # Control arm proportion practicing KMC (e.g. 0.081)
    p2,                    # Intervention arm proportion practicing KMC (e.g. 0.20)
    clusters,              # Number of hospital clusters
    weeks_per_period,      # Weeks in each period
    births_per_week,       # Avg births per hospital per week
    icc,                   # Intra-cluster correlation
    alpha = 0.05,          # Two-sided significance level
    sweep_weeks = FALSE,   # If TRUE, returns power for 1:max_weeks
    max_weeks = 6          # Max weeks to sweep over when sweep_weeks=TRUE
) {
  # 1) Compute cluster size and design effect
  m   <- weeks_per_period * births_per_week
  DE  <- 1 + (m - 1) * icc
  
  # 2) Effective N per arm
  n_eff <- (clusters * m) / DE
  
  # 3) Compute effect size & variance
  delta    <- abs(p2 - p1)
  var_diff <- p1 * (1 - p1) / n_eff + p2 * (1 - p2) / n_eff
  se_diff  <- sqrt(var_diff)
  
  # 4) Choose critical value: t‐critical if few clusters, else z‐critical
  if (clusters < 10) {
    df   <- 2 * (clusters - 1)         # degrees of freedom
    crit <- qt(1 - alpha/2, df = df)   # two-sided t critical value
  } else {
    crit <- qnorm(1 - alpha/2)         # normal z critical value
  }
  
  # 5) Compute power
  power <- pnorm(delta / se_diff - crit)
  
  # 6) Prepare result list
  result <- list(
    weeks_per_period    = weeks_per_period,
    births_per_week     = births_per_week,
    cluster_size        = m,
    design_effect       = DE,
    raw_N_per_arm       = clusters * m,
    effective_N_per_arm = n_eff,
    estimated_power     = power
  )
  
  # 7) Optionally sweep period lengths
  if (sweep_weeks) {
    weeks_vec <- seq(1, max_weeks)
    power_vec <- sapply(weeks_vec, function(w) {
      m2   <- w * births_per_week
      DE2  <- 1 + (m2 - 1) * icc
      n2   <- (clusters * m2) / DE2
      se2  <- sqrt(p1 * (1 - p1) / n2 + p2 * (1 - p2) / n2)
      pnorm(delta / se2 - crit)
    })
    result$sweep <- data.frame(
      weeks_per_period = weeks_vec,
      power_estimate   = round(power_vec, 3)
    )
  }
  
  return(result)
}

#===========================
# Example usage:
#===========================
res <- cluster_power(
  p1               = 0.31,    # baseline practice rate
  p2               = 0.56,    # target practice rate
  clusters         = 3,       # number of hospitals
  weeks_per_period = 2,       # weeks per arm
  births_per_week  = 60,      # avg births per week
  icc              = 0.02,    # intra-cluster correlation
  alpha            = 0.05,    # significance level
  sweep_weeks      = TRUE,    # include sweep of 1:max_weeks
  max_weeks        = 6        # max weeks to sweep
)

print(res$estimated_power)  # single power estimate
print(res$sweep)            # power by 1:6 week scenarios


#===========================
# Plot the power sweep
#===========================
library(ggplot2)
library(scales)

ggplot(res$sweep, aes(x = weeks_per_period, y = power_estimate)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = res$sweep$weeks_per_period) +
  labs(
    x     = "Weeks per Period",
    y     = "Estimated Power",
    title = "Power vs. Period Length in Clustered CRT"
  ) +
  theme_minimal(base_size = 14)
