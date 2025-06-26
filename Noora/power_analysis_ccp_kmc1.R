#==========================================
# One-Step Cluster Crossover Design Power Calculator
#  (2w control → 2w washout → 2w intervention, simutaneously happening within all selected hospitals)
#==========================================
cluster_power <- function(
    p1,                     # baseline proportion (0–1)
    p2,                     # treatment proportion (0–1)
    clusters,               # number of hospitals (clusters)
    births_per_week,        # avg births per hospital per week
    weeks_control      = 2, # length of control period (weeks)
    weeks_washout      = 2, # length of washout period (weeks; data excluded)
    weeks_intervention = 2, # length of intervention period (weeks)
    icc                = 0.02,  # intra-cluster correlation
    alpha              = 0.05,  # two-sided significance level
    use_t              = TRUE   # use t‐crit when clusters < 10
) {
  # 1. Input validation
  if (!is.numeric(p1) || p1 < 0 || p1 > 1) stop("p1 must be between 0 and 1")
  if (!is.numeric(p2) || p2 < 0 || p2 > 1) stop("p2 must be between 0 and 1")
  if (clusters < 2)                   stop("Need at least 2 clusters")
  if (!is.numeric(icc) || icc < 0 || icc > 1) stop("icc must be between 0 and 1")
  
  # 2. Cluster‐period sizes (washout dropped)
  m0    <- weeks_control      * births_per_week
  m1    <- weeks_intervention * births_per_week
  m_avg <- (m0 + m1) / 2       # average cluster‐period size
  
  # 3. Hussey & Hughes variance formula for T=2 active periods (i.e., one treatment and one control period)
  T_active <- 2
  # design effect
  DE       <- 1 + (m_avg - 1) * icc
  # outcome variance under H₀
  sigma2   <- p1 * (1 - p1)
  # Var(β̂) = [σ² · DE] / [I · m_avg · T_active] × [T_active/(T_active−1)]
  var_beta <- sigma2 * DE / (clusters * m_avg * T_active) * (T_active / (T_active - 1))
  se_beta  <- sqrt(var_beta)
  
  # 4. Critical value (Normal vs. t)
  if (use_t && clusters < 10) {
    df   <- clusters - 1
    crit <- qt(1 - alpha/2, df = df)
  } else {
    crit <- qnorm(1 - alpha/2)
  }
  
  # 5. Power
  delta <- abs(p2 - p1)
  power <- pnorm(delta / se_beta - crit)
  
  # 6. Return
  list(
    p1                = p1,
    p2                = p2,
    clusters          = clusters,
    weeks_control     = weeks_control,
    weeks_washout     = weeks_washout,
    weeks_intervention= weeks_intervention,
    births_per_week   = births_per_week,
    icc               = icc,
    design_effect     = DE,
    sigma2            = sigma2,
    T_active          = T_active,
    var_beta          = var_beta,
    se_estimate       = se_beta,
    critical_value    = crit,
    delta             = delta,
    power             = power
  )
}


#===========================
# Single power estimate
#===========================

# average birth per week across the 8 participating hospitals in Karnataka
births_per_week_karnataka = mean(c(259.78, 169.56, 145.44, 133.33,
                                   125.00, 114.11, 102.99, 101.11)) / 31 * 7

res <- cluster_power(
  p1               = 0.081, # default KMC practice rate from Noora survey (Subramanian et al., 2020)
  p2               = 0.181, # 10% point expected increase due to the social norms intervention
  clusters         = 8, # 8 hospitals in karnataka
  births_per_week  = births_per_week_karnataka,
  weeks_control      = 2,
  weeks_washout      = 2,
  weeks_intervention = 2,
  icc                = 0.02,
  alpha              = 0.05,
  use_t              = TRUE
)

print(res$power)

#===========================
# Plot different scenarios
#===========================

# Set parameters & build grid of scenarios
p1              <- 0.31
clusters_vec    <- 2:8         # 2 through 8 hospitals
deltas          <- seq(0.10, 0.30, by = 0.05)  # +10pp to +30pp

# create data frame of all combinations
df <- expand.grid(
  clusters = clusters_vec,
  delta    = deltas
)

# Compute power for each row by calling sw_cluster_power()
df$power <- mapply(
  FUN = function(k, d) {
    # call the existing function, extract the 'power' element
    cluster_power(
      p1                = p1,
      p2                = p1 + d,
      clusters          = k,
      births_per_week   = births_per_week_karnataka,
      weeks_control     = 2,
      weeks_washout     = 2,
      weeks_intervention= 2,
      icc               = 0.02,
      alpha             = 0.05,
      use_t             = TRUE
    )$power
  },
  df$clusters,
  df$delta
)

# Plot with ggplot2
ggplot(df, aes(x = clusters, y = power, color = factor(delta))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_brewer(
    palette = "Dark2",
    name    = "Effect size\n(Δ pp)",
    labels  = paste0("+", deltas * 100, "pp")
  ) +
  scale_x_continuous(breaks = clusters_vec) +
  labs(
    x        = "Number of Hospitals (Clusters)",
    y        = "Estimated Power",
    title    = "Power vs. # of Karnataka Hospitals for Various Effect Sizes",
    subtitle = "2w control → 2w washout → 2w intervention\n(Baseline KMC practice rate = 8.1%, ICC = 0.02, births/week = 32)"
  ) +
  ylim(0, 1) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
