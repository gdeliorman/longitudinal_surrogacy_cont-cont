## Biometrics revision 
## Additional analysis 
## New numerical example with increasingly dense grids on fixed interval [1, 20].


#rm(list = ls())

library(ggplot2)    ##for plot
library(matrixcalc) ##for positive definite matrix

##You can write the results your own path
# out_dir <- "/Users/gokcedeliorman/Desktop/PHD/Third_paper_multivariate_surrogates/after_revision/simulation_2"
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

p_grid_fixed <- 2:20
p_grid_old <- 1:20 

gamma <- -log(0.5)
rho_u_values <- c(0, 0.1, 0.3, 0.5, 0.7)

## ICA general formula
r2_lambda_from_cov <- function(Sigma_delta, p, tol = 1e-10) {
  
  Sigma_delta <- as.matrix(Sigma_delta)
  
  Sigma_T <- Sigma_delta[1:p, 1:p, drop = FALSE]
  Sigma_S <- Sigma_delta[(p + 1):(2 * p), (p + 1):(2 * p), drop = FALSE]
  
  if (!is.positive.definite(Sigma_delta, tol = tol)) {
    stop("Full Sigma_delta is not positive definite.")
  }
  if (!is.positive.definite(Sigma_T, tol = tol)) {
    stop("Sigma_T is not positive definite.")
  }
  if (!is.positive.definite(Sigma_S, tol = tol)) {
    stop("Sigma_S is not positive definite.")
  }
  
  lambda <- det(Sigma_delta) / (det(Sigma_T) * det(Sigma_S))
  r2 <- 1 - lambda
  
  return(r2)
}

## Galecki model closed form
r2_galecki_closed <- function(p, rho_u) {
  1 - (1 - rho_u^2)^p
}

## Random intercept closed form
r2_ri_closed <- function(p, alpha_p, Sigma_u, Sigma_D) {
  
  u11 <- Sigma_u[1, 1]
  u12 <- Sigma_u[1, 2]
  u22 <- Sigma_u[2, 2]
  
  d11 <- Sigma_D[1, 1]
  d12 <- Sigma_D[1, 2]
  d22 <- Sigma_D[2, 2]
  
  rho_u2 <- u12^2 / (u11 * u22)
  
  rho_ud2 <- (u12 + alpha_p * d12)^2 /
    ((u11 + alpha_p * d11) * (u22 + alpha_p * d22))
  
  r2 <- 1 - (1 - rho_u2)^(p - 1) * (1 - rho_ud2)
  
  return(r2)
}


##galecki model simulation 1 / old version

galecki_old_df <- data.frame()

for (rho_u in rho_u_values) {
  for (p in p_grid_old) {
    
    time_grid <- 1:p
    Rp <- exp(-gamma * abs(outer(time_grid, time_grid, "-")))
    
    Sigma_u <- matrix(
      c(1, rho_u,
        rho_u, 1),
      nrow = 2,
      byrow = TRUE
    )
    
    Sigma_delta <- kronecker(Sigma_u, Rp)
    
    R2_general <- r2_lambda_from_cov(Sigma_delta, p)
    R2_closed <- r2_galecki_closed(p, rho_u)
    
    galecki_old_df <- rbind(
      galecki_old_df,
      data.frame(
        grid_type = "Expanding follow-up",
        p = p,
        delta_p = 1,
        rho_adj = ifelse(p == 1, NA, 0.5),
        rho_u = rho_u,
        R2_general = R2_general,
        R2_closed = R2_closed,
        difference = R2_general - R2_closed
      )
    )
  }
}

galecki_old_df <- galecki_old_df[order(galecki_old_df$p, galecki_old_df$rho_u), ]
rownames(galecki_old_df) <- NULL

# write.csv(
#   galecki_old_df,
#   file = file.path(out_dir, "results_expanding_followup_galecki.csv"),
#   row.names = FALSE
# )


##galecki model new simulation
galecki_fixed_df <- data.frame()

for (rho_u in rho_u_values) {
  for (p in p_grid_fixed) {
    
    delta_p <- 19 / (p - 1)
    time_grid <- 1 + (0:(p - 1)) * delta_p
    
    Rp <- exp(-gamma * abs(outer(time_grid, time_grid, "-")))
    rho_adj <- exp(-gamma * delta_p)
    
    Sigma_u <- matrix(
      c(1, rho_u,
        rho_u, 1),
      nrow = 2,
      byrow = TRUE
    )
    
    Sigma_delta <- kronecker(Sigma_u, Rp)
    
    R2_general <- r2_lambda_from_cov(Sigma_delta, p)
    R2_closed <- r2_galecki_closed(p, rho_u)
    
    galecki_fixed_df <- rbind(
      galecki_fixed_df,
      data.frame(
        p = p,
        delta_p = delta_p,
        rho_adj = rho_adj,
        rho_u = rho_u,
        R2_general = R2_general,
        R2_closed = R2_closed,
        difference = R2_general - R2_closed
      )
    )
  }
}

galecki_fixed_df <- galecki_fixed_df[order(galecki_fixed_df$p, galecki_fixed_df$rho_u), ]
rownames(galecki_fixed_df) <- NULL


# write.csv(
#   galecki_fixed_df,
#   file = file.path(out_dir, "results_fixed_interval_galecki.csv"),
#   row.names = FALSE
# )


##the results are the same 
galecki_compare <- merge(
  galecki_fixed_df,
  galecki_old_df[, c("p", "rho_u", "R2_closed")],
  by = c("p", "rho_u"),
  suffixes = c("_fixed", "_old")
)


###random intercept model 
ri_scenarios <- data.frame(
  scenario = c(
    "Persistent association only",
    "Occasion-specific association only",
    "Both"
  ),
  u12 = c(0, 0.3, 0.3),
  d12 = c(0.7, 0, 0.7)
)

ri_old_df <- data.frame()
ri_fixed_df <- data.frame()

for (i in 1:nrow(ri_scenarios)) {
  
  scenario_name <- ri_scenarios$scenario[i]
  u12 <- ri_scenarios$u12[i]
  d12 <- ri_scenarios$d12[i]
  
  for (p in p_grid_fixed) {
    
    ## Old expanding grid
    time_grid_old <- 1:p
    Rp_old <- exp(-gamma * abs(outer(time_grid_old, time_grid_old, "-")))
    one_p <- rep(1, p)
    Jp <- one_p %*% t(one_p)
    
    alpha_direct_old <- as.numeric(t(one_p) %*% solve(Rp_old, one_p))
    alpha_closed_old <- (p - p * 0.5 + 2 * 0.5) / (1 + 0.5)
    
    Sigma_u <- matrix(
      c(1, u12,
        u12, 1),
      nrow = 2,
      byrow = TRUE
    )
    
    Sigma_D <- matrix(
      c(1, d12,
        d12, 1),
      nrow = 2,
      byrow = TRUE
    )
    
    Sigma_delta_old <- kronecker(Sigma_u, Rp_old) + kronecker(Sigma_D, Jp)
    
    R2_general_old <- r2_lambda_from_cov(Sigma_delta_old, p)
    R2_closed_old <- r2_ri_closed(
      p = p,
      alpha_p = alpha_direct_old,
      Sigma_u = Sigma_u,
      Sigma_D = Sigma_D
    )
    
    ri_old_df <- rbind(
      ri_old_df,
      data.frame(
        grid_type = "Expanding follow-up",
        scenario = scenario_name,
        p = p,
        delta_p = 1,
        rho_adj = 0.5,
        alpha_direct = alpha_direct_old,
        alpha_closed = alpha_closed_old,
        u12 = u12,
        d12 = d12,
        R2_general = R2_general_old,
        R2_closed = R2_closed_old,
        difference = R2_general_old - R2_closed_old
      )
    )
    
    ## New fixed interval [1, 20]
    delta_p <- 19 / (p - 1)
    time_grid_fixed <- 1 + (0:(p - 1)) * delta_p
    
    Rp_fixed <- exp(-gamma * abs(outer(time_grid_fixed, time_grid_fixed, "-")))
    rho_adj <- exp(-gamma * delta_p)
    
    alpha_direct <- as.numeric(t(one_p) %*% solve(Rp_fixed, one_p))
    alpha_closed <- (p - p * rho_adj + 2 * rho_adj) / (1 + rho_adj)
    
    Sigma_delta_fixed <- kronecker(Sigma_u, Rp_fixed) + kronecker(Sigma_D, Jp)
    
    R2_general <- r2_lambda_from_cov(Sigma_delta_fixed, p)
    R2_closed <- r2_ri_closed(
      p = p,
      alpha_p = alpha_direct,
      Sigma_u = Sigma_u,
      Sigma_D = Sigma_D
    )
    
    ri_fixed_df <- rbind(
      ri_fixed_df,
      data.frame(
        scenario = scenario_name,
        p = p,
        delta_p = delta_p,
        rho_adj = rho_adj,
        alpha_direct = alpha_direct,
        alpha_closed = alpha_closed,
        u12 = u12,
        d12 = d12,
        R2_general = R2_general,
        R2_closed = R2_closed,
        difference = R2_general - R2_closed
      )
    )
  }
}

ri_old_df <- ri_old_df[
  order(ri_old_df$p, match(ri_old_df$scenario, ri_scenarios$scenario)),
]
rownames(ri_old_df) <- NULL

ri_fixed_df <- ri_fixed_df[
  order(ri_fixed_df$p, match(ri_fixed_df$scenario, ri_scenarios$scenario)),
]
rownames(ri_fixed_df) <- NULL

# write.csv(
#   ri_fixed_df,
#   file = file.path(out_dir, "results_fixed_interval_random_intercept.csv"),
#   row.names = FALSE
# )


ri_compare_p20 <- merge(
  ri_fixed_df[ri_fixed_df$p == 20, ],
  ri_old_df[ri_old_df$p == 20, c("scenario", "R2_closed")],
  by = "scenario",
  suffixes = c("_fixed", "_old")
)



## plots

fig_galecki_expanding_followup <- ggplot(galecki_old_df,
                                         aes(x = p, y = R2_closed,
                                             color = factor(rho_u))) +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_vline(xintercept = 5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Number of measurement occasions",
    y = expression(R[Lambda]^2),
    color = expression(rho[u])
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

fig_galecki_expanding_followup

fig_galecki_fixed_interval <- ggplot(galecki_fixed_df,
                                     aes(x = p, y = R2_closed,
                                         color = factor(rho_u))) +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_vline(xintercept = 5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Number of measurement occasions",
    y = expression(R[Lambda]^2),
    color = expression(rho[u])
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

fig_galecki_fixed_interval


##long format

ri_plot_df <- rbind(ri_old_df, data.frame(grid_type = "Fixed interval [1,20]", ri_fixed_df))

ri_persistent_df <- ri_plot_df[
  ri_plot_df$scenario == "Persistent association only",
]

fig_grid_comparison_RI_persistent <- ggplot(ri_persistent_df,
                                            aes(x = p, y = R2_closed,
                                                color = grid_type)) +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_vline(xintercept = 5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Number of measurement occasions",
    y = expression(R[Lambda]^2),
    color = "Grid"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

fig_grid_comparison_RI_persistent

##long format
ri_occasion_df <- ri_plot_df[
  ri_plot_df$scenario == "Occasion-specific association only",
]

fig_grid_comparison_RI_occasion <- ggplot(ri_occasion_df,
                                          aes(x = p, y = R2_closed,
                                              color = grid_type)) +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_vline(xintercept = 5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Number of measurement occasions",
    y = expression(R[Lambda]^2),
    color = "Grid"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

fig_grid_comparison_RI_occasion

##long format
ri_both_df <- ri_plot_df[
  ri_plot_df$scenario == "Both",
]

fig_grid_comparison_RI_both <- ggplot(ri_both_df,
                                      aes(x = p, y = R2_closed,
                                          color = grid_type)) +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_vline(xintercept = 5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Number of measurement occasions",
    y = expression(R[Lambda]^2),
    color = "Grid"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

fig_grid_comparison_RI_both

##to see differences in ri models expanding vs fixed
ri_R2_closed_comparison <- merge(
  ri_fixed_df[, c("scenario", "p", "rho_adj", "alpha_direct", "R2_closed")],
  ri_old_df[, c("scenario", "p", "rho_adj", "alpha_direct", "R2_closed")],
  by = c("scenario", "p"),
  suffixes = c("_fixed", "_expanding")
)

ri_R2_closed_comparison$scenario_p <- paste0(
  ri_R2_closed_comparison$scenario,
  " | p = ",
  ri_R2_closed_comparison$p
)

ri_R2_closed_comparison$R2_closed_difference <-
  ri_R2_closed_comparison$R2_closed_fixed -
  ri_R2_closed_comparison$R2_closed_expanding

ri_R2_closed_comparison$alpha_difference <-
  ri_R2_closed_comparison$alpha_direct_fixed -
  ri_R2_closed_comparison$alpha_direct_expanding

ri_R2_closed_comparison <- ri_R2_closed_comparison[, c(
  "scenario_p",
  "scenario",
  "p",
  "rho_adj_expanding",
  "rho_adj_fixed",
  "alpha_direct_expanding",
  "alpha_direct_fixed",
  "alpha_difference",
  "R2_closed_expanding",
  "R2_closed_fixed",
  "R2_closed_difference"
)]

ri_R2_closed_comparison



##alpha_p vs p 
alpha_old <- unique(ri_old_df[, c("grid_type", "p", "alpha_direct")])
alpha_fixed <- unique(data.frame(
  grid_type = "Fixed interval [1,20]",
  ri_fixed_df[, c("p", "alpha_direct")]
))
alpha_df <- rbind(alpha_old, alpha_fixed)

fig_alpha_grid_comparison <- ggplot(alpha_df,
                                    aes(x = p, y = alpha_direct,
                                        color = grid_type)) +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_vline(xintercept = 5, linetype = "dashed") +
  labs(
    x = "Number of measurement occasions",
    y = expression(alpha[p]),
    color = "Grid"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")


fig_alpha_grid_comparison


##Validations

##all ICA in [0,1] apart from negligible floating-point error
all_R2_values <- c(
  galecki_old_df$R2_general,
  galecki_old_df$R2_closed,
  galecki_fixed_df$R2_general,
  galecki_fixed_df$R2_closed,
  ri_old_df$R2_general,
  ri_old_df$R2_closed,
  ri_fixed_df$R2_general,
  ri_fixed_df$R2_closed
)

all(round(all_R2_values, 10) >= 0 & round(all_R2_values, 10) <= 1)

##differences between closed form and general formula < 10^-10

all_general_closed_differences <- c(
  galecki_old_df$difference,
  galecki_fixed_df$difference,
  ri_old_df$difference,
  ri_fixed_df$difference
)

all(abs(all_general_closed_differences) < 1e-10)

##in Galecki model the results are the same 
galecki_compare$R2_closed_fixed_vs_expanding <- 
  galecki_compare$R2_closed_fixed - galecki_compare$R2_closed_old

all(abs(galecki_compare$R2_closed_fixed_vs_expanding) == 0)

##in RI model 
ri_compare_p20$R2_closed_fixed_vs_expanding <-
  ri_compare_p20$R2_closed_fixed - ri_compare_p20$R2_closed_old

all(abs(ri_compare_p20$R2_closed_fixed_vs_expanding) < 1e-10)


