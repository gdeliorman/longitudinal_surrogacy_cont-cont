## Biometrics revision - Simulation 1: Population-level numerical calculation
##   - No data are generated.
##   - No longitudinal model is fitted.


library(ggplot2) ##for plots
library(matrixcalc) ##for pdf marix

# out_dir <- "/Users/gokcedeliorman/Desktop/PHD/Third_paper_multivariate_surrogates/after_revision"
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

p_grid <- 1:20
rho_R <- 0.5 
rho_u_values <- c(0, 0.1, 0.3, 0.5, 0.7)


## ICA general formula using determinant definition
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

## Galecki closed form
r2_galecki_closed <- function(p, rho_u) {
  1 - (1 - rho_u^2)^p
}


## First sim:  Galecki simulation
galecki_df <- data.frame()

for (rho_u in rho_u_values) {
  for (p in p_grid) {
    
    ## AR(1)-type temporal correlation matrix
    Rp <- rho_R ^ abs(outer(1:p, 1:p, "-"))
    
    ## Sigma_u matrix
    Sigma_u <- matrix(
      c(1, rho_u,
        rho_u, 1),
      nrow = 2,
      byrow = TRUE
    )
    
    ## Galecki model:
    ## Sigma_delta = Sigma_u %x% Rp
    Sigma_delta <- kronecker(Sigma_u, Rp)
    
    ## ICA from determinant definition
    r2_general <- r2_lambda_from_cov(Sigma_delta, p)
    
    ## ICA from closed form
    r2_closed <- r2_galecki_closed(p, rho_u)
    
    galecki_df <- rbind(
      galecki_df,
      data.frame(
        p = p,
        rho_u = rho_u,
        R2_general = r2_general,
        R2_closed = r2_closed,
        difference = r2_general - r2_closed
      )
    )
  }
}

galecki_df
# write.csv( galecki_df,file = file.path(out_dir, "results_simulation1_galecki.csv"),
# row.names = FALSE)

ggplot(galecki_df, aes(x = p, y = R2_closed, color = factor(rho_u))) +
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
  theme(
    legend.position = "bottom"
  )

##Second simulation: Random intercept model

## Random-intercept closed form
r2_ri_closed <- function(p, rho_R, Sigma_u, Sigma_D) {
  
  u11 <- Sigma_u[1, 1]
  u12 <- Sigma_u[1, 2]
  u22 <- Sigma_u[2, 2]
  
  d11 <- Sigma_D[1, 1]
  d12 <- Sigma_D[1, 2]
  d22 <- Sigma_D[2, 2]
  
  ## alpha_p = 1_p' R_p^{-1} 1_p
  alpha_p <- (p - p * rho_R + 2 * rho_R) / (1 + rho_R)
  
  rho_u2 <- u12^2 / (u11 * u22)
  
  rho_ud2 <- (u12 + alpha_p * d12)^2 /
    ((u11 + alpha_p * d11) * (u22 + alpha_p * d22))
  
  r2 <- 1 - (1 - rho_u2)^(p - 1) * (1 - rho_ud2)
  
  return(r2)
}


## Three random-intercept scenarios
ri_scenarios <- data.frame(
  scenario = c(
    "Persistent association only",
    "Occasion-specific association only",
    "Both"
  ),
  u12 = c(0, 0.3, 0.3),
  d12 = c(0.7, 0, 0.7)
)

## Random-intercept simulation
ri_df <- data.frame()

for (i in 1:nrow(ri_scenarios)) {
  
  scenario_name <- ri_scenarios$scenario[i]
  u12 <- ri_scenarios$u12[i]
  d12 <- ri_scenarios$d12[i]
  
  for (p in p_grid) {
    
    ## AR(1)-type temporal correlation matrix
    Rp <- rho_R ^ abs(outer(1:p, 1:p, "-"))
    
    ## Random-intercept model settings:
    ## u11 = u22 = d11 = d22 = 1
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
    
    one_p <- rep(1, p)
    Jp <- one_p %*% t(one_p)
    
    ## Sigma_delta = Sigma_u %x% Rp + Sigma_D %x% Jp
    Sigma_delta <- kronecker(Sigma_u, Rp) + kronecker(Sigma_D, Jp)
    
    ## ICA from determinant definition
    r2_general <- r2_lambda_from_cov(Sigma_delta, p)
    
    ## ICA from random-intercept closed form
    r2_closed <- r2_ri_closed(
      p = p,
      rho_R = rho_R,
      Sigma_u = Sigma_u,
      Sigma_D = Sigma_D
    )
    
    ri_df <- rbind(
      ri_df,
      data.frame(
        scenario = scenario_name,
        p = p,
        u12 = u12,
        d12 = d12,
        R2_general = r2_general,
        R2_closed = r2_closed,
        difference = r2_general - r2_closed
      )
    )
  }
}

ri_df
# write.csv( ri_df,file = file.path(out_dir, "results_simulation1_random_intercept.csv"),
#       row.names = FALSE)

ggplot(ri_df, aes(x = p, y = R2_closed, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point() +
  geom_vline(xintercept = 5, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Number of measurement occasions",
    y = expression(R[Lambda]^2),
    color = "Scenario"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )
