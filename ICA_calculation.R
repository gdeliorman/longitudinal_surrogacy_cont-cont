## ICA calculation for three models: Galecki, endpoint specific random intercept and shared random intercept

##Packages
## It requires Surrogate version 3.4.4 or later.

install.packages("Surrogate")
library(Surrogate)
packageVersion("Surrogate")


# If the default installation fails to provide Surrogate version 3.4.4 or later, 
# install the package from source and then load it.

# install.packages(
#   "Surrogate",
#   repos = "https://cran.r-project.org",
#   type = "source"
# )
# 
# library(Surrogate)

# Set this to the full path of the folder where you extracted the ZIP file.
OUTPUT_DIR <- "path/to/extracted/folder"


##read model fit files RDS files
fit_galecki<- readRDS(file.path(OUTPUT_DIR, "fit_galecki.rds"))
fit_ri<- readRDS(file.path(OUTPUT_DIR, "fit_ri.rds"))
fit_rislope<- readRDS(file.path(OUTPUT_DIR, "fit_rislope.rds"))
fit_shared_ri<- readRDS(file.path(OUTPUT_DIR, "fit_shared_ri.rds"))
fit_shared_rislope<- readRDS(file.path(OUTPUT_DIR, "fit_shared_rislope.rds"))

##GALECKI ICA
fit_galecki$model <- "galecki"

ICA_galecki_05<- Surrogate::longitudinal_ica( object= fit_galecki,
                                   evaluation_times=c(1,2,4,6,8),
                                   sensitivity = list(residual = seq(-1, 1, by = 0.05)),
                                   extrapolation = "warn",
                                   return_candidates = FALSE,
                                   progress = FALSE,
                                   tol = 1e-12)


##Example Error: vector memory limit of 16.0 Gb reached, see mem.maxVSize()

# ICA_galecki_01<- Surrogate::longitudinal_ica( object= fit_galecki,
#                                               evaluation_times=c(1,2,4,6,8),
#                                               sensitivity = list(residual = seq(-1, 1, by = 0.01)),
#                                               extrapolation = "warn",
#                                               return_candidates = FALSE,
#                                               progress = FALSE,
#                                               tol = 1e-12)


## Endpoint specific random intercept ICA
fit_ri$model <- "endpoint_ri"

##it takes time, file is 7 GB!

# ICA_ri<- Surrogate::longitudinal_ica(object= fit_ri,
#                           evaluation_times=c(1,2,4,6,8),
#                           sensitivity = list(random_effects = seq(-1, 1, by = 0.05),
#                                              residual = seq(-1, 1, by = 0.05)),
#                           extrapolation = "warn",
#                           return_candidates = FALSE,
#                           progress = FALSE,
#                           tol = 1e-12)


ICA_ri_01<- Surrogate::longitudinal_ica(object= fit_ri,
                                     evaluation_times=c(1,2,4,6,8),
                                     sensitivity = list(random_effects = seq(-1, 1, by = 0.1),
                                                        residual = seq(-1, 1, by = 0.1)),
                                     extrapolation = "warn",
                                     return_candidates = FALSE,
                                     progress = FALSE,
                                     tol = 1e-12)


##shared random intercept ICA

fit_shared_ri$model <- "shared_ri"

ICA_shared_ri<- Surrogate::longitudinal_ica(object= fit_shared_ri,
                                            evaluation_times=c(1,2,4,6,8),
                                            sensitivity = list(random_effects = seq(-1, 1, by = 0.001),
                                                               residual = seq(-1, 1, by = 0.1)),
                                            extrapolation = "warn",
                                            return_candidates = FALSE,
                                            progress = FALSE,
                                            tol = 1e-12)


## ExampleError: vector memory limit of 16.0 Gb reached, see mem.maxVSize() when residual 0.001
# ICA_shared_ri<- Surrogate::longitudinal_ica(object= fit_shared_ri,
#                                  evaluation_times=c(1,2,4,6,8),
#                                  sensitivity = list(random_effects = seq(-1, 1, by = 0.001),
#                                                     residual = seq(-1, 1, by = 0.01)),
#                                  extrapolation = "warn",
#                                  return_candidates = FALSE,
#                                  progress = FALSE,
#                                  tol = 1e-12)





##summarized the results
summary_longitudinal_ica(ICA_galecki_05)
summary_longitudinal_ica(ICA_ri_01)
summary_longitudinal_ica(ICA_shared_ri)

##summarized the positive-definite completions 
counts <- function(x) {
  with(x, c(V = Num.Pos.Def.V, D = Num.Pos.Def.D, Total = Num.Pos.Def.Pairs))
}

rbind(
  Galecki = counts(ICA_galecki_05),
  `Endpoint-specific RI` = counts(ICA_ri_01),
  `Shared RI` = counts(ICA_shared_ri)
)


##Full ECDF and ECDF zoom plots

ica_list <- list(
  "Galecki" = ICA_galecki_05$R2_Lambda,
  "Endpoint-specific RI" = ICA_ri_01$R2_Lambda,
  "Shared RI" = ICA_shared_ri$R2_Lambda
)
cols <- c("Galecki" = "#0072B2",
          "Endpoint-specific RI" = "#D55E00",
          "Shared RI" = "#009E73")

ltys <- c(
  "Galecki" = "dashed",
  "Endpoint-specific RI" = "solid",
  "Shared RI" = "13"
)

thresholds <- c(0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95)
lwd <- 2.5

## Full ECDF color
x_grid <- sort(unique(c(
  seq(0, 1, length.out = 1500),
  1 - 10^seq(0, -8, length.out = 1500)
)))

plot(NA, type = "n",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ECDF of ICA",
     xlab = expression(R[Lambda]^2),
     ylab = "ECDF")

for (k in seq_along(ica_list)) {
  x <- sort(ica_list[[k]][is.finite(ica_list[[k]])])
  y <- findInterval(x_grid, x) / length(x)

  lines(x_grid, y, type = "s",
        col = cols[names(ica_list)[k]],
        lwd = 2.8,
        lty = ltys[names(ica_list)[k]])
}

abline(v = thresholds, col = "grey85", lty = 3)

legend("topleft",
       legend = names(ica_list),
       col = cols[names(ica_list)],
       lwd = 2.3,
       seg.len = 4,
       lty = ltys[names(ica_list)],
       bty = "o",
       bg = "white",
       cex = 0.85)



## Lower-tail zoom
zoom_xlim <- c(0.02, 0.71)
zoom_grid <- sort(unique(c(
  seq(zoom_xlim[1], zoom_xlim[2], length.out = 2500),
  thresholds[thresholds >= zoom_xlim[1] &
               thresholds <= zoom_xlim[2]]
)))

galecki <- ica_list[["Galecki"]]
galecki <- galecki[is.finite(galecki)]

plot(ecdf(galecki),
     do.points = FALSE,
     verticals = TRUE,
     col = "#0072B2",
     lwd = lwd,
     lty = 2,
     main = "Lower-tail zoom",
     xlab = expression(R[Lambda]^2),
     ylab = "ECDF",
     xlim = zoom_xlim,
     ylim = c(0, 0.005))


for (nm in c("Endpoint-specific RI", "Shared RI")) {
  x <- sort(ica_list[[nm]][is.finite(ica_list[[nm]])])
  y <- findInterval(zoom_grid, x) / length(x)

  lines(zoom_grid, y,
        col = cols[nm],
        lwd = lwd,
        lty = ltys[nm])
}

abline(v = thresholds, col = "grey85", lty = 3)


legend("topleft",
       legend = names(ica_list),
       col = cols[names(ica_list)],
       lwd = 2.3,
       seg.len = 4,
       lty = ltys[names(ica_list)],
       bty = "o",
       bg = "white",
       cex = 0.85)


##PLOTS BLACK AND WHITE

## Full ECDF
plot(NA, type = "n",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "ECDF of ICA",
     xlab = expression(R[Lambda]^2),
     ylab = "ECDF")

for (k in seq_along(ica_list)) {
  x <- sort(ica_list[[k]][is.finite(ica_list[[k]])])
  y <- findInterval(x_grid, x) / length(x)

  lines(x_grid, y, type = "s",
        lwd = 2.5,
        lty = ltys[names(ica_list)[k]])
}

abline(v = thresholds, col = "grey85", lty = 3)

legend("topleft",
       legend = names(ica_list),
       lwd = 2,
       seg.len = 3.3,
       lty = ltys[names(ica_list)],
       bty = "o",
       bg = "white",
       cex = 0.85)



## Lower-tail zoom black and white
plot(ecdf(galecki),
     do.points = FALSE,
     verticals = TRUE,
     lwd = lwd,
     lty = "dashed",
     main = "Lower-tail zoom",
     xlab = expression(R[Lambda]^2),
     ylab = "ECDF",
     xlim = zoom_xlim,
     ylim = c(0, 0.005))


for (nm in c("Endpoint-specific RI", "Shared RI")) {
  x <- sort(ica_list[[nm]][is.finite(ica_list[[nm]])])
  y <- findInterval(zoom_grid, x) / length(x)

  lines(zoom_grid, y,
        lwd = lwd,
        lty = ltys[nm])
}

abline(v = thresholds, col = "grey85", lty = 3)

legend("topleft",
       legend = names(ica_list),
       lwd = 2,
       seg.len = 3.3,
       lty = ltys[names(ica_list)],
       bty = "o",
       bg = "white",
       cex = 0.85)



tab <- t(sapply(ica_list, function(x) {
  sapply(thresholds, function(t) mean(x < t))
}))

colnames(tab) <- sprintf("<%.2f", thresholds)

round(tab, 4)

