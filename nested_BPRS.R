##Nested BPRS analysis

# Load required packages
install.packages("Surrogate")
library(Surrogate)
library(dplyr)


# Set this to the full path of the folder where you extracted the ZIP file.
OUTPUT_DIR <- "path/to/extracted/folder"


BPRS_grids <- list(
  week1 = c(1),
  week1_2 = c(1, 2),
  week1_2_4 = c(1, 2, 4),
  week1_2_4_6 = c(1, 2, 4, 6),
  week1_2_4_6_8 = c(1, 2, 4, 6, 8)
)


##NESTED BPRS for random intercept model
fit_ri<- readRDS(file.path(OUTPUT_DIR, "fit_ri.rds"))
est_ri <- extract_longitudinal_parameters(fit_ri)


ica_nested<- ICA_contcont_long_ri_nested(V0=est_ri$V0,
                                         V1=est_ri$V1,
                                         D0=est_ri$D0,
                                         D1=est_ri$D1,
                                         times_T=c(1, 2, 4, 6, 8),
                                         times_S=BPRS_grids,
                                         temporal_rho=est_ri$rho,
                                         grid_V = seq(-1, 1, by = 0.05),
                                         grid_D = seq(-1, 1, by = 0.05),
                                         tol = 1e-12,
                                         progress = FALSE)

weeks <- colnames(ica_nested$R2_Lambda)
summaries <- lapply(weeks, function(weeks) {
  summary_longitudinal_ica(
    list(R2_Lambda = ica_nested$R2_Lambda[, weeks]),
    threshold = 0.95
  )
})

names(summaries) <- weeks
do.call(rbind, summaries)
