##fit the model

install.packages("Surrogate")
library(Surrogate)

data(Schizo_PANSS_BPRS)
dat_long<-Schizo_PANSS_BPRS

# Set this to the full path of the folder where you extracted the ZIP file.
OUTPUT_DIR <- "path/to/extracted/folder"

#OUTPUT_DIR <- "/Users/gokcedeliorman/Desktop/PHD/Third_paper_multivariate_surrogates/3rd_paper_new_analysis_ariel_package/fit_outputs_for_package"

dat <- Surrogate::prepare_longitudinal_surrogate_data(
  data = dat_long,
  id = "Id",
  treatment = "Treat",
  time = "Time_cont",
  endpoint = "Endpoint",
  outcome = "Response_log",
  surrogate_label = "BPRS",
  true_label = "PANSS",
  treatment_levels = c(0, 1))


mean_formula <- Y ~ 0 + endpoint:visitF:treatF

spec_galecki <- Surrogate::make_longitudinal_re_none(name = "Galecki's marginal model")

spec_ri <- Surrogate::make_longitudinal_re_endpoint_specific(
  name = "M2: endpoint-specific random intercepts",
  S_0 = "ri",
  T_0 = "ri",
  S_1 = "ri",
  T_1 = "ri"
)

## M3: endpoint-specific random intercepts and slopes, AR(1)
spec_rislope <- Surrogate::make_longitudinal_re_endpoint_specific(
  name = "M3: endpoint-specific random intercepts and slopes, AR1",
  S_0 = "rislope",
  T_0 = "rislope",
  S_1 = "rislope",
  T_1 = "rislope"
)


## M4: shared random intercept within treatment arm
spec_shared_ri <- Surrogate::make_longitudinal_re_shared_by_arm(
  name = "M4: shared random intercept within arm",
  arm0 = "ri",
  arm1 = "ri"
)

## M5: shared random intercept and slope within treatment arm
spec_shared_rislope <- Surrogate::make_longitudinal_re_shared_by_arm(
  name = "M5: shared random intercept and slope within arm",
  arm0 = "rislope",
  arm1 = "rislope"
)



start_time<- Sys.time()

#M1
fit_galecki <- Surrogate::fit_longitudinal_surrogate_model(dat,
                                                mean_formula,
                                                model_spec = spec_galecki,
                                                residual_type = "exponential")


saveRDS(fit_galecki, file.path(OUTPUT_DIR, "fit_galecki.rds"))

## M2
fit_ri <- Surrogate::fit_longitudinal_surrogate_model(dat,
                                           mean_formula,
                                           model_spec = spec_ri,
                                           residual_type = "exponential",
                                           start_fit = fit_galecki
)

saveRDS(fit_ri, file.path(OUTPUT_DIR, "fit_ri.rds"))

## M3
fit_rislope <- Surrogate::fit_longitudinal_surrogate_model(dat,
                                                mean_formula,
                                                model_spec = spec_rislope,
                                                residual_type = "exponential",
                                                start_fit = fit_ri
)

saveRDS(fit_rislope, file.path(OUTPUT_DIR, "fit_rislope.rds"))


## M4
fit_shared_ri <- Surrogate::fit_longitudinal_surrogate_model(dat,
                                                  mean_formula,
                                                  model_spec = spec_shared_ri,
                                                  residual_type = "exponential",
                                                  start_fit = fit_galecki
)

saveRDS(fit_shared_ri, file.path(OUTPUT_DIR, "fit_shared_ri.rds"))


## M5
fit_shared_rislope <- Surrogate::fit_longitudinal_surrogate_model(dat,
                                                       mean_formula,
                                                       model_spec = spec_shared_rislope,
                                                       residual_type = "exponential",
                                                       start_fit = fit_shared_ri
)

saveRDS(fit_shared_rislope, file.path(OUTPUT_DIR, "fit_shared_rislope.rds"))



fits_all <- list(
  fit_galecki,
  fit_ri,
  fit_rislope,
  fit_shared_ri,
  fit_shared_rislope
)

end_time<- Sys.time()


##model comparison
fit_stats <- do.call(rbind, lapply(fits_all, Surrogate::longitudinal_fit_stats))
fit_stats

write.csv(
  fit_stats,
  file.path(OUTPUT_DIR, "fit_statistics_all_models.csv"),
  row.names = FALSE
)

##check the matrix diagnostics
diagnostics <- do.call(rbind, lapply(fits_all, Surrogate::longitudinal_matrix_diagnostics))
diagnostics


write.csv(
  diagnostics,
  file.path(OUTPUT_DIR, "matrix_diagnostics_all_models.csv"),
  row.names = FALSE
)




