##Individual plots

# Load required packages
install.packages("Surrogate")
library(Surrogate) 
library(dplyr)

# Set this to the full path of the folder where you extracted the ZIP file.
OUTPUT_DIR <- "path/to/extracted/folder"

##Longitudinal Schizophrenia data
data("Schizo_PANSS_BPRS")

plot_individual_profiles <- function(data, endpoint, z, ylim) {
  d <- data[data$Treat == z & data$Endpoint == endpoint, , drop = FALSE]
  d <- d[order(d$Id, d$Time_cont), , drop = FALSE]
  
  # pdf(file.path(OUTPUT_DIR, filename), width = 8, height = 6)
  # on.exit(dev.off(), add = TRUE)
  
  plot(NA_real_, NA_real_,
       xlim = c(0, 8), ylim = ylim, type = "n", xaxt = "n",
       xlab = "Time", ylab = paste("log", endpoint),
       main = paste(endpoint, "Z =", z, "mean profile"))
  axis(1, at = c(0, 1, 2, 4, 6, 8))
  
  for (id_i in unique(d$Id)) {
    subject <- d[d$Id == id_i, , drop = FALSE]
    subject <- subject[order(subject$Time_cont), , drop = FALSE]
    baseline <- subject$Baseline_log[!is.na(subject$Baseline_log)][1]
    
    lines(c(0, subject$Time_cont),
          c(baseline, subject$Response_log),
          col = "hotpink")
  }
  
  baseline_mean <- mean(
    unique(d[, c("Id", "Baseline_log")])$Baseline_log,
    na.rm = TRUE
  )
  visits <- c(1, 2, 4, 6, 8)
  mean_y <- c(
    baseline_mean,
    vapply(visits, function(t) {
      mean(d$Response_log[d$Time_cont == t], na.rm = TRUE)
    }, numeric(1))
  )
  mean_x <- c(0, visits)
  
  points(mean_x, mean_y, pch = 16)
  lines(mean_x, mean_y, col = "black", lwd = 2)
}

plot_individual_profiles(Schizo_PANSS_BPRS, "PANSS", 0, c(3.4, 5.1))
plot_individual_profiles(Schizo_PANSS_BPRS, "PANSS", 1, c(3.4, 5.1))
plot_individual_profiles(Schizo_PANSS_BPRS, "BPRS", 0, c(2.7, 4.8))
plot_individual_profiles(Schizo_PANSS_BPRS, "BPRS", 1, c(2.7, 4.6))

