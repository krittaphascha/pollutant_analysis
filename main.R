set.seed(123)
setwd("~/Projects/hf_pm_analysis")

dat <- read.csv("data/daily_all_weather_30jul.csv")

# Load necessary packages

library(dplyr)
library(tibble)
library(dlnm)
library(splines)
library(tsModel)
library(gnm)
library(tidyverse)
library(xtable)

## SELECT POLLUTANT OF INTEREST
pollutant_col_name <- "PM2.5"
outcome_col_name <- "hf_prim"

lag_number <- 7
output_dir <- "30jul"

# ===========================
# properly format date
dat$date_start <- as.Date(dat$date_start, "%Y-%m-%d")
dat$dow <- as.factor(weekdays(dat$date_start))
dat$month <- as.factor(months(dat$date_start))
dat$year <- as.factor(format(dat$date_start, format = "%Y"))

# get strata for case crossover analysis
dat$stratum <- as.factor(dat$year:dat$month:dat$dow)

prov_name <- as.character(unique(dat$prov_name))

# separate dataframe by provinces
datalist <- lapply(
  provinces,
  function(provname) dat[dat$prov_name == provname, ]
)
names(datalist) <- prov_name

# ============================

results_dir <- paste0("output/", output_dir, "/results_", pollutant_col_name, "-", outcome_col_name, "-", lag_number, "lag")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE) # recursive = TRUE if pollutant_col_name could contain slashes (e.g. "PM2.5/O3")

# For proper name format in the plot

if (outcome_col_name == "hf_prim") {
  outcome_plot_name <- "Heart Failure Hospitalization"
} else {
  outcome_plot_name <- outcome_col_name
}


cat("\nSaving results to directory:", file.path(getwd(), results_dir), "\n")

## RANGE FOR METEOROLOGICAL VARIABLES ##
rangepoll <- t(sapply(
  datalist,
  function(x) {
    range(x[[pollutant_col_name]],
      na.rm = TRUE
    )
  }
))

# ## DLNM PARAMETERS
boundpoll <- c(min(rangepoll), max(rangepoll))

## ADDITIONAL INFORMATION
m <- length(datalist)

####################################################
## FIRST STAGE MODEL (PROVINCE-SPECIFIC ESTIMATE) ##
####################################################

## BUILT OBJECTS WHERE RESULTS WILL BE STORED: # ymat IS THE MATRIX FOR THE OUTCOME PARAMETERS
## Slist IS THE LIST WITH (CO)VARIANCE MATRICES # LINEAR FOR NO2 AND 3 DFs FOR ITS LAG
ymat_pollutant <- matrix(NA, length(datalist), 3, dimnames = list(prov_name, paste("spl", seq(3), sep = ""))) # 3DF
list_vcov_pollutant <- vector("list", length(datalist))
names(list_vcov_pollutant) <- prov_name

## FUNCTION TO COMPUTE THE Q-AIC IN QUASI-POISSON MODELS
fqaic <- function(model) {
  loglik <- sum(dpois(model$y, model$fitted.values, log = TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2 * loglik + 2 * summary(model)$df[3] * phi
  return(qaic)
}

## MATRIX FOR Q-AIC VALUE
qaic <- matrix(NA, length(datalist), 1, dimnames = list(prov_name, paste("Q-AIC")))

## WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED: RUN THE FIRST STAGE ANALYSIS
options(warn = -1)

## add COLUMNS that cannot be nulled
required_cols <- c("temperature", "humidity", "pressure")

#### MODEL SELECTION ######

# create function to return sum of QAIC values given model
# model_predictors <- "cb_pollutant + ns(temperature, 3) + ns(humidity, 3) + ns(pressure, 3) + as.factor(is_holiday)"
model_predictors <- "cb_pollutant + ns(temperature, 3) + as.factor(is_holiday)"
model_formula_string <- paste(outcome_col_name, "~", model_predictors)
model_formula_string <- as.formula(model_formula_string)

create_model_formula <- function(predictors_string) {
  model_formula_string <- paste(outcome_col_name, "~", predictors_string)
  return(as.formula(model_formula_string))
}

model_1_string <- "cb_pollutant + ns(temperature, 3) + ns(humidity, 3) + ns(pressure, 3) + as.factor(is_holiday)"
model_2_string <- "cb_pollutant"
model_3_string <- "cb_pollutant + ns(temperature, 3)"
model_4_string <- "cb_pollutant + ns(humidity, 3)"
model_5_string <- "cb_pollutant + ns(pressure, 3)"
model_6_string <- "cb_pollutant + as.factor(is_holiday)"
model_7_string <- paste(model_6_string, "ns(temperature, 3)", sep = " + ")
model_8_string <- paste(model_6_string, "ns(humidity, 3)", sep = " + ")
model_9_string <- paste(model_6_string, "ns(pressure, 3)", sep = " + ")
model_10_string <- paste(model_7_string, "ns(humidity, 3)", sep = " + ")
model_11_string <- paste(model_7_string, "ns(pressure, 3)", sep = " + ")
model_12_string <- paste(model_7_string, "as.factor(is_lotteryday)", sep = " + ")

model_1 <- create_model_formula(model_1_string)
model_2 <- create_model_formula(model_2_string)
model_3 <- create_model_formula(model_3_string)
model_4 <- create_model_formula(model_4_string)
model_5 <- create_model_formula(model_5_string)
model_6 <- create_model_formula(model_6_string)

# model_6 (holidays) gain the lowest QAIC values
model_7 <- create_model_formula(model_7_string)
model_8 <- create_model_formula(model_8_string)
model_9 <- create_model_formula(model_9_string)

# model_7 (temperature) gains the lowest QAIC values
model_10 <- create_model_formula(model_10_string)
model_11 <- create_model_formula(model_11_string)
model_12 <- create_model_formula(model_12_string)

test_model <- function(model_string, datalist) {
  # Initialize a matrix to store results
  tmp_qaic <- matrix(NA, length(datalist), 1, dimnames = list(prov_name, "Q-AIC"))

  # The system.time wrapper can be used here to time the whole loop
  system.time({
    for (i in seq_along(datalist)) { # Using seq_along is safer than seq(m)

      # PRINT ITERATION
      cat(i, "")

      # Get data for this iteration
      original_data <- datalist[[i]]

      # --- Data Cleaning Step ---
      # Use a NEW variable name here to avoid overwriting the input list
      cleaned_data <- original_data[complete.cases(original_data[, required_cols]), ]

      if (nrow(cleaned_data) == 0) {
        warning(paste("Iteration", i, "skipped: No complete cases after removing NAs."))
        next
      }

      options(na.action = "na.exclude")

      # CREATE THE SPLINE
      cb_pollutant <- crossbasis(cleaned_data[[pollutant_col_name]],
        lag = c(0, 7),
        argvar = list(type = "lin", cen = FALSE),
        arglag = list(fun = "ns", df = 3)
      )

      # Add the crossbasis object to the cleaned data
      cleaned_data$cb_pollutant <- cb_pollutant

      # RUN THE CONDITIONAL-POISSON MODEL
      model <- gnm(model_string,
        family = quasipoisson(),
        eliminate = factor(stratum),
        data = cleaned_data
      ) # Use the new variable here

      # Q-AIC COMPUTATION
      tmp_qaic[i, ] <- fqaic(model)

      # EXTRACT AND SAVE THE RELATED COEF AND VCOV
      pred_stage1_pollutant <- crosspred(cb_pollutant, model, cen = FALSE, cumul = TRUE)
      ymat_pollutant[i, ] <- pred_stage1_pollutant$coef
      list_vcov_pollutant[[i]] <- pred_stage1_pollutant$vcov
    } # End of for loop
  }) # End of system.time

  return(sum(tmp_qaic, na.rm = TRUE))
}

# Run the test for each model
model_1_qaic <- test_model(model_1, datalist)
model_2_qaic <- test_model(model_2, datalist)
model_3_qaic <- test_model(model_3, datalist)
model_4_qaic <- test_model(model_4, datalist)
model_5_qaic <- test_model(model_5, datalist)
model_6_qaic <- test_model(model_6, datalist)
model_7_qaic <- test_model(model_7, datalist)
model_8_qaic <- test_model(model_8, datalist)
model_9_qaic <- test_model(model_9, datalist)
model_10_qaic <- test_model(model_10, datalist)
model_11_qaic <- test_model(model_11, datalist)
model_12_qaic <- test_model(model_12, datalist)

print(paste("Model 1 QAIC:", model_1_qaic))
print(paste("Model 2 QAIC:", model_2_qaic))
print(paste("Model 3 QAIC:", model_3_qaic))
print(paste("Model 4 QAIC:", model_4_qaic))
print(paste("Model 5 QAIC:", model_5_qaic))
print(paste("Model 6 QAIC:", model_6_qaic))
print(paste("Model 7 QAIC:", model_7_qaic))
print(paste("Model 8 QAIC:", model_8_qaic))
print(paste("Model 9 QAIC:", model_9_qaic))
print(paste("Model 10 QAIC:", model_10_qaic))
print(paste("Model 11 QAIC:", model_11_qaic))
print(paste("Model 12 QAIC:", model_12_qaic))

# Save all in list
model_qaic_list <- list(
  model_1 = model_1_qaic,
  model_2 = model_2_qaic,
  model_3 = model_3_qaic,
  model_4 = model_4_qaic,
  model_5 = model_5_qaic,
  model_6 = model_6_qaic,
  model_7 = model_7_qaic,
  model_8 = model_8_qaic,
  model_9 = model_9_qaic,
  model_10 = model_10_qaic,
  model_11 = model_11_qaic,
  model_12 = model_12_qaic
)
# Save the model QAIC values to a file
model_qaic_df <- as.data.frame(model_qaic_list)
model_string <- c(
  model_1_string,
  model_2_string,
  model_3_string,
  model_4_string,
  model_5_string,
  model_6_string,
  model_7_string,
  model_8_string,
  model_9_string,
  model_10_string,
  model_11_string,
  model_12_string
)
model_qaic_df <- rbind(model_string, model_qaic_df)

write.csv(model_qaic_df, file.path(results_dir, "model_qaic_values.csv"), row.names = TRUE)

# MODEL7 Gain lowest QAIC values
# cb_pollutant + as.factor(is_holiday) + ns(temperature, 3)

## Run the SELECTED MODEL

# update the required_cols parameters as the model change
# in this case holiday has non-null value
required_cols <- c("temperature")

system.time(
  for (i in seq(m)) {
    # PRINT ITERATION
    cat(i, "")

    # Get data for this iteration
    original_data <- datalist[[i]]

    # --- Data Cleaning Step ---
    # Keep only rows where ALL required columns have non-NA values
    data <- original_data[complete.cases(original_data[, required_cols]), ]

    if (nrow(data) == 0) {
      warning(paste("Iteration", i, "skipped: No complete cases after removing NAs."))
      next # Skip to the next iteration
    }

    # MISSING EXCLUDED IN ESTIMATION BUT RE-INSERTED IN PREDICTION/RESIDUALS
    options(na.action = "na.exclude")

    # CREATE THE SPLINE
    cb_pollutant <- crossbasis(data[[pollutant_col_name]],
      lag = lag_number,
      argvar = list(type = "lin", cen = FALSE),
      arglag = list(fun = "ns", df = 3)
    )

    # Define the model formula
    model_formula_string <- model_7

    # RUN THE CONDITIONAL-POISSON MODEL
    model <- gnm(model_formula_string,
      family = quasipoisson(),
      eliminate = factor(stratum),
      data = data
    )

    # Q-AIC COMPUTATION
    qaic[i, ] <- fqaic(model)

    # EXTRACT AND SAVE THE RELATED COEF AND VCOV
    pred_stage1_pollutant <- crosspred(cb_pollutant, model, cen = FALSE, cumul = TRUE)
    ymat_pollutant[i, ] <- pred_stage1_pollutant$coef
    list_vcov_pollutant[[i]] <- pred_stage1_pollutant$vcov
  }
)

# RESET WARNING
options(warn = 0)

##########################################
## SECOND STAGE MODEL (POOLED ESTIMATE) ##
##########################################

## LOAD THE PACKAGES (mvmeta PACKAGE IS ASSUMED TO BE INSTALLED)
library(mvmeta)
## Drop Na provinces
valid_rows <- !apply(is.na(ymat_pollutant), 1, all)
ymat_pollutant_dropna <- ymat_pollutant[valid_rows, , drop = FALSE]

list_vcov_pollutant_dropna <- list_vcov_pollutant[valid_rows]
non_null_indices <- !sapply(list_vcov_pollutant_dropna, is.null)
list_vcov_pollutant_dropna <- list_vcov_pollutant_dropna[non_null_indices]
ymat_pollutant_dropna <- ymat_pollutant_dropna[non_null_indices, , drop = FALSE]

# print name of province dropped
drop_provinces <- prov_name[!valid_rows]
if (length(drop_provinces) > 0) {
  cat("Dropped provinces due to NA values:", paste(drop_provinces, collapse = ", "), "\n")
} else {
  cat("No provinces dropped due to NA values.\n")
}

## MULTIVARIATE META-ANALYSIS
mvmeta_results_pollutant <- mvmeta(ymat_pollutant_dropna, list_vcov_pollutant_dropna, method = "reml")
summary(mvmeta_results_pollutant)

## BASIS USED TO PREDICT THE ASSOCIATION, EQUAL TO THAT USED FOR ESTIMATION
pollutant_pred_sequence <- seq(boundpoll[1], boundpoll[2], length = 50)
cb_pred_basis_pollutant <- crossbasis(pollutant_pred_sequence,
  lag = lag_number,
  argvar = list(type = "lin", cen = FALSE),
  arglag = list(fun = "ns", df = 3)
)
summary(cb_pred_basis_pollutant)

## EXPLORE EFFECTS PER 10 PPB INCREASE
cpred_pooled_pollutant <- crosspred(cb_pred_basis_pollutant,
  coef = coef(mvmeta_results_pollutant),
  vcov = vcov(mvmeta_results_pollutant),
  model.link = "log",
  bylag = 1,
  at = 10,
  lag = lag_number,
  cumul = TRUE
)

summary(cpred_pooled_pollutant)

with(cpred_pooled_pollutant, cbind(allRRfit, allRRlow, allRRhigh))

## Plot to RR association-lag curve
png(
  file.path(results_dir, paste0("plot_lag_response", pollutant_col_name, "-", outcome_col_name, ".png")),
  width = 10,
  height = 6,
  units = "in",
  res = 300
)
plot(cpred_pooled_pollutant, "slices",
  var = 10, ylab = "RR and 95% CI", ci.arg = list(density = 15, lwd = 2),
  main = paste(outcome_plot_name, "\nAssociation with a 10-unit increase in", pollutant_col_name)
)
dev.off()
svg(
  file.path(results_dir, paste0("plot_lag_response", pollutant_col_name, "-", outcome_col_name, ".svg")),
  width = 10,
  height = 6,
)
plot(cpred_pooled_pollutant, "slices",
  var = 10, ylab = "RR and 95% CI", ci.arg = list(density = 15, lwd = 2),
  main = paste(outcome_plot_name, "\nAssociation with a 10-unit increase in", pollutant_col_name)
)
dev.off()

# Cumulative curve
png(
  file.path(
    results_dir,
    paste0(
      "plot_cumulative_lag_response",
      pollutant_col_name, "-", outcome_col_name,
      ".png"
    )
  ),
  width = 10,
  height = 6,
  units = "in",
  res = 300
)
plot(cpred_pooled_pollutant, "slices",
  var = 10, cumul = TRUE, ylab = "Cumulative RR and 95% CI",
  main = paste(
    outcome_plot_name,
    "\nCumulative association with a 10-unit increase in",
    pollutant_col_name
  )
)
dev.off()
svg(
  file.path(
    results_dir,
    paste0(
      "plot_cumulative_lag_response",
      pollutant_col_name, "-", outcome_col_name,
      ".svg"
    )
  ),
  width = 10,
  height = 6,
)
plot(cpred_pooled_pollutant, "slices",
  var = 10, cumul = TRUE, ylab = "Cumulative RR and 95% CI",
  main = paste(
    outcome_plot_name,
    "\nCumulative association with a 10-unit increase in",
    pollutant_col_name
  )
)
dev.off()


# create another crosspred obj that have prediction across pollutant ranges
cpred_pooled_pollutant_range <- crosspred(cb_pred_basis_pollutant,
  coef = coef(mvmeta_results_pollutant),
  vcov = vcov(mvmeta_results_pollutant),
  model.link = "log",
  bylag = 1,
  at = pollutant_pred_sequence,
  lag = lag_number,
  cumul = TRUE
)

png(
  file.path(results_dir, paste0("plot_range_association", pollutant_col_name, "-", outcome_col_name, ".png")),
  width = 10,
  height = 6,
  units = "in",
  res = 300
)
plot(cpred_pooled_pollutant_range, "overall",
  lag = 0, ylab = "RR and 95% CI",
  main = paste(outcome_plot_name, "\nAssociation with", pollutant_col_name, "across its range"),
  xlab = paste(pollutant_col_name, "concentration ((μg/m³))"),
)
dev.off()
svg(
  file.path(results_dir, paste0("plot_range_association", pollutant_col_name, "-", outcome_col_name, ".svg")),
  width = 10,
  height = 6,
)
plot(cpred_pooled_pollutant_range, "overall",
  lag = 0, ylab = "RR and 95% CI",
  main = paste(outcome_plot_name, "\nAssociation with", pollutant_col_name, "across its range"),
  xlab = paste(pollutant_col_name, "concentration ((μg/m³))"),
)
dev.off()

## ESTIMATED EFFECTS AT EACH LAG
lagged_RR_pooled_pollutant <- with(cpred_pooled_pollutant, t(rbind(matRRfit, matRRlow, matRRhigh)))
colnames(lagged_RR_pooled_pollutant) <- c("RR", "Lower", "Upper")

## EXPLORE EFFECTS OF THE SELECTED POLLUTANT (CUMULATIVE LAG) FOR EVERY SINGLE PROVINCE (PROVINCE-SPECIFIC)
# Initialize matrix with names of provinces that have valid estimates
valid_provinces <- rownames(ymat_pollutant_dropna)
RR_overall_prov_pollutant <- matrix(NA, length(valid_provinces), 3, dimnames = list(valid_provinces, c("RR", "Lower", "Upper")))
system.time(
  for (i in seq_along(list_vcov_pollutant_dropna)) { # Iterate over the cleaned list
    cat(i, "") # PRINT ITERATION

    # CREATE THE SPLINE
    # Use the ymat_pollutant_dropna and list_vcov_pollutant_dropna which are aligned
    cpred_prov_overall_pollutant <- crosspred(cb_pred_basis_pollutant,
      coef = ymat_pollutant_dropna[i, ],
      vcov = list_vcov_pollutant_dropna[[i]],
      model.link = "log",
      bylag = 1,
      at = 10,
      cumul = TRUE
    )

    RR <- round(with(cpred_prov_overall_pollutant, cbind(allRRfit, allRRlow, allRRhigh)), 4)

    # EXTRACT AND SAVE THE RELATED COEF AND VCOV
    RR_overall_prov_pollutant[i, 1] <- RR[, 1]
    RR_overall_prov_pollutant[i, 2] <- RR[, 2]
    RR_overall_prov_pollutant[i, 3] <- RR[, 3]
  }
)

## EXPLORE EFFECTS OF THE SELECTED POLLUTANT (DISTRIBUTED LAG 0-7) FOR EVERY SINGLE PROVINCE (PROVINCE-SPECIFIC)
# Using valid_provinces for dimnames

# Loop from 0 to the maximum lag number
for (i in 0:lag_number) {
  # --- Create single lag matrices ---

  # 1. Construct the variable name as a string (e.g., "RR_prov_lag0_pollutant")
  var_name_lag <- paste0("RR_prov_lag", i, "_pollutant")

  # 2. Create the matrix that will be assigned to that variable
  matrix_lag <- matrix(NA,
    nrow = length(valid_provinces),
    ncol = 4,
    dimnames = list(valid_provinces, c("Lag", "RR", "Lower", "Upper"))
  )

  # 3. Assign the matrix to the variable name you constructed
  assign(var_name_lag, matrix_lag)


  # --- Create cumulative lag matrices ---

  # 1. Construct the variable name (e.g., "RR_prov_cumulative_lag0_pollutant")
  var_name_cumulative <- paste0("RR_prov_cumulative_lag", i, "_pollutant")

  # 2. Create the matrix
  matrix_cumulative <- matrix(NA,
    nrow = length(valid_provinces),
    ncol = 4,
    dimnames = list(valid_provinces, c("cumulative_Lag", "RR", "Lower", "Upper"))
  )

  # 3. Assign the matrix to the constructed name
  assign(var_name_cumulative, matrix_cumulative)
}

system.time(
  for (i in seq_along(list_vcov_pollutant_dropna)) { # Iterate over the cleaned list
    cat(i, "") # PRINT ITERATION
    prov_name <- names(list_vcov_pollutant_dropna)[i]
    cat(" - ", prov_name, "\n") # PRINT PROVINCE NAME

    # CREATE THE SPLINE
    cpred_prov_lagged_pollutant <- crosspred(cb_pred_basis_pollutant,
      coef = ymat_pollutant_dropna[i, ],
      vcov = list_vcov_pollutant_dropna[[i]],
      model.link = "log",
      bylag = 1,
      at = 10,
      cumul = TRUE
    )

    RR_lagged <- round(with(cpred_prov_lagged_pollutant, t(rbind(matRRfit, matRRlow, matRRhigh))), 4)
    RR_cumul_lagged <- round(with(cpred_prov_lagged_pollutant, t(rbind(cumRRfit, cumRRlow, cumRRhigh))), 4)

    # assign respective value to matrix_lag and matrix_cumulative
    for (lag in 0:lag_number) {
      # Construct the variable names
      var_name_lag <- paste0("RR_prov_lag", lag, "_pollutant")
      var_name_cumulative <- paste0("RR_prov_cumulative_lag", lag, "_pollutant")

      # Get the matrices
      matrix_lag <- get(var_name_lag)
      matrix_cumulative <- get(var_name_cumulative)

      # Assign values to the matrices
      matrix_lag[i, ] <- c(paste("Lag", lag), RR_lagged[lag + 1, ])
      matrix_cumulative[i, ] <- c(paste("cumulative_Lag", lag), RR_cumul_lagged[lag + 1, ])

      # Reassign the updated matrices back to their variable names
      assign(var_name_lag, matrix_lag)
      assign(var_name_cumulative, matrix_cumulative)
    }
  }
)

##################################################
## SAVE RESULTS TO FILES FOR LATER COMPARISON ##
##################################################

# Save Q-AIC values for each province
if (exists("qaic")) {
  qaic_df <- as.data.frame(qaic)
  qaic_df <- cbind(Province = rownames(qaic_df), qaic_df)
  write.csv(qaic_df, file.path(results_dir, "qaic_values.csv"), row.names = FALSE)
  cat("Saved: qaic_values.csv\n")
}

# Save first-stage model coefficients and vcov matrices
# These are ymat_pollutant and list_vcov_pollutant
if (exists("ymat_pollutant") && exists("list_vcov_pollutant")) {
  save(ymat_pollutant, list_vcov_pollutant, file = file.path(results_dir, "first_stage_estimates.RData"))
  cat("Saved: first_stage_estimates.RData\n")
}

# Save the multivariate meta-analysis results object
if (exists("mvmeta_results_pollutant")) {
  saveRDS(mvmeta_results_pollutant, file = file.path(results_dir, "mvmeta_results_pollutant.rds"))
  cat("Saved: mvmeta_results_pollutant.rds\n")

  # Save the summary of the multivariate meta-analysis to a text file
  capture.output(summary(mvmeta_results_pollutant), file = file.path(results_dir, "mvmeta_summary.txt"))
  cat("Saved: mvmeta_summary.txt\n")
}

# Save the pooled overall cumulative Relative Risks from meta-analysis
if (exists("cpred_pooled_pollutant")) {
  pooled_cumulative_RR <- with(cpred_pooled_pollutant, t(rbind(cumRRfit, cumRRlow, cumRRhigh)))
  colnames(pooled_cumulative_RR) <- c("RR", "Lower", "Upper")
  # Convert to data frame
  pooled_cumulative_RR_df <- as.data.frame(pooled_cumulative_RR)

  # Add a lag column
  pooled_cumulative_RR_df$Lag <- 0:(nrow(pooled_cumulative_RR_df) - 1)

  # Reorder columns
  pooled_cumulative_RR_df <- pooled_cumulative_RR_df[, c("Lag", "RR", "Lower", "Upper")]

  # save to csv
  write.csv(pooled_cumulative_RR_df, file.path(results_dir, "pooled_cumulative_RR.csv"), row.names = FALSE)
  cat("Saved: pooled_cumulative_RR.csv\n")

  # 6b. Save the pooled_crosspred object itself for more detailed inspection later
  saveRDS(cpred_pooled_pollutant, file = file.path(results_dir, "cpred_pooled_pollutant.rds"))
  cat("Saved: cpred_pooled_pollutant.rds\n")
}

# 7. Save the pooled lagged Relative Risks from meta-analysis
if (exists("lagged_RR_pooled_pollutant")) {
  # Ensure lagged_RR_pooled_pollutant is a data frame with proper column names
  lagged_RR_pooled_pollutant_df <- as.data.frame(lagged_RR_pooled_pollutant)
  if (ncol(lagged_RR_pooled_pollutant_df) == 3 && !identical(colnames(lagged_RR_pooled_pollutant_df), c("RR", "Lower", "Upper"))) {
    colnames(lagged_RR_pooled_pollutant_df) <- c("RR", "Lower", "Upper") # Ensure correct names if not set by script
  }
  # Add a lag column
  lagged_RR_pooled_pollutant_df$Lag <- 0:(nrow(lagged_RR_pooled_pollutant_df) - 1)
  # Reorder columns
  lagged_RR_pooled_pollutant_df <- lagged_RR_pooled_pollutant_df[, c("Lag", "RR", "Lower", "Upper")]

  write.csv(lagged_RR_pooled_pollutant_df, file.path(results_dir, "pooled_lagged_RR.csv"), row.names = FALSE)
  cat("Saved: pooled_lagged_RR.csv\n")
}

# 8. Save province-specific overall cumulative Relative Risks
if (exists("RR_overall_prov_pollutant")) {
  RR_overall_prov_pollutant_df <- as.data.frame(RR_overall_prov_pollutant)
  RR_overall_prov_pollutant_df <- cbind(Province = rownames(RR_overall_prov_pollutant_df), RR_overall_prov_pollutant_df)
  write.csv(RR_overall_prov_pollutant_df, file.path(results_dir, "province_specific_overall_RR.csv"), row.names = FALSE)
  cat("Saved: province_specific_overall_RR.csv\n")
}

# 9. Save province-specific lagged Relative Risks (Lag 0 to 7)

# First, combine them into a single data frame for easier handling
all_prov_lagged_RR_list <- list()

for (lag_val in 0:lag_number) { # Use lag_number
  obj_name <- paste0("RR_prov_lag", lag_val, "_pollutant")
  if (exists(obj_name)) {
    current_lag_data <- get(obj_name)
    # Convert to data frame, ensure Province is a column
    current_lag_df <- as.data.frame(current_lag_data)
    current_lag_df <- cbind(Province = rownames(current_lag_df), current_lag_df)
    rownames(current_lag_df) <- NULL # Remove row names after putting them in a column

    # Save individual file for this lag
    file_name_individual <- paste0("province_specific_lag", lag_val, "_RR.csv")
    write.csv(current_lag_df, file.path(results_dir, file_name_individual), row.names = FALSE)
    cat("Saved:", file_name_individual, "\n")

    # Add to the list for a combined file (make sure column names are consistent)
    # The 'Lag' column in your current objects is like "Lag 0", "Lag 1", etc.
    # For a combined file, it's better to have a numeric lag column.
    # We'll extract the numeric part or assume based on lag_val.
    # For simplicity, we'll assume the structure is [Province, Lag_Char, RR, Lower, Upper]
    colnames(current_lag_df) <- c("Province", "Lag_Description", "RR", "Lower", "Upper") # Adjust if your columns are named differently
    current_lag_df$Lag_Numeric <- lag_val # Add a numeric lag column
    all_prov_lagged_RR_list[[paste0("lag", lag_val)]] <- current_lag_df
  }
}

for (lag_val in 0:lag_number) { # Use lag_number (which is 7 in your script)
  obj_name <- paste0("RR_prov_cumulative_lag", lag_val, "_pollutant")
  if (exists(obj_name)) {
    current_lag_data <- get(obj_name)
    # Convert to data frame, ensure Province is a column
    current_lag_df <- as.data.frame(current_lag_data)
    current_lag_df <- cbind(Province = rownames(current_lag_df), current_lag_df)
    rownames(current_lag_df) <- NULL # Remove row names after putting them in a column

    # Save individual file for this lag
    file_name_individual <- paste0("province_specific_cumulative_lag", lag_val, "_RR.csv")
    write.csv(current_lag_df, file.path(results_dir, file_name_individual), row.names = FALSE)
    cat("Saved:", file_name_individual, "\n")

    # Add to the list for a combined file
    # The 'Lag' column in your current objects is like "Lag 0", "Lag 1", etc.
    colnames(current_lag_df) <- c("Province", "Lag_Description", "RR", "Lower", "Upper") # Adjust if your columns are named differently
    current_lag_df$Lag_Numeric <- lag_val # Add a numeric lag column
    all_prov_lagged_RR_list[[paste0("cumulative_lag", lag_val)]] <- current_lag_df
  }
}

# Combine all province-specific lagged RRs into one CSV file
if (length(all_prov_lagged_RR_list) > 0) {
  combined_prov_lagged_RR <- do.call(rbind, all_prov_lagged_RR_list)
  # Reorder columns for the combined file
  combined_prov_lagged_RR <- combined_prov_lagged_RR[, c("Province", "Lag_Numeric", "Lag_Description", "RR", "Lower", "Upper")]
  write.csv(combined_prov_lagged_RR, file.path(results_dir, "province_specific_all_lags_RR_combined.csv"), row.names = FALSE)
  cat("Saved: province_specific_all_lags_RR_combined.csv\n")
}


# 10. Save the prediction basis and related parameters
if (exists("cb_pred_basis_pollutant")) {
  saveRDS(cb_pred_basis_pollutant, file = file.path(results_dir, "cb_pred_basis_pollutant.rds"))
  cat("Saved: cb_pred_basis_pollutant.rds\n")
}
if (exists("boundpoll")) {
  saveRDS(boundpoll, file = file.path(results_dir, "boundpoll.rds"))
  cat("Saved: boundpoll.rds\n")
}
# Save key parameters used for this run
params <- list(
  pollutant_col_name = pollutant_col_name,
  lag_number = lag_number,
  date_of_run = Sys.time(),
  ns_temp_df = 3
)
saveRDS(params, file = file.path(results_dir, "run_parameters.rds"))
cat("Saved: run_parameters.rds\n")

cat("\nAll specified results have been saved.\n")

# 11. Plotting Boxplot for lag0, lag1, lag2, lag3, lag0-1, lag 0-2, lag0-3

RR_pooled_lag_filter <- lagged_RR_pooled_pollutant_df %>%
  filter(Lag <= 3) %>%
  mutate(
    Effect_Type = "Individual",
    X_label = paste("lag", Lag, sep = "")
  )

RR_pooled_cum_filter <- pooled_cumulative_RR_df %>%
  filter(between(Lag, 0, 3)) %>%
  mutate(
    Effect_Type = "Cumulative Lag",
    X_label = paste("lag0-", Lag, sep = "")
  )

RR_pooled_plot_df <- rbind(RR_pooled_lag_filter, RR_pooled_cum_filter)
print(RR_pooled_plot_df)

# Restructure the RR_pooled_plot_df for reporting
RR_pooled_df <- RR_pooled_plot_df %>%
  remove_rownames() %>%
  column_to_rownames(var = "X_label") %>%
  mutate_if(is.numeric, round, 5) %>%
  mutate("95% CI" = paste(Lower, "-", Upper, sep = " ")) %>%
  select(RR, "95% CI")
write.csv(RR_pooled_df, file.path(results_dir, "RR_pooled_df.csv"), row.names = TRUE)

xtable(RR_pooled_df, digits = 2, align = "l|cc") %>%
  print(file = file.path(results_dir, "RR_pooled_df.tex"), include.rownames = TRUE)

# Select only relevant cols
label_order <- c("lag0", "lag1", "lag2", "lag3", "lag0-1", "lag0-2", "lag0-3")
RR_pooled_plot_df <- RR_pooled_plot_df %>%
  filter(X_label %in% label_order) %>%
  select(RR, Lower, Upper, X_label)
print(RR_pooled_plot_df)

ggplot(RR_pooled_lag_filter, aes(x = X_label, y = RR)) +
  geom_point(size = 3, aes(color = "RR")) + # Plots the central RR points
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = "95% Confidence Interval"), width = 0.2) + # Adds the upper and lower bounds as error bars
  geom_line(aes(y = RR, color = "RR"), linetype = "dashed") +
  geom_hline(yintercept = 1, color = "red", linetype = "dotted", linewidth = 1, alpha = 0.4) +
  annotate(geom = "text", x = 1, y = 0.9995, label = "RR = 1", color = "black", size = 4) +
  scale_color_manual(values = c("RR" = "blue", "95% Confidence Interval" = "black")) +
  labs(
    x = "Lag",
    y = "Relative Risk (RR)",
    color = "Legend"
  ) +
  labs(
    subtitle = outcome_plot_name,
    title = paste("Relative Risk by Lag Day for", pollutant_col_name),
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )
ggsave(
  file.path(
    results_dir,
    paste0("RR_lag_day", pollutant_col_name, ".png")
  ),
  width = 8,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(
    results_dir,
    paste0("RR_lag_day", pollutant_col_name, ".svg")
  ),
  width = 8,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

# Plot the RR with lag and lag_cumu in X
ggplot(RR_pooled_cum_filter, aes(x = X_label, y = RR)) +
  geom_point(size = 3, aes(color = "RR")) + # Plots the central RR points
  geom_errorbar(aes(ymin = Lower, ymax = Upper, color = "95% Confidence Interval"), width = 0.2) + # Adds the upper and lower bounds as error bars
  geom_line(aes(y = RR, color = "RR"), linetype = "dashed") +
  geom_hline(yintercept = 1, color = "red", linetype = "dotted", linewidth = 1, alpha = 0.4) +
  annotate(geom = "text", x = 1, y = 0.9995, label = "RR = 1", color = "black", size = 4) +
  scale_color_manual(values = c("RR" = "blue", "95% Confidence Interval" = "black")) +
  labs(
    title = paste("Cumulative Relative Risk by Lag for", pollutant_col_name),
    subtitle = outcome_plot_name,
    x = "Lag",
    y = "Cumulative Relative Risk (RR)",
    color = "Legend"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )
ggsave(
  file.path(
    results_dir,
    paste0("RR_cumu_lag", pollutant_col_name, ".png")
  ),
  width = 8,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(
    results_dir,
    paste0("RR_cumu_lag", pollutant_col_name, ".svg")
  ),
  width = 8,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

rr_cum_line_plot <- RR_pooled_cum_filter %>%
  mutate(across(
    c("RR", "Lower", "Upper"),
    as.numeric
  )) %>%
  mutate(
    row_number = row_number()
  ) %>%
  ggplot(aes(x = row_number, y = RR, ymin = Lower, ymax = Upper)) +
  geom_ribbon(fill = "grey80", alpha = 0.5) +
  geom_line(color = "blue", linewidth = 0.8) +
  scale_x_continuous(
    breaks = seq(1, nrow(RR_pooled_cum_filter), by = 1),
    labels = c("Lag0", "Lag0-1", "Lag0-2", "Lag0-3")
  ) +
  labs(
    subtitle = outcome_plot_name,
    title = "Cumulative association with a 10-unit increase in PM2.5",
    x = "Lag",
    y = "Cumulative RR and 95% CI"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  )
ggsave(
  file.path(
    results_dir,
    paste0("RR_cumu_lag_line_", pollutant_col_name, ".png") # Added _new to avoid overwriting
  ),
  width = 8,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(
    results_dir,
    paste0("RR_cumu_lag_line_", pollutant_col_name, ".svg") # Added _new to avoid overwriting
  ),
  width = 8,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
# Plot province specific relative risk of incrase in 10 PM2.5 values as a forest plot

forest_plot_df <- as.data.frame(RR_prov_lag0_pollutant) %>%
  rownames_to_column(var = "Province") %>%
  mutate(
    pooled = FALSE
  ) %>%
  select(-Lag)

new_pooled_row <- lagged_RR_pooled_pollutant_df %>%
  filter(Lag == 0) %>%
  mutate(
    RR = round(RR, 4),
    Lower = round(Lower, 4),
    Upper = round(Upper, 4),
    Province = "Pooled",
    pooled = TRUE
  ) %>%
  select(Province, RR, Lower, Upper, pooled)

forest_plot_df <- rbind(forest_plot_df, new_pooled_row)
rownames(forest_plot_df) <- NULL

forest_plot_province <- forest_plot_df %>%
  mutate(sort_order = ifelse(pooled, 1, 2)) %>%
  mutate(across(c(RR, Upper, Lower), as.numeric)) %>%
  arrange(sort_order, Province) %>%
  mutate(Province = fct_rev(fct_inorder(Province))) %>%
  ggplot(aes(
    x = RR,
    y = Province,
    color = pooled
  )) +
  geom_errorbarh(
    aes(
      xmin = Lower,
      xmax = Upper,
    )
  ) +
  geom_point(data = . %>% filter(pooled == FALSE), size = 2, shape = 16) +
  geom_point(data = . %>% filter(pooled == TRUE), size = 3, shape = 18) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5, alpha = 0.6) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
  labs(
    subtitle = outcome_plot_name,
    title = paste(
      "Relative Risk per 10 μg/m³ increase of",
      pollutant_col_name,
      sep = " "
    ),
    x = "Relative Risk (RR) with 95% Confidence Interval",
    y = "Province"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "none"
  )
ggsave(
  file.path(results_dir, paste0("forest_plot_province_", pollutant_col_name, ".png")),
  plot = forest_plot_province,
  width = 8,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(results_dir, paste0("forest_plot_province_", pollutant_col_name, ".svg")),
  plot = forest_plot_province,
  width = 8,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)
################################################
##                  END OF SAVING               ##
##################################################
