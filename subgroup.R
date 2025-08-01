# ===================================================================
# Subgroup Analysis for DLNM Models
#
# This script performs a two-stage distributed lag non-linear model
# (DLNM) analysis on specified subgroups of a dataset.
#
# Structure:
# 1. A main function `run_analysis_for_subgroup` encapsulates the
#    entire analysis pipeline.
# 2. The script loads and prepares the full dataset.
# 3. It iterates through unique values of a specified subgrouping
#    column (e.g., 'region_code').
# 4. In each iteration, it calls the main function to run the analysis
#    for that subgroup and saves results in a separate folder.
# ===================================================================

# ============================
# Load Necessary Packages
# ============================
library(dplyr)
library(dlnm)
library(splines)
library(tsModel)
library(gnm)
library(tidyverse)
library(mvmeta)
library(xtable)

# Set seed for reproducibility
set.seed(123)


# ===================================================================
# MAIN ANALYSIS FUNCTION
# ===================================================================
#' @title Run Two-Stage DLNM Analysis for a Data Subgroup
#' @description This function performs a first-stage province-specific
#'              time-series analysis followed by a second-stage meta-analysis
#'              to pool the estimates for a given data subgroup.
#' @param subgroup_data A data frame containing the data for one subgroup.
#' @param subgroup_name A character string used to name the output directory.
#' @param pollutant_col_name The name of the pollutant variable column.
#' @param outcome_col_name The name of the outcome variable column.
#' @param lag_number The maximum lag in days.
#' @param base_output_dir The root directory where the subgroup folder will be created.
#'
run_analysis_for_subgroup <- function(subgroup_data,
                                      subgroup_name,
                                      pollutant_col_name,
                                      outcome_col_name,
                                      lag_number,
                                      base_output_dir) {
  cat("\n========================================================\n")
  cat("RUNNING ANALYSIS FOR SUBGROUP:", subgroup_name, "\n")
  cat("========================================================\n\n")

  # Use the subgroup data for the analysis
  dat <- subgroup_data

  # --- 1. Setup and Preparation ---

  # Create a dedicated directory for the subgroup's results
  results_dir <- file.path(base_output_dir, subgroup_name)
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  cat("Saving results to directory:", file.path(getwd(), results_dir), "\n")

  # Get province names within the subgroup
  prov_name <- as.character(unique(dat$prov_name))

  # Separate dataframe by provinces
  datalist <- lapply(
    prov_name, # Use prov_name which is specific to the subgroup
    function(provname) dat[dat$prov_name == provname, ]
  )
  names(datalist) <- prov_name

  m <- length(datalist)

  # --- 2. First Stage: Province-Specific Models ---

  # Define the model formula.
  # This model was pre-selected based on QAIC in the original script.
  # Model: cb_pollutant + ns(temperature, 3) + as.factor(is_holiday)
  model_formula <- as.formula(paste(outcome_col_name, "~ cb_pollutant + ns(temperature, 3) + as.factor(is_holiday)"))
  required_cols <- c("temperature", "is_holiday", pollutant_col_name, outcome_col_name)

  # Initialize objects to store results
  ymat_pollutant <- matrix(NA, m, 3, dimnames = list(prov_name, paste0("spl", seq(3))))
  list_vcov_pollutant <- vector("list", m)
  names(list_vcov_pollutant) <- prov_name

  # Function to compute Q-AIC
  fqaic <- function(model) {
    loglik <- sum(dpois(model$y, model$fitted.values, log = TRUE))
    phi <- summary(model)$dispersion
    qaic <- -2 * loglik + 2 * summary(model)$df[3] * phi
    return(qaic)
  }

  qaic_matrix <- matrix(NA, m, 1, dimnames = list(prov_name, "Q-AIC"))

  options(warn = -1) # Suppress warnings for prediction beyond boundaries

  # Loop through each province in the subgroup
  for (i in seq_along(datalist)) {
    cat("Processing province", i, "of", m, ":", names(datalist)[i], "\n")

    # Get data for the province and remove rows with NAs in required columns
    data_prov <- datalist[[i]]
    data_prov_clean <- data_prov[complete.cases(data_prov[, required_cols]), ]

    if (nrow(data_prov_clean) == 0) {
      warning(paste("Province", names(datalist)[i], "skipped: No complete cases."))
      next
    }

    # Create the cross-basis matrix
    cb_pollutant <- crossbasis(data_prov_clean[[pollutant_col_name]],
      lag = lag_number,
      argvar = list(type = "lin", cen = FALSE),
      arglag = list(fun = "ns", df = 3)
    )

    # Run the conditional quasi-Poisson model
    model <- gnm(model_formula,
      family = quasipoisson(),
      eliminate = factor(stratum),
      data = data_prov_clean,
      na.action = "na.exclude"
    ) # Use cleaned data

    # Store results
    qaic_matrix[i, ] <- fqaic(model)
    pred_stage1 <- crosspred(
      cb_pollutant,
      model,
      cen = FALSE,
      model.link = "log",
      by = 1,
      at = 10,
      cumul = TRUE
    )
    ymat_pollutant[i, ] <- pred_stage1$coef
    list_vcov_pollutant[[i]] <- pred_stage1$vcov
  }

  options(warn = 0) # Reset warnings

  # --- 3. Second Stage: Multivariate Meta-Analysis ---

  # Remove provinces with NA results from the first stage
  valid_rows <- !apply(is.na(ymat_pollutant), 1, all)
  ymat_pollutant_dropna <- ymat_pollutant[valid_rows, , drop = FALSE]
  list_vcov_pollutant_dropna <- list_vcov_pollutant[valid_rows]

  if (nrow(ymat_pollutant_dropna) == 1) {
    # in case of insufficient data, save the first-stage esimates in the same pattern

    cat("\nSkipping meta-analysis for subgroup", subgroup_name, "due to < 2 valid provinces.\n")

    # a. Plot and save lag-response curve
    png(file.path(results_dir, "plot_lag_response.png"))
    plot(pred_stage1, "slices", var = 10, ylab = "RR", main = paste(subgroup_name, "- Single Lag-Response"))
    dev.off()

    # b. Plot and save cumulative lag-response curve
    png(file.path(results_dir, "plot_cumulative_lag_response.png"))
    plot(pred_stage1, "slices", var = 10, cumul = TRUE, ylab = "Cumulative RR", main = paste(subgroup_name, "- Single Cumulative Lag-Response"))
    dev.off()

    # c. Save single RR data (lagged and cumulative)

    RR_pooled_lagged <- with(pred_stage1, t(rbind(matRRfit, matRRlow, matRRhigh)))
    colnames(RR_pooled_lagged) <- c("RR", "Lower", "Upper")

    RR_pooled_cumulative <- with(pred_stage1, t(rbind(cumRRfit, cumRRlow, cumRRhigh)))
    colnames(RR_pooled_cumulative) <- c("CumRR", "Lower", "Upper")

    RR_pooled_lagged <- as.data.frame(RR_pooled_lagged)
    RR_pooled_cumulative <- as.data.frame(RR_pooled_cumulative)

    RR_pooled_lagged$Lag <- 0:lag_number
    RR_pooled_cumulative$Lag <- 0:lag_number

    # properly sort the columns order
    RR_pooled_lagged <- RR_pooled_lagged[, c("Lag", "RR", "Lower", "Upper")]
    RR_pooled_cumulative <- RR_pooled_cumulative[, c("Lag", "CumRR", "Lower", "Upper")]

    write.csv(RR_pooled_lagged, file.path(results_dir, "pooled_lagged_RR.csv"), row.names = FALSE)
    write.csv(RR_pooled_cumulative, file.path(results_dir, "pooled_cumulative_RR.csv"), row.names = FALSE)

    # d. Save first-stage estimates
    save(ymat_pollutant, list_vcov_pollutant, file = file.path(results_dir, "first_stage_estimates.RData"))

    cat("\nAnalysis for subgroup", subgroup_name, "is complete.\n")

    return() # Exit the function if not enough data to pool
  }

  # Run meta-analysis
  mvmeta_results <- mvmeta(ymat_pollutant_dropna, list_vcov_pollutant_dropna, method = "reml")

  # --- 4. Generate and Save Pooled Results & Plots ---

  # Define prediction range and create prediction basis
  poll_range <- range(dat[[pollutant_col_name]], na.rm = TRUE)
  pred_sequence <- seq(poll_range[1], poll_range[2], length = 50)
  cb_pred_basis <- crossbasis(pred_sequence,
    lag = lag_number,
    argvar = list(type = "lin", cen = FALSE),
    arglag = list(fun = "ns", df = 3)
  )

  # Predict pooled exposure-lag-response for a 10-unit increase
  cpred_pooled <- crosspred(cb_pred_basis,
    coef = coef(mvmeta_results),
    vcov = vcov(mvmeta_results),
    model.link = "log",
    by = 1, # Ensure bylag is set for compatibility, 1 is a safe default
    at = 10,
    cumul = TRUE
  )

  # a. Save meta-analysis summary
  capture.output(summary(mvmeta_results), file = file.path(results_dir, "mvmeta_summary.txt"))

  # b. Plot and save lag-response curve
  png(file.path(results_dir, "plot_lag_response.png"))
  plot(cpred_pooled, "slices", var = 10, ylab = "RR", main = paste(subgroup_name, "- Pooled Lag-Response"))
  dev.off()

  # c. Plot and save cumulative lag-response curve
  png(file.path(results_dir, "plot_cumulative_lag_response.png"))
  plot(cpred_pooled, "slices", var = 10, cumul = TRUE, ylab = "Cumulative RR", main = paste(subgroup_name, "- Pooled Cumulative Lag-Response"))
  dev.off()

  # d. Save pooled RR data (lagged and cumulative)
  RR_pooled_lagged <- with(cpred_pooled, t(rbind(matRRfit, matRRlow, matRRhigh)))
  colnames(RR_pooled_lagged) <- c("RR", "Lower", "Upper")

  RR_pooled_cumulative <- with(cpred_pooled, t(rbind(cumRRfit, cumRRlow, cumRRhigh)))
  colnames(RR_pooled_cumulative) <- c("CumRR", "Lower", "Upper")

  RR_pooled_lagged <- as.data.frame(RR_pooled_lagged)
  RR_pooled_cumulative <- as.data.frame(RR_pooled_cumulative)

  RR_pooled_lagged$Lag <- 0:lag_number
  RR_pooled_cumulative$Lag <- 0:lag_number

  # properly sort the columns order
  RR_pooled_lagged <- RR_pooled_lagged[, c("Lag", "RR", "Lower", "Upper")]
  RR_pooled_cumulative <- RR_pooled_cumulative[, c("Lag", "CumRR", "Lower", "Upper")]

  write.csv(RR_pooled_lagged, file.path(results_dir, "pooled_lagged_RR.csv"), row.names = FALSE)
  write.csv(RR_pooled_cumulative, file.path(results_dir, "pooled_cumulative_RR.csv"), row.names = FALSE)

  # e. Save first-stage estimates
  save(ymat_pollutant, list_vcov_pollutant, file = file.path(results_dir, "first_stage_estimates.RData"))

  cat("\nAnalysis for subgroup", subgroup_name, "is complete.\n")
}


# ===================================================================
# SCRIPT EXECUTION
# ===================================================================

# --- 1. Load and Prepare Full Dataset ---
cat("Loading and preparing the main dataset...\n")

# Make sure the working directory is set to your project root
# setwd("~/Projects/hf_pm_analysis")
dat_full <- read.csv("data/daily_all_weather_30jul.csv")

# Properly format date and create time-stratification variables
dat_full$date_start <- as.Date(dat_full$date_start, "%Y-%m-%d")
dat_full$dow <- as.factor(weekdays(dat_full$date_start))
dat_full$month <- as.factor(months(dat_full$date_start))
dat_full$year <- as.factor(format(dat_full$date_start, format = "%Y"))
dat_full$stratum <- as.factor(dat_full$year:dat_full$month:dat_full$dow)

# --- 2. Define Analysis Parameters ---

## !! IMPORTANT: SPECIFY YOUR PARAMETERS HERE !! ##
POLLUTANT_VAR <- "PM2.5"
OUTCOME_VAR <- "hf_prim"
LAG <- 7
BASE_OUTPUT_DIR <- "output/regional_subgroup_analysis_30jul"
# The column name in 'dat_full' used to split the data into subgroups.
# Ensure this column exists in your dataframe.
SUBGROUPING_COLUMN <- "nhso_region_code"


# --- 3. Run Analysis Loop for Each Subgroup ---

# Check if the subgrouping column exists
if (!SUBGROUPING_COLUMN %in% names(dat_full)) {
  stop(paste("Error: The specified subgrouping column '", SUBGROUPING_COLUMN, "' does not exist in the data frame.", sep = ""))
}

# Get unique subgroup identifiers
subgroups <- unique(dat_full[[SUBGROUPING_COLUMN]])
subgroups <- subgroups[!is.na(subgroups)] # Remove NA from the list of subgroups

cat(paste("\nFound", length(subgroups), "subgroups to analyze from column '", SUBGROUPING_COLUMN, "'.\n", sep = " "))

# Loop through each subgroup, filter data, and run the analysis function
for (group in subgroups) {
  # Create a clean name for the directory (e.g., "region_1")
  subgroup_folder_name <- paste0(SUBGROUPING_COLUMN, "_", group)

  # Filter the main dataframe to get data for the current subgroup
  data_for_group <- dat_full %>%
    filter(.data[[SUBGROUPING_COLUMN]] == group)

  # Call the main analysis function
  run_analysis_for_subgroup(
    subgroup_data = data_for_group,
    subgroup_name = subgroup_folder_name,
    pollutant_col_name = POLLUTANT_VAR,
    outcome_col_name = OUTCOME_VAR,
    lag_number = LAG,
    base_output_dir = BASE_OUTPUT_DIR
  )
}

cat("\n\nAll subgroup analyses are complete.\n")

cat("\nPrinting NHSO region codes mapping.")

nhso_code_map <- dat %>%
  select(nhso_region_code, prov_name) %>%
  unique() %>%
  group_by(nhso_region_code) %>%
  mutate(prov_name = paste(prov_name, collapse = ", ")) %>%
  unique() %>%
  arrange(nhso_region_code)
write.csv(nhso_code_map, file.path(BASE_OUTPUT_DIR, "nhso_region_codes_mapping.csv"), row.names = FALSE)

cat("=====Finished=====")
