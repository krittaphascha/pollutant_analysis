setwd("D:\\Tiger\\afhf_pm")
set.seed(123)

# read data
data_path <- "data\\output\\daily_all_weather_30May.csv"
dat <- read.csv(data_path)

options(na.action = "na.exclude")


# properly format date
dat$date_start <- as.Date(dat$date_start, "%Y-%m-%d")
dat$dow <- as.factor(weekdays(dat$date_start))
dat$month <- as.factor(months(dat$date_start))
dat$year <- as.factor(format(dat$date_start, format = "%Y"))

# get strata for case crossover analysis
dat$stratum <- as.factor(dat$year:dat$month:dat$dow)

provinces <- as.character(unique(dat$prov_name))

# separate dataframe by provinces
datalist <- lapply(provinces, function(provname) dat[dat$prov_name == provname, ])
names(datalist) <- provinces

library(dplyr); library(dlnm); library(splines); library(tsModel); library(gnm)

## SELECT POLLUTANT OF INTEREST
pollutant_col_name <- "PM2.5"
outcome_col_name <- "hf_prim_hos_count"

# Create a directory to store the results
# This directory will be named based on the pollutant and lag number for easy identification
results_dir <- paste0("R//output//results_", pollutant_col_name, '-' ,outcome_col_name)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE) # recursive = TRUE if pollutant_col_name could contain slashes (e.g. "PM2.5/O3")

cat("\nSaving results to directory:", file.path(getwd(), results_dir), "\n")

## RANGE FOR METEOROLOGICAL VARIABLES ##
rangepoll <- t(sapply(datalist, function(x) range(x[[pollutant_col_name]],  na.rm = T)))

# ## DLNM PARAMETERS 
boundpoll <- c(min(rangepoll), max(rangepoll))

## ADDITIONAL INFORMATION
m <- length(datalist)

lag_number <- 7

####################################################
## FIRST STAGE MODEL (PROVINCE-SPECIFIC ESTIMATE) ##
####################################################

## BUILT OBJECTS WHERE RESULTS WILL BE STORED: # ymat IS THE MATRIX FOR THE OUTCOME PARAMETERS
## Slist IS THE LIST WITH (CO)VARIANCE MATRICES # LINEAR FOR NO2 AND 3 DFs FOR ITS LAG
ymat_pollutant <- matrix(NA, length(datalist), 3, dimnames = list(provinces, paste("spl", seq(3), sep = ""))) # 3DF
list_vcov_pollutant <- vector("list", length(datalist))
names(list_vcov_pollutant) <- provinces

## FUNCTION TO COMPUTE THE Q-AIC IN QUASI-POISSON MODELS
fqaic <- function(model) {
  loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
  phi <- summary(model)$dispersion
  qaic <- -2*loglik + 2*summary(model)$df[3]*phi
  return(qaic)
}

## MATRIX FOR Q-AIC VALUE
qaic <- matrix(NA, length(datalist), 1, dimnames = list(provinces, paste("Q-AIC")))

## WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED: RUN THE FIRST STAGE ANALYSIS
options(warn = -1)

## add COLUMNS that cannot be nulled
required_cols <- c("temperature", "humidity", "pressure")

## Stepwise regression for model selection first
## Testing on BKK data
# library(MuMIn)
# 
# data_test <- datalist[[1]]
# cb_pollutant_test <- crossbasis(data_test[[pollutant_col_name]], lag = 7, argvar = list(type = "lin", cen = FALSE), arglag = list(fun = "ns", df = 3))
# 
# options(na.action='na.exclude')
# full.model <- gnm(hf_prim_hos_count ~ cb_pollutant_test+ ns(temperature, 3) + ns(humidity, 3) + ns(pressure, 3) + as.factor(is_holiday) + as.factor(is_lotteryday),
#                  family = quasipoisson(),
#                  eliminate = factor(stratum),
#                  data = data_test)
# 
# options(na.action='na.pass')
# dd <- dredge(full.model, rank=fqaic)
# print(dd)

# options(na.action='na.exclude')

system.time(
  for(i in seq(m)) {
    
    # PRINT ITERATION
    cat(i,"")
    
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
    options(na.action="na.exclude")
    
    # CREATE THE SPLINE
    cb_pollutant <- crossbasis(data[[pollutant_col_name]], lag = lag_number, argvar = list(type = "lin", cen = FALSE), arglag = list(fun = "ns", df = 3))
    
    # Define the model formula
    # 1. Define the right-hand side (predictors) of the formula
    # model_predictors <- "cb_pollutant + ns(temperature, 3) + ns(humidity, 3) + ns(pressure, 3) + as.factor(is_holiday)"
    model_predictors <- "cb_pollutant + ns(temperature, 3) + as.factor(is_holiday)"
    
    # 2. Paste the outcome and predictors together to form the full formula string
    model_formula_string <- paste(outcome_col_name, "~", model_predictors) 
    model_formula_string <- as.formula(model_formula_string)
    
    # RUN THE CONDITIONAL-POISSON MODEL
    model <- gnm(model_formula_string,
                 family = quasipoisson(),
                 eliminate = factor(stratum),
                 data = data)

    # Q-AIC COMPUTATION
    qaic[i,] <- fqaic(model)
    
    # EXTRACT AND SAVE THE RELATED COEF AND VCOV
    pred_stage1_pollutant <- crosspred(cb_pollutant, model, cen = FALSE, cumul=TRUE)
    ymat_pollutant[i,] <- pred_stage1_pollutant$coef
    list_vcov_pollutant[[i]] <- pred_stage1_pollutant$vcov
    
  })

mean(qaic, na.rm=TRUE)

# RESET WARNING
options(warn = 0)

##########################################
## SECOND STAGE MODEL (POOLED ESTIMATE) ##
##########################################

## LOAD THE PACKAGES (mvmeta PACKAGE IS ASSUMED TO BE INSTALLED)
library(mvmeta);

## Drop Na provinces
valid_rows <- !apply(is.na(ymat_pollutant), 1, all) 
ymat_pollutant_dropna <- ymat_pollutant[valid_rows, , drop = FALSE]

list_vcov_pollutant_dropna <- list_vcov_pollutant[valid_rows]
non_null_indices <- !sapply(list_vcov_pollutant_dropna, is.null)
list_vcov_pollutant_dropna <- list_vcov_pollutant_dropna[non_null_indices]
ymat_pollutant_dropna <- ymat_pollutant_dropna[non_null_indices, , drop = FALSE]

# Check if there are enough studies for meta-analysis
if (nrow(ymat_pollutant_dropna) < 2) {
  stop("Insufficient data for meta-analysis after removing NAs/NULLs. Need at least 2 provinces with valid estimates.")
}

## MULTIVARIATE META-ANALYSIS
mvmeta_results_pollutant <- mvmeta(ymat_pollutant_dropna, list_vcov_pollutant_dropna, method = "reml")
summary(mvmeta_results_pollutant)

## BASIS USED TO PREDICT THE ASSOCIATION, EQUAL TO THAT USED FOR ESTIMATION
pollutant_pred_sequence <- seq(boundpoll[1], boundpoll[2], length = 30)
cb_pred_basis_pollutant <- crossbasis(pollutant_pred_sequence, 
                                      lag = lag_number, 
                                      argvar = list(type = "lin", cen = FALSE), 
                                      arglag = list(fun = "ns", df = 3))

## EXPLORE EFFECTS PER 10 PPB INCREASE
cpred_pooled_pollutant <- crosspred(cb_pred_basis_pollutant, 
                                    coef = coef(mvmeta_results_pollutant), 
                                    vcov = vcov(mvmeta_results_pollutant), 
                                    model.link = "log", 
                                    bylag = 1, 
                                    at = 10, 
                                    lag = lag_number,
                                    cumul = TRUE)

## OVEARLL CUMULATIVE EFFECT
with(cpred_pooled_pollutant, cbind(allRRfit, allRRlow, allRRhigh))

## Plot to RR association-lag curve
png(file.path(results_dir, paste0('plot_lag_response', pollutant_col_name, '-',outcome_col_name, '.png')))
plot(cpred_pooled_pollutant, "slices", var=10, col=3, ylab="Overlall RR", ci.arg=list(density=15,lwd=2),
     main=paste(outcome_col_name, "\nAssociation with a 10-unit increase in", pollutant_col_name))
dev.off()

# Cumulative curve
png(file.path(results_dir, paste0('plot_cumulative_lag_response', pollutant_col_name, '-',outcome_col_name, '.png')))
plot(cpred_pooled_pollutant, "slices", var=10, col=2, cumul=TRUE, ylab="Overall Cumulative RR",
     main=paste(outcome_col_name, "\nCumulative association with a 10-unit increase in", pollutant_col_name))
dev.off()

## ESTIMATED EFFECTS AT EACH LAG
lagged_RR_pooled_pollutant <- with(cpred_pooled_pollutant, t(rbind(matRRfit, matRRlow, matRRhigh)))
colnames(lagged_RR_pooled_pollutant) <- c("RR", "Lower", "Upper")

## EXPLORE EFFECTS OF THE SELECTED POLLUTANT (CUMULATIVE LAG) FOR EVERY SINGLE PROVINCE (PROVINCE-SPECIFIC)
# Initialize matrix with names of provinces that have valid estimates
valid_provinces <- rownames(ymat_pollutant_dropna)
RR_overall_prov_pollutant <- matrix(NA, length(valid_provinces), 3, dimnames = list(valid_provinces, c("RR", "Lower", "Upper")))
system.time(
  for(i in seq_along(list_vcov_pollutant_dropna)) { # Iterate over the cleaned list
    cat(i, "") # PRINT ITERATION
    
    # CREATE THE SPLINE
    # Use the ymat_pollutant_dropna and list_vcov_pollutant_dropna which are aligned
    cpred_prov_overall_pollutant <- crosspred(cb_pred_basis_pollutant, 
                                              coef = ymat_pollutant_dropna[i,], 
                                              vcov = list_vcov_pollutant_dropna[[i]], 
                                              model.link = "log", 
                                              bylag = 1, 
                                              at = 10,
                                              cumul = TRUE)
    
    RR <- round(with(cpred_prov_overall_pollutant, cbind(allRRfit, allRRlow, allRRhigh)), 4)
    
    # EXTRACT AND SAVE THE RELATED COEF AND VCOV
    RR_overall_prov_pollutant[i,1] <- RR[,1]
    RR_overall_prov_pollutant[i,2] <- RR[,2]
    RR_overall_prov_pollutant[i,3] <- RR[,3]
  })

## EXPLORE EFFECTS OF THE SELECTED POLLUTANT (DISTRIBUTED LAG 0-7) FOR EVERY SINGLE PROVINCE (PROVINCE-SPECIFIC)
# Using valid_provinces for dimnames
RR_prov_lag0_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("Lag","RR","Lower","Upper")))
RR_prov_lag1_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("Lag","RR","Lower","Upper")))
RR_prov_lag2_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("Lag","RR","Lower","Upper")))
RR_prov_lag3_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("Lag","RR","Lower","Upper")))
RR_prov_lag4_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("Lag","RR","Lower","Upper")))
RR_prov_lag5_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("Lag","RR","Lower","Upper")))
RR_prov_lag6_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("Lag","RR","Lower","Upper")))
RR_prov_lag7_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("Lag","RR","Lower","Upper")))

RR_prov_cumulative_lag0_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("cumulative_Lag","RR","Lower","Upper")))
RR_prov_cumulative_lag1_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("cumulative_Lag","RR","Lower","Upper")))
RR_prov_cumulative_lag2_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("cumulative_Lag","RR","Lower","Upper")))
RR_prov_cumulative_lag3_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("cumulative_Lag","RR","Lower","Upper")))
RR_prov_cumulative_lag4_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("cumulative_Lag","RR","Lower","Upper")))
RR_prov_cumulative_lag5_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("cumulative_Lag","RR","Lower","Upper")))
RR_prov_cumulative_lag6_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("cumulative_Lag","RR","Lower","Upper")))
RR_prov_cumulative_lag7_pollutant <- matrix(NA, length(valid_provinces), 4, dimnames = list(valid_provinces, c("cumulative_Lag","RR","Lower","Upper")))

system.time(
  for(i in seq_along(list_vcov_pollutant_dropna)) { # Iterate over the cleaned list
    cat(i,"") # PRINT ITERATION
    
    # CREATE THE SPLINE
    cpred_prov_lagged_pollutant <- crosspred(cb_pred_basis_pollutant, 
                                             coef = ymat_pollutant_dropna[i,], 
                                             vcov = list_vcov_pollutant_dropna[[i]], 
                                             model.link = "log", 
                                             bylag = 1, 
                                             at = 10,
                                             cumul = TRUE)
    
    RR_lagged <- round(with(cpred_prov_lagged_pollutant, t(rbind(matRRfit, matRRlow, matRRhigh))), 4)
    RR_cumul_lagged <- round(with(cpred_prov_lagged_pollutant, t(rbind(cumRRfit, cumRRlow, cumRRhigh))), 4)
    
    # EXTRACT AND SAVE THE RELATED COEF AND VCOV
    # Ensure RR_lagged has 8 rows (for lag 0 to 7)
    if(nrow(RR_lagged) == 8) {
      RR_prov_lag0_pollutant[i,] <- c("Lag 0", RR_lagged[1,1], RR_lagged[1,2], RR_lagged[1,3])
      RR_prov_lag1_pollutant[i,] <- c("Lag 1", RR_lagged[2,1], RR_lagged[2,2], RR_lagged[2,3])
      RR_prov_lag2_pollutant[i,] <- c("Lag 2", RR_lagged[3,1], RR_lagged[3,2], RR_lagged[3,3])
      RR_prov_lag3_pollutant[i,] <- c("Lag 3", RR_lagged[4,1], RR_lagged[4,2], RR_lagged[4,3])
      RR_prov_lag4_pollutant[i,] <- c("Lag 4", RR_lagged[5,1], RR_lagged[5,2], RR_lagged[5,3])
      RR_prov_lag5_pollutant[i,] <- c("Lag 5", RR_lagged[6,1], RR_lagged[6,2], RR_lagged[6,3])
      RR_prov_lag6_pollutant[i,] <- c("Lag 6", RR_lagged[7,1], RR_lagged[7,2], RR_lagged[7,3])
      RR_prov_lag7_pollutant[i,] <- c("Lag 7", RR_lagged[8,1], RR_lagged[8,2], RR_lagged[8,3])
      RR_prov_cumulative_lag0_pollutant[i,] <- c("cumulative_Lag 0", RR_cumul_lagged[1,1], RR_cumul_lagged[1,2], RR_cumul_lagged[1,3])
      RR_prov_cumulative_lag1_pollutant[i,] <- c("cumulative_Lag 1", RR_cumul_lagged[2,1], RR_cumul_lagged[2,2], RR_cumul_lagged[2,3])
      RR_prov_cumulative_lag2_pollutant[i,] <- c("cumulative_Lag 2", RR_cumul_lagged[3,1], RR_cumul_lagged[3,2], RR_cumul_lagged[3,3])
      RR_prov_cumulative_lag3_pollutant[i,] <- c("cumulative_Lag 3", RR_cumul_lagged[4,1], RR_cumul_lagged[4,2], RR_cumul_lagged[4,3])
      RR_prov_cumulative_lag4_pollutant[i,] <- c("cumulative_Lag 4", RR_cumul_lagged[5,1], RR_cumul_lagged[5,2], RR_cumul_lagged[5,3])
      RR_prov_cumulative_lag5_pollutant[i,] <- c("cumulative_Lag 5", RR_cumul_lagged[6,1], RR_cumul_lagged[6,2], RR_cumul_lagged[6,3])
      RR_prov_cumulative_lag6_pollutant[i,] <- c("cumulative_Lag 6", RR_cumul_lagged[7,1], RR_cumul_lagged[7,2], RR_cumul_lagged[7,3])
      RR_prov_cumulative_lag7_pollutant[i,] <- c("cumulative_Lag 7", RR_cumul_lagged[8,1], RR_cumul_lagged[8,2], RR_cumul_lagged[8,3])
      
    } else {
      warning(paste("Lagged RR matrix for province", valid_provinces[i], "does not have 8 rows. Skipping detailed lag assignment."))
    }
  })

##################################################
## SAVE RESULTS TO FILES FOR LATER COMPARISON ##
##################################################

# 1. Save Q-AIC values for each province
if (exists("qaic")) {
  qaic_df <- as.data.frame(qaic)
  qaic_df <- cbind(Province = rownames(qaic_df), qaic_df)
  write.csv(qaic_df, file.path(results_dir, "qaic_values.csv"), row.names = FALSE)
  cat("Saved: qaic_values.csv\n")
}

# 2. Save first-stage model coefficients and vcov matrices
# These are ymat_pollutant and list_vcov_pollutant
if (exists("ymat_pollutant") && exists("list_vcov_pollutant")) {
  save(ymat_pollutant, list_vcov_pollutant, file = file.path(results_dir, "first_stage_estimates.RData"))
  cat("Saved: first_stage_estimates.RData\n")
}

# 3. Save the multivariate meta-analysis results object
if (exists("mvmeta_results_pollutant")) {
  saveRDS(mvmeta_results_pollutant, file = file.path(results_dir, "mvmeta_results_pollutant.rds"))
  cat("Saved: mvmeta_results_pollutant.rds\n")
  
  # 3b. Save the summary of the multivariate meta-analysis to a text file
  capture.output(summary(mvmeta_results_pollutant), file = file.path(results_dir, "mvmeta_summary.txt"))
  cat("Saved: mvmeta_summary.txt\n")
}

# 4. Save the pooled overall cumulative Relative Risks from meta-analysis
if (exists("cpred_pooled_pollutant")) {
  pooled_overall_RR <- with(cpred_pooled_pollutant, cbind(RR = allRRfit, Lower = allRRlow, Upper = allRRhigh))
  # If there's only one estimate (e.g. at = 10), it might be a vector, convert to data.frame
  if (is.vector(pooled_overall_RR)) {
    pooled_overall_RR_df <- data.frame(t(pooled_overall_RR))
    colnames(pooled_overall_RR_df) <- c("RR", "Lower", "Upper")
    # Add the 'at' value if it's singular and known (here it is 10)
    pooled_overall_RR_df$at <- 10 # Assuming 'at = 10' as in your script
  } else {
    pooled_overall_RR_df <- as.data.frame(pooled_overall_RR)
  }
  write.csv(pooled_overall_RR_df, file.path(results_dir, "pooled_overall_cumulative_RR.csv"), row.names = FALSE)
  cat("Saved: pooled_overall_cumulative_RR.csv\n")
  
  # 4b. Save the pooled_crosspred object itself for more detailed inspection later
  saveRDS(cpred_pooled_pollutant, file = file.path(results_dir, "cpred_pooled_pollutant.rds"))
  cat("Saved: cpred_pooled_pollutant.rds\n")
}

# 5. Save the pooled lagged Relative Risks from meta-analysis
if (exists("lagged_RR_pooled_pollutant")) {
  # Ensure lagged_RR_pooled_pollutant is a data frame with proper column names
  lagged_RR_pooled_pollutant_df <- as.data.frame(lagged_RR_pooled_pollutant)
  if(ncol(lagged_RR_pooled_pollutant_df) == 3 && !identical(colnames(lagged_RR_pooled_pollutant_df), c("RR", "Lower", "Upper"))){
    colnames(lagged_RR_pooled_pollutant_df) <- c("RR", "Lower", "Upper") # Ensure correct names if not set by script
  }
  # Add a lag column
  lagged_RR_pooled_pollutant_df$Lag <- 0:(nrow(lagged_RR_pooled_pollutant_df)-1)
  # Reorder columns
  lagged_RR_pooled_pollutant_df <- lagged_RR_pooled_pollutant_df[, c("Lag", "RR", "Lower", "Upper")]
  
  write.csv(lagged_RR_pooled_pollutant_df, file.path(results_dir, "pooled_lagged_RR.csv"), row.names = FALSE)
  cat("Saved: pooled_lagged_RR.csv\n")
}

# 6. Save province-specific overall cumulative Relative Risks
if (exists("RR_overall_prov_pollutant")) {
  RR_overall_prov_pollutant_df <- as.data.frame(RR_overall_prov_pollutant)
  RR_overall_prov_pollutant_df <- cbind(Province = rownames(RR_overall_prov_pollutant_df), RR_overall_prov_pollutant_df)
  write.csv(RR_overall_prov_pollutant_df, file.path(results_dir, "province_specific_overall_RR.csv"), row.names = FALSE)
  cat("Saved: province_specific_overall_RR.csv\n")
}

# 7. Save province-specific lagged Relative Risks (Lag 0 to 7)
# Note: Your current script creates RR_prov_lagX_pollutant matrices where numeric RR values might be
# stored as characters due to c("Lag X", ...). For robust numeric analysis later,
# consider storing the 'Lag' column as numeric and other columns (RR, Lower, Upper) as numeric.
# The code below will save them as they are currently structured.

# First, combine them into a single data frame for easier handling (optional, but often useful)
all_prov_lagged_RR_list <- list()

for(lag_val in 0:lag_number) { # Use lag_number (which is 7 in your script)
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

for(lag_val in 0:lag_number) { # Use lag_number (which is 7 in your script)
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
    
    # Add to the list for a combined file (make sure column names are consistent)
    # The 'Lag' column in your current objects is like "Lag 0", "Lag 1", etc.
    # For a combined file, it's better to have a numeric lag column.
    # We'll extract the numeric part or assume based on lag_val.
    # For simplicity, we'll assume the structure is [Province, Lag_Char, RR, Lower, Upper]
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


# 8. Save the prediction basis and related parameters
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
  data_path = data_path,
  date_of_run = Sys.time()
  # Add any other parameters you vary, e.g., df for splines if you change them
  # ns_temp_df = 3, # example if you made this a variable
  # ns_humid_df = 3, # example
  # ns_pressure_df = 3, # example
  # arglag_df = 3 # example
)
saveRDS(params, file = file.path(results_dir, "run_parameters.rds"))
cat("Saved: run_parameters.rds\n")

cat("\nAll specified results have been saved.\n")

##################################################
##                  END OF SAVING               ##
##################################################