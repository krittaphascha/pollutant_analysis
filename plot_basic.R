set.seed(123)
setwd("~/Projects/hf_pm_analysis")

library(xtable)
library(dplyr)
library(zoo)
library(ggplot2)
library(patchwork)
library(reshape2)
library(tidyr)

dat <- read.csv("data/daily_all_weather_30jul.csv")

# properly format date
dat$date_start <- as.Date(dat$date_start, "%Y-%m-%d")
dat$dow <- as.factor(weekdays(dat$date_start))
dat$month <- as.factor(months(dat$date_start))
dat$year <- as.factor(format(dat$date_start, format = "%Y"))

# get strata for case crossover analysis
dat$stratum <- as.factor(dat$year:dat$month:dat$dow)

results_dir <- "output/demographic_pollutants"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

invalid_prov <- c("Saraburi", "Samut Sakhon", "Nonthaburi")
dat <- dat %>%
  filter(!prov_name %in% invalid_prov)

# ==============================
# Plotting all basic demographic data and pollutant data

# 1. Table 1
# Select the numeric columns you want to summarize
cols_to_summarize <- c(
  "O3", "PM2.5", "PM10", "NO2", "NOX", "CO", "SO2",
  "hf_prim", "humidity", "temperature", "pressure"
)

# Apply the summary function to each column and transpose the result
tab1 <- t(apply(dat[, cols_to_summarize], 2, function(x) {
  c(
    mean(x, na.rm = TRUE),
    sd(x, na.rm = TRUE),
    min(x, na.rm = TRUE),
    quantile(x, c(0.25, 0.50, 0.75), na.rm = TRUE),
    max(x, na.rm = TRUE)
  )
}))

# Assign column and row names
dimnames(tab1) <- list(
  c(
    "Ozone (O3)", "Fine Particulate Matter (PM2.5)", "Particulate Matter (PM10)",
    "Nitrogen Dioxide (NO2)", "Nitrogen Oxides (NOX)", "Carbon Monoxide (CO)",
    "Sulfur Dioxide (SO2)", "HF Hospitalization", "Humidity (%)",
    "Temperature (Celsius)", "Air Pressure (hPa)"
  ),
  c("Mean", "SD", "Min", "25%", "Median", "75%", "Max")
)

# Print the resulting table
capture.output(
  {
    print(tab1, digits = 2)
    print(xtable(tab1, digits = 1, align = "l|rrrrrrr"))
  },
  file = file.path(results_dir, "table1_summary_statistics.txt")
)


# 2. Spearman correlation matrix between each pollutant and HF hospitalization
correlation_matrix <- cor(dat[, c("O3", "PM2.5", "PM10", "NO2", "NOX", "CO", "SO2")],
  use = "pairwise.complete.obs", method = "spearman"
)
correlation_matrix[upper.tri(correlation_matrix)] <- NA # Set upper triangle to NA for better visualization
# Print the correlation matrix
capture.output(
  {
    print(correlation_matrix, digits = 2)
  },
  file = file.path(results_dir, "spearman_correlation_matrix.txt")
)

# create heatmap correlation matrix
melted_correlation_matrix <- melt(correlation_matrix)
ggplot(
  data = melted_correlation_matrix,
  aes(x = Var1, y = Var2, fill = value)
) +
  geom_tile() +
  scale_fill_distiller() +
  geom_text(aes(Var1, Var2, label = round(value, 2)),
    color = "black", size = 4
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)
  ) +
  ggtitle("Spearman Correlation Matrix") +
  theme(plot.title = element_text(hjust = 0.7))
ggsave(file.path(results_dir, "correlation_matrix_heatmap.png"), width = 8, height = 6, dpi = 300)
ggsave(file.path(results_dir, "correlation_matrix_heatmap.svg"), width = 8, height = 6, dpi = 300)

# Plot Heatmap to show Pollutant distribution value in each province
monthly_pol <- dat %>%
  group_by(prov_name, month) %>%
  summarise(across(c("O3", "PM2.5", "PM10", "NO2", "NOX", "CO", "SO2"), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(month = factor(month, levels = month.name))

ggplot(monthly_pol, aes(x = month, y = prov_name, fill = PM2.5)) +
  geom_tile(color = "white") + # geom_tile is the function for heatmaps
  scale_fill_viridis_c(name = "Monthly Avg. PM2.5 (μg/m³)") +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "bottom") +
  labs(
    title = "Monthly PM2.5 Concentration by Province",
    x = "Month",
    y = "Province"
  ) +
  theme_minimal() + # A clean theme
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1.1, vjust = 0.95), # Angle the month text
    panel.grid = element_blank() # Remove grid lines
  )
ggsave(
  file.path(
    results_dir,
    paste0("heatmap_monthly_avg_", pollutant_col_name, ".png")
  ),
  width = 8,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(
    results_dir,
    paste0("heatmap_monthly_avg_", pollutant_col_name, ".svg")
  ),
  width = 8,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)


# Plot average PM2.5 level by date

daily_avg <- dat %>%
  group_by(date_start) %>%
  summarise(
    avg_pm25 = mean(PM2.5, na.rm = TRUE),
    sum_hf_prim = sum(hf_prim, na.rm = TRUE),
    avg_temp = mean(temperature, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pm25_ma7 = rollmean(
      avg_pm25, 7,
      fill = NA, align = "right"
    ),
    temp_ma7 = rollmean(
      avg_temp, 7,
      fill = NA, align = "right"
    ),
    sum_hf_prim_ma7 = rollmean(
      sum_hf_prim, 7,
      fill = NA, align = "right"
    )
  )

daily_pm25 <- ggplot(daily_avg, aes(x = date_start, y = pm25_ma7)) +
  geom_line(color = "blue") +
  labs(
    title = "Daily Average PM2.5 Levels (7-day Moving Average)",
    x = "Date",
    y = "Average PM2.5 (μg/m³)"
  ) +
  scale_x_date(
    date_breaks = "4 months",
    date_labels = "%b %Y"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
  )
ggsave(
  file.path(results_dir, "daily_avg_pm25.png"),
  width = 12,
  height = 4,
  unit = "in",
  bg = "white",
  dpi = 300
)
ggsave(
  file.path(results_dir, "daily_avg_pm25.svg"),
  width = 12,
  height = 4,
  unit = "in",
  bg = "white",
  dpi = 300
)

daily_temp <- ggplot(daily_avg, aes(x = date_start, y = temp_ma7)) +
  geom_line(color = "darkgreen") +
  labs(
    title = "Daily Average Temperature Levels (7-day Moving Average)",
    x = "Date",
    y = "Average Temperature (Celcius)"
  ) +
  scale_x_date(
    date_breaks = "4 months",
    date_labels = "%b %Y"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
  )
ggsave(
  file.path(results_dir, "daily_avg_temp.png"),
  width = 12,
  height = 4,
  unit = "in",
  bg = "white",
  dpi = 300
)
ggsave(
  file.path(results_dir, "daily_avg_temp.svg"),
  width = 12,
  height = 4,
  unit = "in",
  bg = "white",
  dpi = 300
)

daily_hf <- ggplot(daily_avg, aes(x = date_start, y = sum_hf_prim_ma7)) +
  geom_line(color = "red") +
  labs(
    title = "Daily Total HF Hospitalizations (7-day Moving Average)",
    x = "Date",
    y = "Sum of HF Hospitalizations"
  ) +
  scale_x_date(
    date_breaks = "4 months",
    date_labels = "%b %Y"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
  )
ggsave(
  file.path(results_dir, "daily_sum_hf_prim.png"),
  width = 12,
  height = 4,
  unit = "in",
  bg = "white",
  dpi = 300
)
ggsave(
  file.path(results_dir, "daily_sum_hf_prim.svg"),
  width = 12,
  height = 4,
  unit = "in",
  bg = "white",
  dpi = 300
)

# =======================
# Plot time series data of Pm2.5 level and Total Hospitalization Count
# Process data and create the plot

# Process the data for plotting
plot_data_hf_pm2.5 <- dat %>%
  select(date_start, hf_prim, PM2.5) %>%
  group_by(date_start) %>%
  summarise(
    total_hf = sum(hf_prim, na.rm = TRUE),
    average_pm25 = mean(PM2.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_hf = rollmean(total_hf, 7, fill = NA, align = "right"),
    average_pm25 = rollmean(average_pm25, 7, fill = NA, align = "right")
  ) %>%
  drop_na()

# Determine the scaling factor for the secondary y-axis.
# This is used to overlay the two series.
scaling_factor_hf <- max(plot_data_hf_pm2.5$total_hf) / max(plot_data_hf_pm2.5$average_pm25)

# Create the dual-axis plot
plot_dual_axis_hf_pm <- plot_data_hf_pm2.5 %>%
  ggplot(aes(x = date_start)) +
  # Plot the first series: Total Hospitalizations
  geom_line(aes(y = total_hf, color = "Total Hospitalizations"), size = 0.4) +
  # Plot the second series: Average PM2.5, scaled to the primary y-axis
  geom_line(aes(y = average_pm25 * scaling_factor_hf, color = "Average PM2.5"), size = 0.4, alpha = 0.5) +
  # Create the secondary y-axis by reversing the scaling transformation
  scale_y_continuous(
    name = "Total Hospitalization Count",
    sec.axis = sec_axis(~ . / scaling_factor_hf, name = "Average PM2.5 (µg/m³)")
  ) +
  # Add informative labels and a title
  labs(
    title = "Hospitalizations with Heart Failure and PM2.5 Levels Over Time",
    subtitle = "Daily Aggregated Values with 7-day Moving Average",
    x = "Date",
    y = "Total Hospitalization Count"
  ) +
  # Customize colors and legend title
  scale_color_manual(
    name = "Metric",
    values = c("Total Hospitalizations" = "blue", "Average PM2.5" = "firebrick")
  ) +
  scale_x_date(
    breaks = scales::date_breaks("4 months"),
    date_labels = "%b %Y"
  ) +
  # Apply a clean theme and adjust legend position
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.y.left = element_text(color = "blue", face = "bold"),
    axis.text.y.left = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "firebrick"),
    axis.text.y.right = element_text(color = "firebrick"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(results_dir, "hf_pm25_dual_axis_time_series.png"),
  plot = plot_dual_axis_hf_pm,
  height = 6, width = 12, dpi = 300, units = "in",
  bg = "white"
)
ggsave(
  filename = file.path(results_dir, "hf_pm25_dual_axis_time_series.svg"),
  plot = plot_dual_axis_hf_pm,
  height = 6, width = 12, dpi = 300, units = "in",
  bg = "white"
)

# Process the data for plotting Atrial Fibrillation
plot_data_af_pm2.5 <- dat %>%
  select(date_start, af_prim, PM2.5) %>%
  group_by(date_start) %>%
  summarise(
    total_af = sum(af_prim, na.rm = TRUE),
    average_pm25 = mean(PM2.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_af = rollmean(total_af, 7, fill = NA, align = "right"),
    average_pm25 = rollmean(average_pm25, 7, fill = NA, align = "right")
  ) %>%
  drop_na()

# Determine the scaling factor for the secondary y-axis.
# This is used to overlay the two series.
scaling_factor_af <- max(plot_data_af_pm2.5$total_af) / max(plot_data_af_pm2.5$average_pm25)

# Create the dual-axis plot
plot_dual_axis_af_pm <- plot_data_af_pm2.5 %>%
  ggplot(aes(x = date_start)) +
  # Plot the first series: Total Hospitalizations
  geom_line(aes(y = total_af, color = "Total Hospitalizations"), size = 0.4) +
  # Plot the second series: Average PM2.5, scaled to the primary y-axis
  geom_line(aes(y = average_pm25 * scaling_factor_af, color = "Average PM2.5"), size = 0.4, alpha = 0.5) +
  # Create the secondary y-axis by reversing the scaling transformation
  scale_y_continuous(
    name = "Total Hospitalization Count",
    sec.axis = sec_axis(~ . / scaling_factor_af, name = "Average PM2.5 (µg/m³)")
  ) +
  # Add informative labels and a title
  labs(
    title = "Hospitalizations with Atrial Fibrillation and PM2.5 Levels Over Time",
    subtitle = "Daily Aggregated Values with 7-day Moving Average",
    x = "Date",
    y = "Total Hospitalization Count"
  ) +
  # Customize colors and legend title
  scale_color_manual(
    name = "Metric",
    values = c("Total Hospitalizations" = "blue", "Average PM2.5" = "firebrick")
  ) +
  scale_x_date(
    breaks = scales::date_breaks("4 months"),
    date_labels = "%b %Y"
  ) +
  # Apply a clean theme and adjust legend position
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.y.left = element_text(color = "blue", face = "bold"),
    axis.text.y.left = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "firebrick"),
    axis.text.y.right = element_text(color = "firebrick"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(results_dir, "af_pm25_dual_axis_time_series.png"),
  plot = plot_dual_axis_af_pm,
  height = 6, width = 12, dpi = 300, units = "in",
  bg = "white"
)
ggsave(
  filename = file.path(results_dir, "af_pm25_dual_axis_time_series.svg"),
  plot = plot_dual_axis_af_pm,
  height = 6, width = 12, dpi = 300, units = "in",
  bg = "white"
)

# Process the data for plotting
plot_data_hf_all_pm2.5 <- dat %>%
  select(date_start, hf, PM2.5) %>%
  group_by(date_start) %>%
  summarise(
    total_hf = sum(hf, na.rm = TRUE),
    average_pm25 = mean(PM2.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_hf = rollmean(total_hf, 7, fill = NA, align = "right"),
    average_pm25 = rollmean(average_pm25, 7, fill = NA, align = "right")
  ) %>%
  drop_na()

# Determine the scaling factor for the secondary y-axis.
# This is used to overlay the two series.
scaling_factor_hf <- max(plot_data_hf_all_pm2.5$total_hf) / max(plot_data_hf_all_pm2.5$average_pm25)

# Create the dual-axis plot
plot_dual_axis_hf_all_pm <- plot_data_hf_all_pm2.5 %>%
  ggplot(aes(x = date_start)) +
  # Plot the first series: Total Hospitalizations
  geom_line(aes(y = total_hf, color = "Total Hospitalizations"), size = 0.4) +
  # Plot the second series: Average PM2.5, scaled to the primary y-axis
  geom_line(aes(y = average_pm25 * scaling_factor_hf, color = "Average PM2.5"), size = 0.4, alpha = 0.5) +
  # Create the secondary y-axis by reversing the scaling transformation
  scale_y_continuous(
    name = "Total Hospitalization Count",
    sec.axis = sec_axis(~ . / scaling_factor_hf, name = "Average PM2.5 (µg/m³)")
  ) +
  # Add informative labels and a title
  labs(
    title = "All Hospitalizations associated with Heart Failure and PM2.5 Levels Over Time",
    subtitle = "Daily Aggregated Values with 7-day Moving Average",
    x = "Date",
    y = "Total Hospitalization Count"
  ) +
  # Customize colors and legend title
  scale_color_manual(
    name = "Metric",
    values = c("Total Hospitalizations" = "blue", "Average PM2.5" = "firebrick")
  ) +
  scale_x_date(
    breaks = scales::date_breaks("4 months"),
    date_labels = "%b %Y"
  ) +
  # Apply a clean theme and adjust legend position
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.y.left = element_text(color = "blue", face = "bold"),
    axis.text.y.left = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "firebrick"),
    axis.text.y.right = element_text(color = "firebrick"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(results_dir, "hf_all_pm25_dual_axis_time_series.png"),
  plot = plot_dual_axis_hf_all_pm,
  height = 6, width = 12, dpi = 300, units = "in",
  bg = "white"
)
ggsave(
  filename = file.path(results_dir, "hf_all_pm25_dual_axis_time_series.svg"),
  plot = plot_dual_axis_hf_all_pm,
  height = 6, width = 12, dpi = 300, units = "in",
  bg = "white"
)

# Process the data for plotting Atrial Fibrillation
plot_data_af_all_pm2.5 <- dat %>%
  select(date_start, af, PM2.5) %>%
  group_by(date_start) %>%
  summarise(
    total_af = sum(af, na.rm = TRUE),
    average_pm25 = mean(PM2.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_af = rollmean(total_af, 7, fill = NA, align = "right"),
    average_pm25 = rollmean(average_pm25, 7, fill = NA, align = "right")
  ) %>%
  drop_na()

# Determine the scaling factor for the secondary y-axis.
# This is used to overlay the two series.
scaling_factor_af <- max(plot_data_af_all_pm2.5$total_af) / max(plot_data_af_all_pm2.5$average_pm25)

# Create the dual-axis plot
plot_dual_axis_af_all_pm <- plot_data_af_all_pm2.5 %>%
  ggplot(aes(x = date_start)) +
  # Plot the first series: Total Hospitalizations
  geom_line(aes(y = total_af, color = "Total Hospitalizations"), size = 0.4) +
  # Plot the second series: Average PM2.5, scaled to the primary y-axis
  geom_line(aes(y = average_pm25 * scaling_factor_af, color = "Average PM2.5"), size = 0.4, alpha = 0.5) +
  # Create the secondary y-axis by reversing the scaling transformation
  scale_y_continuous(
    name = "Total Hospitalization Count",
    sec.axis = sec_axis(~ . / scaling_factor_af, name = "Average PM2.5 (µg/m³)")
  ) +
  # Add informative labels and a title
  labs(
    title = "All Hospitalizations associated with Atrial Fibrillation and PM2.5 Levels Over Time",
    subtitle = "Daily Aggregated Values with 7-day Moving Average",
    x = "Date",
    y = "Total Hospitalization Count"
  ) +
  # Customize colors and legend title
  scale_color_manual(
    name = "Metric",
    values = c("Total Hospitalizations" = "blue", "Average PM2.5" = "firebrick")
  ) +
  scale_x_date(
    breaks = scales::date_breaks("4 months"),
    date_labels = "%b %Y"
  ) +
  # Apply a clean theme and adjust legend position
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.y.left = element_text(color = "blue", face = "bold"),
    axis.text.y.left = element_text(color = "blue"),
    axis.title.y.right = element_text(color = "firebrick"),
    axis.text.y.right = element_text(color = "firebrick"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(results_dir, "af_all_pm25_dual_axis_time_series.png"),
  plot = plot_dual_axis_af_all_pm,
  height = 6, width = 12, dpi = 300, units = "in",
  bg = "white"
)
ggsave(
  filename = file.path(results_dir, "af_all_pm25_dual_axis_time_series.svg"),
  plot = plot_dual_axis_af_all_pm,
  height = 6, width = 12, dpi = 300, units = "in",
  bg = "white"
)
