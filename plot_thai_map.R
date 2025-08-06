set.seed(123)
setwd("~/Projects/hf_pm_analysis")

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

# =================

# Load necessary libraries
library(ggplot2)
library(sf)
library(dplyr)
library(geodata)

# Download province-level (level=1) data for Thailand
thailand_provinces <- gadm(country = "THA", level = 1, path = tempdir()) |>
  st_as_sf()

# Create total hospitalization count by province
avg_hf_count <- dat %>%
  group_by(prov_name) %>%
  summarise(
    hf_count = mean(hf_prim, na.rm = TRUE), # Average hospitalization count
    .groups = "drop"
  )

##  Merge Spatial Data with Your Data
thailand_map_with_hf <- thailand_provinces |>
  left_join(avg_hf_count, by = c("NAME_1" = "prov_name"))

##  Plot the Heat Map (Choropleth)
ggplot(data = thailand_map_with_hf) +
  geom_sf(aes(fill = hf_count), color = "white", size = 0.1) +
  scale_fill_gradient(
    low = "white",
    high = "blue",
    na.value = "darkgrey"
  ) + # Use a nice color scale
  labs(
    title = "Average Daily Hospitalization by Heart Failure by Province in Thailand",
    fill = "Hospitalization Count",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_blank(), # Removes both x and y axis titles
    axis.text = element_blank(), # Removes the axis text (numbers)
    axis.ticks = element_blank(), # Removes the axis tick marks
    panel.grid = element_blank()
  )
ggsave(
  filename = file.path(results_dir, "thai_map_hf.png"),
  height = 8, width = 10, dpi = 300, units = "in",
  bg = "white"
)
ggsave(
  filename = file.path(results_dir, "thai_map_hf.svg"),
  height = 8, width = 10, dpi = 300, units = "in",
  bg = "white"
)

# Save the aggregated data to csv
thailand_map_with_hf_save <- thailand_map_with_hf %>%
  select(NAME_1, VARNAME_1, NL_NAME_1, hf_count)

write.csv(thailand_map_with_hf_save,
  file = file.path(results_dir, "thai_map_hf.csv"),
  row.names = FALSE
)

# =============================================
# Create average map of PM2.5 by provinces
avg_prov_poll <- dat %>%
  group_by(prov_name) %>%
  summarise(
    average_pm25 = mean(PM2.5, na.rm = TRUE), # Average hospitalization count
    .groups = "drop"
  )

##  Merge Spatial Data with Your Data
thailand_map_with_pm25 <- thailand_provinces |>
  left_join(avg_prov_poll, by = c("NAME_1" = "prov_name"))

##  Plot the Heat Map (Choropleth)
ggplot(data = thailand_map_with_pm25) +
  geom_sf(aes(fill = average_pm25), color = "white", size = 0.1) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    na.value = "darkgrey"
  ) + # Use a nice color scale
  labs(
    title = "Avearge PM2.5 levels by Province in Thailand",
    fill = "Average PM2.5",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_blank(), # Removes both x and y axis titles
    axis.text = element_blank(), # Removes the axis text (numbers)
    axis.ticks = element_blank(), # Removes the axis tick marks
    panel.grid = element_blank()
  )
ggsave(
  filename = file.path(results_dir, "thai_map_pm25.png"),
  height = 8, width = 10, dpi = 300, units = "in",
  bg = "white"
)
ggsave(
  filename = file.path(results_dir, "thai_map_pm25.svg"),
  height = 8, width = 10, dpi = 300, units = "in",
  bg = "white"
)

# Save the aggregated data to csv
thailand_map_with_pm25_save <- thailand_map_with_pm25 %>%
  select(NAME_1, VARNAME_1, NL_NAME_1, average_pm25) %>%
  select(-geometry)

write.csv(thailand_map_with_pm25_save,
  file = file.path(results_dir, "thai_map_pm25.csv"),
  row.names = FALSE
)
