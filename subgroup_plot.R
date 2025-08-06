# This file is for plotting the subgroup analysis results for easier visualization.

library(ggplot2)
library(dplyr)

# Set the plotting function
plot_across_subgroup <- function(
    df,
    subgroup_col,
    max_lag = 4) {
  # The '{{}}' syntax is used to correctly handle the subgroup_col argument
  rr_line_plot <- df %>%
    mutate(across(
      c("Lag", "RR", "Lower", "Upper"),
      as.numeric
    )) %>%
    filter(Lag < max_lag) %>%
    ggplot(
      # Map color, fill, and group to your subgroup column
      aes(
        x = Lag,
        y = RR,
        ymin = Lower,
        ymax = Upper,
        color = {{ subgroup_col }},
        fill = {{ subgroup_col }},
        group = {{ subgroup_col }}
      )
    ) +
    # Use alpha to make ribbons transparent, the fill color is from aes()
    geom_errorbar(width = 0.2, linewidth = 0.8, alpha = 0.4, position = "dodge") +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
    scale_x_continuous(breaks = 0:(max_lag - 1)) +
    labs(
      # You may want to make the title and subtitle more dynamic
      subtitle = "RR association by subgroup",
      title = "RR association with a 10-unit increase in PM2.5",
      x = "Lag",
      y = "RR and 95% CI",
      # These lines will rename the legend
      color = "Subgroup",
      fill = "Subgroup"
    ) +
    theme_minimal(base_size = 14) + # A clean theme to start with
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right" # Move legend to the bottom
    )

  return(rr_line_plot)
}
# ====================== RUN THE CODE =====================
#
# set the output dir
BASE_OUTPUT_DIR <- "output/regional_subgroup_analysis_30jul"

# load the subgroup analysis data

# Gender
pooled_rr_male <- read.csv("/Users/tiger/Projects/hf_pm_analysis/output/30jul/results_PM2.5-hf_prim_male-7lag/pooled_lagged_RR.csv")
pooled_rr_male$gender <- "male"
pooled_rr_female <- read.csv("/Users/tiger/Projects/hf_pm_analysis/output/30jul/results_PM2.5-hf_prim_female-7lag/pooled_lagged_RR.csv")
pooled_rr_female$gender <- "female"

# Age 65+
pooled_rr_under65 <- read.csv("/Users/tiger/Projects/hf_pm_analysis/output/30jul/results_PM2.5-hf_prim_under65-7lag/pooled_lagged_RR.csv")
pooled_rr_under65$age_group <- "under65"
pooled_rr_over65 <- read.csv("/Users/tiger/Projects/hf_pm_analysis/output/30jul/results_PM2.5-hf_prim_over65-7lag/pooled_lagged_RR.csv")
pooled_rr_over65$age_group <- "over65"

# Gender DF plot
pooled_rr_gender <- rbind(pooled_rr_male, pooled_rr_female)
write.csv(pooled_rr_gender, file.path(BASE_OUTPUT_DIR, "pooled_rr_gender.csv"))

gender_subplot <- plot_across_subgroup(pooled_rr_gender, gender)
ggsave(
  file.path(BASE_OUTPUT_DIR, "gender_RR_lagged_plot.png"),
  plot = gender_subplot,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(BASE_OUTPUT_DIR, "gender_RR_lagged_plot.svg"),
  plot = gender_subplot,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

# Age DF plot
pooled_rr_age <- rbind(pooled_rr_under65, pooled_rr_over65)
write.csv(pooled_rr_age, file.path(BASE_OUTPUT_DIR, "pooled_rr_age.csv"))

age_subplot <- plot_across_subgroup(pooled_rr_age, age_group)
ggsave(
  file.path(BASE_OUTPUT_DIR, "age_RR_lagged_plot.png"),
  plot = age_subplot,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(BASE_OUTPUT_DIR, "age_RR_lagged_plot.svg"),
  plot = age_subplot,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)
