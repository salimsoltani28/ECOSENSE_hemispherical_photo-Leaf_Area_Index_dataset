# F07 Code: LT vs DHP LAI Comparison with GLMM Analysis

library(ggplot2)
library(dplyr)
library(patchwork)
library(lme4) 
library(MuMIn)

# Set working directory
setwd("/mnt/gsdata/users/lotz/LittervsLens/")

merged_data_cumulative <- read.csv("/mnt/gsdata/users/lotz/LittervsLens/results/LAI_results_wood_removal_full.csv", header = TRUE)

# Define plot IDs without conifer influence
plot_ids_wo_conifer <- c('LT11', 'LT13', 'LT14', 'LT23', 'LT24', 'LT33', 'LT34', 'LT41', 'LT44', 'LT51', 'LT52', 'LT53', 'LT62', 'LT61', 'LT63')

# Define the time spans for each phase
begin_start <- as.Date("2024-09-20")
begin_end <- as.Date("2024-10-15")
peak_start <- as.Date("2024-10-16")
peak_end <- as.Date("2024-11-12")
end_start <- as.Date("2024-11-13")
end_end <- as.Date("2024-12-16")

# Add the `phase` column based on the defined timespans
merged_data_cumulative <- merged_data_cumulative %>%
  mutate(
    phase = case_when(
      date >= begin_start & date <= begin_end ~ "onset",
      date >= peak_start & date <= peak_end ~ "peak",
      date >= end_start & date <= end_end ~ "end",
      TRUE ~ "other"
    )
  )

plot_data <- merged_data_cumulative %>%
  filter(endVZA == 20) %>%
  filter(ID %in% plot_ids_wo_conifer) %>%
  filter(phase != "other")  # Remove "other" phase for GLMM

# Define custom phase colors
custom_colors <- c("onset" = "#31688E", "peak" = "#35B779", "end" = "#C6801A")

# ===== FIT GLMM MODEL =====
# Fit the GLMM: LT LAI ~ LAI_LXG1_nowood * phase + (1|ID)
glmm_model <- lmer(LTLAI ~ LAI_LXG1_nowood * phase + (1|ID), data = plot_data)

# Print model summary
cat("GLMM Model Summary:\n")
summary(glmm_model)

# Generate predictions
plot_data$predicted_LTLAI <- predict(glmm_model, plot_data)

# Calculate R² for GLMM (marginal and conditional)

r2_glmm <- r.squaredGLMM(glmm_model)
cat("GLMM R² - Marginal (fixed effects):", round(r2_glmm[1], 3), "\n")
cat("GLMM R² - Conditional (fixed + random):", round(r2_glmm[2], 3), "\n")

# Function to create the original plot
create_plot <- function(plot_data, x_var, y_var, method_name, title_label) {
  plot_data$phase <- factor(plot_data$phase, levels = c("onset", "peak", "end"))
  
  x_range <- range(plot_data$LTLAI, na.rm = TRUE)
  y_range <- range(plot_data[[y_var]], na.rm = TRUE)
  max_range <- max(c(x_range, y_range))
  
  lm_model <- lm(as.formula(paste(y_var, "~ LTLAI")), data = plot_data)
  intercept <- coef(lm_model)["(Intercept)"]
  slope <- coef(lm_model)["LTLAI"]
  r2 <- summary(lm_model)$r.squared
  
  eq_label <- sprintf("y = %.2f + %.2fx,\nR² = %.2f", intercept, slope, r2)
  
  # Calculate MAE by phase for this plot
  plot_data$abs_error <- abs(plot_data$LTLAI - plot_data[[y_var]])
  mae_by_phase <- plot_data %>%
    group_by(phase) %>%
    summarise(MAE = mean(abs_error, na.rm = TRUE), .groups = 'drop')
  
  # Create MAE label
  mae_text <- paste(
    paste0("onset: MAE = ", sprintf("%.2f", mae_by_phase$MAE[mae_by_phase$phase == "onset"])),
    paste0("peak: MAE = ", sprintf("%.2f", mae_by_phase$MAE[mae_by_phase$phase == "peak"])),
    paste0("end: MAE = ", sprintf("%.2f", mae_by_phase$MAE[mae_by_phase$phase == "end"])),
    sep = "\n"
  )
  
  p <- ggplot(data = plot_data, aes(x = LTLAI, y = .data[[y_var]], color = phase)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", color = "darkgrey", se = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    geom_text(
      x = max_range * 0.05, 
      y = max_range * 0.95, 
      label = eq_label, 
      color = "black", 
      size = 6, 
      hjust = 0
    ) +
    geom_text(
      x = max_range * 0.95, 
      y = max_range * 0.25, 
      label = mae_text, 
      color = "black", 
      size = 5, 
      hjust = 1
    ) +
    coord_fixed(ratio = 1) +  
    xlim(0, max_range) + ylim(0, max_range) +  
    labs(
      x = "LT LAI",
      y = "DHP LAI",
      color = NULL
    ) +
    scale_color_manual(values = custom_colors) +
    theme_minimal() +  
    theme(
      text = element_text(family = "Helvetica", size = 25),
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 25),
      axis.title = element_text(size = 25, face = "bold"),
      axis.text = element_text(size = 25),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  return(p)
}

# Function to create GLMM prediction plot
create_glmm_plot <- function(plot_data, title_label) {
  plot_data$phase <- factor(plot_data$phase, levels = c("onset", "peak", "end"))
  
  x_range <- range(plot_data$LTLAI, na.rm = TRUE)
  y_range <- range(plot_data$predicted_LTLAI, na.rm = TRUE)
  max_range <- max(c(x_range, y_range))
  
  # Calculate overall R² for predicted vs observed
  lm_pred <- lm(predicted_LTLAI ~ LTLAI, data = plot_data)
  intercept <- coef(lm_pred)["(Intercept)"]
  slope <- coef(lm_pred)["LTLAI"]
  r2_pred <- summary(lm_pred)$r.squared
  
  eq_label <- sprintf("y = %.2f + %.2fx,\nR² = %.2f", intercept, slope, r2_pred)
  
  # Calculate MAE by phase for GLMM predictions
  plot_data$abs_error_glmm <- abs(plot_data$LTLAI - plot_data$predicted_LTLAI)
  mae_by_phase <- plot_data %>%
    group_by(phase) %>%
    summarise(MAE = mean(abs_error_glmm, na.rm = TRUE), .groups = 'drop')
  
  # Create MAE label
  mae_text <- paste(
    paste0("onset: MAE = ", sprintf("%.2f", mae_by_phase$MAE[mae_by_phase$phase == "onset"])),
    paste0("peak: MAE = ", sprintf("%.2f", mae_by_phase$MAE[mae_by_phase$phase == "peak"])),
    paste0("end: MAE = ", sprintf("%.2f", mae_by_phase$MAE[mae_by_phase$phase == "end"])),
    sep = "\n"
  )
  
  p <- ggplot(data = plot_data, aes(x = LTLAI, y = predicted_LTLAI, color = phase)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", color = "darkgrey", se = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    geom_text(
      x = max_range * 0.05, 
      y = max_range * 0.95, 
      label = eq_label, 
      color = "black", 
      size = 6, 
      hjust = 0
    ) +
    geom_text(
      x = max_range * 0.95, 
      y = max_range * 0.25, 
      label = mae_text, 
      color = "black", 
      size = 5, 
      hjust = 1
    ) +
    coord_fixed(ratio = 1) +  
    xlim(0, max_range) + ylim(0, max_range) +  
    labs(
      x = "LT LAI",
      y = "Predicted DHP LAI",
      color = NULL
    ) +
    scale_color_manual(values = custom_colors) +
    theme_minimal() +  
    theme(
      text = element_text(family = "Helvetica", size = 25),
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 25),
      axis.title = element_text(size = 25, face = "bold"),
      axis.text = element_text(size = 25),
      plot.margin = margin(20, 20, 20, 20),
      legend.position = "none"  # Remove legend from individual plots
    )
  
  return(p)
}

# Create both plots
plot_a <- create_plot(plot_data, "LTLAI", "LAI_LXG1_nowood", "LAI.LXG1", "a")
plot_b <- create_glmm_plot(plot_data, "b")

# Combine plots side by side with legend below and spacing
combined_plot <- plot_a + plot_spacer() + plot_b + 
  plot_layout(ncol = 3, widths = c(1, 0.05, 1), guides = "collect") &
  theme(legend.position = "bottom")

# Save individual and combined plots
ggsave("results/LTvsDHP_raw.png", plot_a, width = 10, height = 8, dpi = 300)
ggsave("results/LTvsDHP_glmm.png", plot_b, width = 10, height = 8, dpi = 300)
ggsave("results/LTvsDHP_combined.png", combined_plot, width = 16, height = 8, dpi = 300)

# ===== CALCULATE ABSOLUTE ERRORS BY PHASE =====
# Calculate absolute errors for raw data
plot_data$abs_error_raw <- abs(plot_data$LTLAI - plot_data$LAI_LXG1_nowood)

# Calculate absolute errors for GLMM predictions
plot_data$abs_error_glmm <- abs(plot_data$LTLAI - plot_data$predicted_LTLAI)

# Calculate mean absolute errors by phase
mae_by_phase <- plot_data %>%
  group_by(phase) %>%
  summarise(
    n_obs = n(),
    MAE_raw = mean(abs_error_raw, na.rm = TRUE),
    MAE_glmm = mean(abs_error_glmm, na.rm = TRUE),
    RMSE_raw = sqrt(mean((LTLAI - LAI_LXG1_nowood)^2, na.rm = TRUE)),
    RMSE_glmm = sqrt(mean((LTLAI - predicted_LTLAI)^2, na.rm = TRUE)),
    .groups = 'drop'
  )

# Print absolute error results
cat("\n=== ABSOLUTE ERRORS BY PHASE ===\n")
print(mae_by_phase)

# Overall errors
overall_mae_raw <- mean(plot_data$abs_error_raw, na.rm = TRUE)
overall_mae_glmm <- mean(plot_data$abs_error_glmm, na.rm = TRUE)
overall_rmse_raw <- sqrt(mean((plot_data$LTLAI - plot_data$LAI_LXG1_nowood)^2, na.rm = TRUE))
overall_rmse_glmm <- sqrt(mean((plot_data$LTLAI - plot_data$predicted_LTLAI)^2, na.rm = TRUE))

cat("\n=== OVERALL ERRORS ===\n")
cat("Raw data - MAE:", round(overall_mae_raw, 3), ", RMSE:", round(overall_rmse_raw, 3), "\n")
cat("GLMM predictions - MAE:", round(overall_mae_glmm, 3), ", RMSE:", round(overall_rmse_glmm, 3), "\n")

# Print GLMM results
cat("\n=== GLMM Results ===\n")
cat("Model: LTLAI ~ LAI_LXG1_nowood * phase + (1|ID)\n")
cat("Marginal R² (fixed effects only):", round(r2_glmm[1], 3), "\n")
cat("Conditional R² (fixed + random effects):", round(r2_glmm[2], 3), "\n")