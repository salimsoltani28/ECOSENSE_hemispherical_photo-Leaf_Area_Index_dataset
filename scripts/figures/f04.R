# F04 Code: VZA Analysis with Clumping Comparison and Woody Component Removal

library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(viridis)
library(cowplot)

# Set working directory
setwd("/mnt/gsdata/users/lotz/LittervsLens/")

# Filter for specific LT IDs
selected_ids <- c('LT11', 'LT13', 'LT14', 'LT23', 'LT24', 'LT33', 'LT34', 'LT41', 'LT44', 'LT51', 'LT52', 'LT53', 'LT62', 'LT61', 'LT63')

# Read the datasets
lt_results <- read.csv("results/LT_results.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
lai_results <- read.csv("results/LAI_results.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")
gap_results <- read.csv("results/LAI_results_gapfrac.csv", header = TRUE, stringsAsFactors = FALSE, sep = ",")

plot_dir <- "plots/"

# Filter and prepare LT data
lt_results <- lt_results %>%
  filter(plot %in% selected_ids) %>%
  mutate(date = as.Date(as.character(date)))

# Prepare LAI data for all VZA
lai_all_vza <- lai_results %>%
  filter(ID %in% selected_ids) %>%
  mutate(date = as.Date(as.character(date), format = "%Y%m%d")) %>%
  select(ID, date, endVZA, Le, LX, LXG1, LXG2, MTA_ell)

# Prepare gap fraction data for all VZA
gap_all_vza <- gap_results %>%
  filter(ID %in% selected_ids) %>%
  mutate(date = as.Date(as.character(date), format = "%Y%m%d")) %>%
  select(ID, date, endVZA, GF0_360) %>%
  rename(gap_fraction = GF0_360)

# Merge all data
model_data <- lt_results %>%
  left_join(lai_all_vza, by = c("plot" = "ID", "date")) %>%
  left_join(gap_all_vza, by = c("plot" = "ID", "date", "endVZA")) %>%
  filter(!is.na(LTLAI) & !is.na(gap_fraction) & !is.na(Le))

# Calculate LAI estimates using different clumping indices
model_data <- model_data %>%
  mutate(
    # Calculate different LAI estimates
    LAI_Le = Le,  # Effective LAI (no clumping correction)
    LAI_LX = Le / LX,
    LAI_LXG1 = Le / LXG1,  
    LAI_LXG2 = Le / LXG2
  )

# ===== WOODY COMPONENT REMOVAL FOR ALL LAI ESTIMATIONS =====
# Create wood-corrected versions for all LAI methods using 2024-12-16 baseline
model_data <- model_data %>%
  group_by(plot, endVZA) %>%
  mutate(
    # Find winter baselines for each LAI method (2024-12-16)
    winter_baseline_Le = ifelse(any(date == as.Date("2024-12-16")), 
                                LAI_Le[date == as.Date("2024-12-16")][1], NA),
    winter_baseline_LX = ifelse(any(date == as.Date("2024-12-16")), 
                                LAI_LX[date == as.Date("2024-12-16")][1], NA),
    winter_baseline_LXG1 = ifelse(any(date == as.Date("2024-12-16")), 
                                  LAI_LXG1[date == as.Date("2024-12-16")][1], NA),
    winter_baseline_LXG2 = ifelse(any(date == as.Date("2024-12-16")), 
                                  LAI_LXG2[date == as.Date("2024-12-16")][1], NA),
    
    # Create wood-corrected LAI estimates
    LAI_Le_nowood = pmax(LAI_Le - winter_baseline_Le, 0, na.rm = TRUE),
    LAI_LX_nowood = pmax(LAI_LX - winter_baseline_LX, 0, na.rm = TRUE),
    LAI_LXG1_nowood = pmax(LAI_LXG1 - winter_baseline_LXG1, 0, na.rm = TRUE),
    LAI_LXG2_nowood = pmax(LAI_LXG2 - winter_baseline_LXG2, 0, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    # Calculate errors for wood-corrected methods
    error_Le_nowood = LTLAI - LAI_Le_nowood,
    error_LX_nowood = LTLAI - LAI_LX_nowood,
    error_LXG1_nowood = LTLAI - LAI_LXG1_nowood,
    error_LXG2_nowood = LTLAI - LAI_LXG2_nowood,
    
    # Calculate absolute errors for wood-corrected methods
    abs_error_Le_nowood = abs(error_Le_nowood),
    abs_error_LX_nowood = abs(error_LX_nowood),
    abs_error_LXG1_nowood = abs(error_LXG1_nowood),
    abs_error_LXG2_nowood = abs(error_LXG2_nowood),
    
    # Also calculate errors for original methods
    error_Le = LTLAI - LAI_Le,
    error_LX = LTLAI - LAI_LX,
    error_LXG1 = LTLAI - LAI_LXG1,
    error_LXG2 = LTLAI - LAI_LXG2,
    
    # Calculate absolute errors for original methods
    abs_error_Le = abs(error_Le),
    abs_error_LX = abs(error_LX),
    abs_error_LXG1 = abs(error_LXG1),
    abs_error_LXG2 = abs(error_LXG2)
  )

cat("Added wood-corrected LAI columns for all methods using 2024-12-16 baseline\n")
cat("Dataset now contains", ncol(model_data), "columns\n")

# Calculate statistics by VZA for each method - CORRECTED VERSION
calculate_stats <- function(data, lai_col, error_col, abs_error_col, method_name) {
  data %>%
    group_by(endVZA) %>%
    do({
      current_data <- .
      lai_values <- current_data[[lai_col]]
      lt_values <- current_data$LTLAI
      
      # Calculate slope - CORRECTED: Reference (LTLAI) on Y-axis
      model <- lm(lai_values ~ lt_values)  # LTLAI ~ estimated_LAI
      slope_val <- coef(model)[2]
      
      # Calculate other stats
      data.frame(
        mean_abs_error = mean(current_data[[abs_error_col]], na.rm = TRUE),
        r_squared = cor(lt_values, lai_values, use = "complete.obs")^2,
        slope = slope_val,
        method = method_name
      )
    }) %>%
    ungroup()
}

# Calculate stats for wood-corrected methods
stats_Le_nowood <- calculate_stats(model_data, "LAI_Le_nowood", "error_Le_nowood", "abs_error_Le_nowood", "Effective LAI (no wood)")
stats_LX_nowood <- calculate_stats(model_data, "LAI_LX_nowood", "error_LX_nowood", "abs_error_LX_nowood", "Clumping LX (no wood)")
stats_LXG1_nowood <- calculate_stats(model_data, "LAI_LXG1_nowood", "error_LXG1_nowood", "abs_error_LXG1_nowood", "Clumping LXG1 (no wood)")
stats_LXG2_nowood <- calculate_stats(model_data, "LAI_LXG2_nowood", "error_LXG2_nowood", "abs_error_LXG2_nowood", "Clumping LXG2 (no wood)")

# Calculate stats for original methods (for comparison)
stats_Le <- calculate_stats(model_data, "LAI_Le", "error_Le", "abs_error_Le", "Effective LAI")
stats_LX <- calculate_stats(model_data, "LAI_LX", "error_LX", "abs_error_LX", "Clumping LX")
stats_LXG1 <- calculate_stats(model_data, "LAI_LXG1", "error_LXG1", "abs_error_LXG1", "Clumping LXG1")
stats_LXG2 <- calculate_stats(model_data, "LAI_LXG2", "error_LXG2", "abs_error_LXG2", "Clumping LXG2")

# Combine wood-corrected stats for main plot
all_stats_nowood <- bind_rows(stats_Le_nowood, stats_LX_nowood, stats_LXG1_nowood, stats_LXG2_nowood)

# Combine all stats for comparison
all_stats_comparison <- bind_rows(
  stats_Le, stats_LX, stats_LXG1, stats_LXG2,
  stats_Le_nowood, stats_LX_nowood, stats_LXG1_nowood, stats_LXG2_nowood
)

# Fix 1: Check what method names are actually in your data
cat("Method names in data:", unique(all_stats_nowood$method), "\n")

# Hardcode colors to specific methods - NO MORE CONFUSION!
p1 <- ggplot(all_stats_nowood, aes(x = endVZA, y = mean_abs_error, color = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Effective LAI (no wood)" = "#440154",     # Dark purple
               "Clumping LX (no wood)" = "#31688E",       # Blue
               "Clumping LXG1 (no wood)" = "#35B779",     # Green  
               "Clumping LXG2 (no wood)" = "#FDE725"),    # Yellow
    labels = c("Effective LAI (no wood)" = "Effective LAI",
               "Clumping LX (no wood)" = "Clumping LX", 
               "Clumping LXG1 (no wood)" = "Clumping LXG1",
               "Clumping LXG2 (no wood)" = "Clumping LXG2")
  ) +
  scale_x_continuous(breaks = seq(10, 90, 10)) +
  labs(x = "View Zenith Angle (VZA°)", y = "Absolute Error", color = "") +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 25),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 25),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 25, hjust = 0.5)
  )

p2 <- ggplot(all_stats_nowood, aes(x = endVZA, y = slope, color = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Effective LAI (no wood)" = "#440154",     # Dark purple
               "Clumping LX (no wood)" = "#31688E",       # Blue
               "Clumping LXG1 (no wood)" = "#35B779",     # Green  
               "Clumping LXG2 (no wood)" = "#FDE725"),    # Yellow
    labels = c("Effective LAI (no wood)" = "Effective LAI",
               "Clumping LX (no wood)" = "Clumping LX", 
               "Clumping LXG1 (no wood)" = "Clumping LXG1",
               "Clumping LXG2 (no wood)" = "Clumping LXG2")
  ) +
  scale_x_continuous(breaks = seq(10, 90, 10)) +
  labs(x = "View Zenith Angle (VZA°)", y = "Slope", color = "") +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 25),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 25),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 25, hjust = 0.5)
  )

p3 <- ggplot(all_stats_nowood, aes(x = endVZA, y = r_squared, color = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Effective LAI (no wood)" = "#440154",     # Dark purple
               "Clumping LX (no wood)" = "#31688E",       # Blue
               "Clumping LXG1 (no wood)" = "#35B779",     # Green  
               "Clumping LXG2 (no wood)" = "#FDE725"),    # Yellow
    labels = c("Effective LAI (no wood)" = "Effective LAI",
               "Clumping LX (no wood)" = "Clumping LX", 
               "Clumping LXG1 (no wood)" = "Clumping LXG1",
               "Clumping LXG2 (no wood)" = "Clumping LXG2")
  ) +
  scale_x_continuous(breaks = seq(10, 90, 10)) +
  scale_y_continuous(limits = c(0.6, 1.0)) +
  labs(x = "View Zenith Angle (VZA°)", y = "R²", color = "") +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 25),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 25),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    plot.title = element_text(size = 25, hjust = 0.5)
  )

# Create legend with hardcoded colors
p_with_legend <- ggplot(all_stats_nowood, aes(x = endVZA, y = mean_abs_error, color = method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("Effective LAI (no wood)" = "#440154",     # Dark purple
               "Clumping LX (no wood)" = "#31688E",       # Blue
               "Clumping LXG1 (no wood)" = "#35B779",     # Green  
               "Clumping LXG2 (no wood)" = "#FDE725"),    # Yellow
    labels = c("Effective LAI (no wood)" = "Effective LAI",
               "Clumping LX (no wood)" = "Clumping LX", 
               "Clumping LXG1 (no wood)" = "Clumping LXG1",
               "Clumping LXG2 (no wood)" = "Clumping LXG2")
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 25),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25),
    legend.margin = margin(t = 30, b = 15),
    legend.spacing.x = unit(1.5, "cm"),
    legend.key.width = unit(2, "cm")
  ) +
  guides(color = guide_legend(
    title = "",
    override.aes = list(size = 5),
    nrow = 1
  ))

# Extract legend using cowplot
legend <- get_legend(p_with_legend)

# Combine plots
plots_combined <- arrangeGrob(p1, p2, p3, ncol = 3, nrow = 1)

# Final combination
final_plot <- arrangeGrob(
  plots_combined,
  legend,
  ncol = 1, 
  nrow = 2,
  heights = c(0.8, 0.2)                                      # ← MORE SPACE FOR LEGEND
)

# Save the plot
ggsave(file.path(plot_dir, "vza_analysis_clumping_comparison_nowood.png"), final_plot, 
       width = 20, height = 7, dpi = 300)                    # ← BIGGER OVERALL SIZE

ggsave(file.path(plot_dir, "vza_analysis_clumping_comparison_nowood.pdf"), final_plot, 
       width = 20, height = 7)

# Save comparison data
write.csv(all_stats_comparison, file.path(plot_dir, "vza_stats_comparison_original_vs_nowood.csv"), row.names = FALSE)

cat("VZA analysis plot (wood-corrected) saved to:", plot_dir, "\n")
cat("Comparison statistics saved to:", plot_dir, "\n")

# Print summary statistics
cat("\n=== WOOD-CORRECTED SUMMARY STATISTICS ===\n")
print(all_stats_nowood)

cat("\n=== COMPARISON: Original vs Wood-Corrected ===\n")
print(all_stats_comparison)