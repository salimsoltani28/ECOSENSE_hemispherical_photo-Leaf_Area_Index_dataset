# ============================================================
#   F05: CANOPY STRUCTURE ANALYSIS — CORRECTED FOR hemispheR G(theta)
# ============================================================

# Set working directory
setwd("/mnt/gsdata/users/lotz/LittervsLens/")

library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(grid)
library(purrr)
library(viridis)

# --- Load datasets ---
lai_csv_path <- "results/LAI_results.csv"
gap_csv_path <- "results/LAI_results_gapfrac.csv"

lai_results <- read.csv(lai_csv_path, stringsAsFactors = FALSE)
gap_results <- read.csv(gap_csv_path, stringsAsFactors = FALSE)

# --- Output directory ---
plot_dir <- "plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# --- Select plots ---
selected_ids <- c('LT11', 'LT13', 'LT14', 'LT23', 'LT24', 'LT33', 'LT34', 
                  'LT41', 'LT44', 'LT51', 'LT52', 'LT53', 'LT62', 'LT61', 'LT63')

filtered_lai_data <- lai_results %>% filter(ID %in% selected_ids)
filtered_gap_data <- gap_results %>% filter(ID %in% selected_ids)

# --- Colorblind-friendly palette ---
colors <- c("#31688E", "#FDE725", "#35B779")

# ============================================================
# 1️⃣ GAP FRACTION
# ============================================================

gap_fraction_summary <- filtered_gap_data %>%
  group_by(endVZA) %>%
  summarise(
    gap_mean = mean(GF0_360, na.rm = TRUE),
    gap_se = sd(GF0_360, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p_gap_fraction <- ggplot(gap_fraction_summary, aes(x = endVZA, y = gap_mean)) +
  geom_line(linewidth = 1.2, color = colors[1]) +
  geom_point(size = 3, color = colors[1]) +
  geom_errorbar(aes(ymin = gap_mean - gap_se, ymax = gap_mean + gap_se),
                width = 2, color = colors[1], alpha = 0.7, linewidth = 0.8) +
  labs(x = "View Zenith Angle (°)", y = "Gap Fraction") +
  theme_minimal(base_size = 25) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.title = element_text(size = 25, face = "bold"),
    axis.text = element_text(size = 25),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(10, 90, 20)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1))

# ============================================================
# 2️⃣ CLUMPING INDICES (LX, LXG1, LXG2)
# ============================================================

clumping_summary <- filtered_lai_data %>%
  group_by(endVZA) %>%
  summarise(
    mean_LX = mean(LX, na.rm = TRUE),
    se_LX = sd(LX, na.rm = TRUE) / sqrt(n()),
    mean_LXG1 = mean(LXG1, na.rm = TRUE),
    se_LXG1 = sd(LXG1, na.rm = TRUE) / sqrt(n()),
    mean_LXG2 = mean(LXG2, na.rm = TRUE),
    se_LXG2 = sd(LXG2, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

clumping_long <- clumping_summary %>%
  pivot_longer(
    cols = c(mean_LX, mean_LXG1, mean_LXG2),
    names_to = "clumping_type",
    values_to = "mean_value",
    names_prefix = "mean_"
  ) %>%
  mutate(
    se_value = case_when(
      clumping_type == "LX" ~ se_LX,
      clumping_type == "LXG1" ~ se_LXG1,
      clumping_type == "LXG2" ~ se_LXG2
    )
  )

p_clumping <- ggplot(clumping_long, aes(x = endVZA, y = mean_value, color = clumping_type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_value - se_value, ymax = mean_value + se_value),
                width = 1.5, alpha = 0.7, linewidth = 0.8) +
  scale_color_manual(
    name = "Clumping Index",
    values = colors,
    labels = c("LX" = "LX", "LXG1" = "LXG1", "LXG2" = "LXG2")
  ) +
  labs(x = "View Zenith Angle (°)", y = "Clumping Index") +
  theme_minimal(base_size = 35) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.title = element_text(size = 25, face = "bold"),
    axis.text = element_text(size = 25),
    legend.position = "top",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(10, 90, 20))

# ============================================================
# 3️⃣ G(θ) — CORRECT hemispheR version
# ============================================================

# Campbell (1986) ellipsoidal LAD model:
# G(θ) = sqrt(cos²(θ) + sin²(θ) / tan²(χ))
# χ = MTA (mean tilt angle) in radians

gtheta_summary <- filtered_lai_data %>%
  filter(!is.na(MTA_ell)) %>%
  mutate(
    theta_rad = endVZA * pi / 180,
    chi_rad = MTA_ell * pi / 180,
    G_theta = sqrt(cos(theta_rad)^2 + (sin(theta_rad)^2) / (tan(chi_rad)^2))
  ) %>%
  group_by(endVZA) %>%
  summarise(
    gtheta = mean(G_theta, na.rm = TRUE),
    sd_gtheta = sd(G_theta, na.rm = TRUE),
    n = sum(!is.na(G_theta)),
    se_gtheta = ifelse(n > 1, sd_gtheta / sqrt(n), 0),
    .groups = "drop"
  )

p_gtheta <- ggplot(gtheta_summary, aes(x = endVZA, y = gtheta)) +
  geom_line(linewidth = 1.2, color = colors[2]) +
  geom_point(size = 3, color = colors[2]) +
  geom_errorbar(aes(ymin = pmax(gtheta - se_gtheta, 0),
                    ymax = pmin(gtheta + se_gtheta, 1.5)), 
                width = 2, color = colors[2], alpha = 0.7, linewidth = 0.8) +
  labs(x = "View Zenith Angle (°)",
       y = "G(theta)") +
  theme_minimal(base_size = 25) +
  theme(
    text = element_text(family = "Helvetica"),
    axis.title = element_text(size = 25, face = "bold"),
    axis.text = element_text(size = 25),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(breaks = seq(10, 90, 20)) +
  scale_y_continuous(limits = c(0, 2))

# ============================================================
# 4️⃣ SAVE PLOTS
# ============================================================

# ggsave(file.path(plot_dir, "gap_fraction_vs_vza.png"), p_gap_fraction, width = 8, height = 6, dpi = 300)
# ggsave(file.path(plot_dir, "clumping_indices_vs_vza.png"), p_clumping, width = 8, height = 6, dpi = 300)
# ggsave(file.path(plot_dir, "gtheta_vs_vza.png"), p_gtheta, width = 8, height = 6, dpi = 300)

# # --- Combined plot ---
# combined_plot <- grid.arrange(p_gap_fraction, p_clumping, p_gtheta, ncol = 3, nrow = 1)

# ggsave(file.path(plot_dir, "combined_canopy_parameters_vs_vza.png"), combined_plot,
#        width = 24, height = 8, dpi = 300)
# ggsave(file.path(plot_dir, "combined_canopy_parameters_vs_vza.pdf"), combined_plot,
#        width = 24, height = 8)

#        # --- Combined plot with better spacing and legend ---
# plots_combined <- grid.arrange(p_gap_fraction, p_clumping, p_gtheta, ncol = 3, nrow = 1)

# Recreate plots combination without legend
plots_combined <- grid.arrange(p_gap_fraction, p_clumping, p_gtheta, ncol = 3, nrow = 1)

# Final combination with more space allocation
final_plot <- arrangeGrob(
  plots_combined,
  ncol = 1, 
  nrow = 2,
  heights = c(0.8, 0.2)  # More space for legend
)

# Save with better dimensions
ggsave(file.path(plot_dir, "combined_canopy_parameters_vs_vza.png"), final_plot,
       width = 20, height = 7, dpi = 300)  # Bigger overall size
ggsave(file.path(plot_dir, "combined_canopy_parameters_vs_vza.pdf"), final_plot,
       width = 20, height = 7)

# ============================================================
# 5️⃣ SUMMARY OUTPUT
# ============================================================

cat("=== SUMMARY STATISTICS ===\n\n")
cat("--- Gap Fraction ---\n")
print(gap_fraction_summary)

cat("\n--- Clumping Indices ---\n")
print(clumping_summary)

cat("\n--- G(theta) ---\n")
print(gtheta_summary)

write.csv(gap_fraction_summary, file.path(plot_dir, "gap_fraction_summary_statistics.csv"), row.names = FALSE)
write.csv(clumping_summary, file.path(plot_dir, "clumping_summary_statistics.csv"), row.names = FALSE)
write.csv(gtheta_summary, file.path(plot_dir, "gtheta_summary_statistics.csv"), row.names = FALSE)

cat("\nPlots and summary statistics saved to:", plot_dir, "\n")
