# F02 and F06 Code: LAI over time with segmented regression and hemispherical photos timeline

library(ggplot2)
library(dplyr)
library(segmented)
library(purrr)

# Set working directory
setwd("/mnt/gsdata/users/lotz/LittervsLens/")

data_raw <- read.csv("results/combined_LAI_wood_removal_with_LTLAI.csv")

str(data_raw)
selected_ids <- c('LT11', 'LT13', 'LT14', 'LT23', 'LT24', 'LT33', 'LT34', 'LT41', 'LT44', 'LT51', 'LT52', 'LT53', 'LT62', 'LT61', 'LT63')

# Filter data for LAI.LXG1 at 20° only
filtered_data <- data_raw %>%
  filter(endVZA == 20) %>%
  filter(ID %in% selected_ids)

# Calculate summary statistics for both actual LAI and LAI.LXG1 (20°)
summary_data <- filtered_data %>%
  group_by(date) %>%
  summarise(
    mean_LT_LAI = mean(LTLAI, na.rm = TRUE),
    p5_LT_LAI = quantile(LTLAI, 0.05, na.rm = TRUE),
    p95_LT_LAI = quantile(LTLAI, 0.95, na.rm = TRUE),
    mean_DHP_LAI = mean(LAI_LXG1_nowood, na.rm = TRUE),
    p5_DHP_LAI = quantile(LAI_LXG1_nowood, 0.05, na.rm = TRUE),
    p95_DHP_LAI = quantile(LAI_LXG1_nowood, 0.95, na.rm = TRUE)
  )

# ===== FIX DATE FORMAT =====
# Convert date column to proper Date format
summary_data$date <- as.Date(summary_data$date)

# ===== CREATE SEGMENTED MODEL =====
# Convert dates to numeric ONLY for segmented regression
summary_data$date_numeric <- as.numeric(summary_data$date)

# Create initial linear model (required by segmented)
lm_model <- lm(mean_LT_LAI ~ date_numeric, data = summary_data)

# Create segmented model
seg_model <- segmented(lm_model, seg.Z = ~date_numeric, npsi = 2)

# Extract breakpoints
breakpoints_numeric <- seg_model$psi[, 2]
cat("Breakpoints (numeric):", breakpoints_numeric, "\n")

# Convert back to dates for plotting
breakpoints_date <- as.Date(breakpoints_numeric, origin = "1970-01-01")
cat("Breakpoints (dates):", as.character(breakpoints_date), "\n")

# Create the plot
p <- ggplot(summary_data, aes(x = date)) +
  # Shaded region for LT LAI 5th-95th percentile range
  geom_ribbon(aes(ymin = p5_LT_LAI, ymax = p95_LT_LAI), fill = "black", alpha = 0.2) +
  # Shaded region for DHP LAI 5th-95th percentile range
  geom_ribbon(aes(ymin = p5_DHP_LAI, ymax = p95_DHP_LAI), fill = "#31688E", alpha = 0.2) +
  # Line for mean LT LAI
  geom_line(aes(y = mean_LT_LAI, color = "Mean LT LAI"), linewidth = 1.2, linetype = "solid") +  # ← FIXED: size to linewidth
  # Line for mean DHP LAI
  geom_line(aes(y = mean_DHP_LAI, color = "Mean DHP LAI"), linewidth = 1, linetype = "solid") +   # ← FIXED: size to linewidth
  # Add breakpoints as dashed vertical lines
  geom_vline(xintercept = breakpoints_date, color = "black", linetype = "dashed") +
  labs(
    x = "Date", 
    y = "LAI",
    color = NULL
  ) +
  scale_color_manual(values = c("Mean LT LAI" = "black", "Mean DHP LAI" = "#31688E")) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 25),
    plot.title = element_text(size = 25, face = "bold"),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.title = element_blank(),
    legend.position = "bottom"
  )

# Save the plot
ggsave("results/LAI_over_time.png", plot = p, width = 10, height = 6, dpi = 300)

# ===== F02 RECONSTRUCTION: Hemispherical Photos Timeline =====
library(jpeg)
library(grid)
library(gridExtra)

# Define the dates for hemispherical photos (matching the timeline)
photo_dates <- c("2024-09-20", "2024-09-30", "2024-10-11", "2024-10-15", 
                 "2024-10-25", "2024-11-12", "2024-12-02", "2024-12-16")

# Convert to Date objects
photo_dates <- as.Date(photo_dates)

# Load hemispherical photos for f02 visualization
photo_files <- c(
  "/mnt/gsdata/users/lotz/Hemi_Photo/f02/LT14_20240920.jpg",
  "/mnt/gsdata/users/lotz/Hemi_Photo/f02/LT14_20240930.jpg",
  "/mnt/gsdata/users/lotz/Hemi_Photo/f02/LT14_20241011.jpg",
  "/mnt/gsdata/users/lotz/Hemi_Photo/f02/LT14_20241015.jpg",
  "/mnt/gsdata/users/lotz/Hemi_Photo/f02/LT14_20241025.jpg",
  "/mnt/gsdata/users/lotz/Hemi_Photo/f02/LT14_20241112.jpg",
  "/mnt/gsdata/users/lotz/Hemi_Photo/f02/LT14_20241202.jpg",
  "/mnt/gsdata/users/lotz/Hemi_Photo/f02/LT14_20241216.jpg"
)

# Function to create simple photo plots using annotation_raster
create_photo_plot <- function(img_path, date_label) {
  if (file.exists(img_path)) {
    img <- readJPEG(img_path)
    
    p <- ggplot() +
      annotation_raster(img, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
      xlim(0, 1) + ylim(0, 1) +
      coord_fixed(ratio = 1) +
      theme_void() +
      theme(
        plot.margin = margin(0, 0, 0, 0),
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    return(p)
  } else {
    # Return empty plot if file doesn't exist
    return(ggplot() + theme_void() + 
           annotate("text", x = 0.5, y = 0.5, label = "Image\nNot Found", size = 4))
  }
}

# Create individual photo plots
photo_plots <- map2(photo_files, format(photo_dates, "%d %b"), create_photo_plot)

# Create LAI data for LT14 specifically
lt14_data <- filtered_data %>%
  filter(ID == "LT14") %>%
  mutate(date = as.Date(date)) %>%  # Convert date to Date format
  arrange(date) %>%
  mutate(lai_difference = c(0, abs(diff(LTLAI))))  # LTLAI differences for LT14

# Debug: Check what dates we have for LT14
cat("Available dates for LT14:\n")
print(lt14_data$date)
cat("LTLAI values for LT14:\n")
print(lt14_data$LTLAI)
cat("LAI difference values for LT14:\n")
print(lt14_data$lai_difference)

# Check which photo dates have actual data
photo_dates_available <- photo_dates[photo_dates %in% lt14_data$date]
cat("Photo dates with available data:\n")
print(photo_dates_available)

# Create plot data with numeric positions for even spacing
plot_positions <- 1:length(photo_dates)
names(plot_positions) <- as.character(photo_dates)

# Map LT14 data to plot positions
lt14_plot_data <- data.frame(
  position = plot_positions[as.character(lt14_data$date)],
  LTLAI = lt14_data$LTLAI,
  lai_difference = lt14_data$lai_difference,
  date = lt14_data$date
) %>%
  filter(!is.na(position))

# Create the main plot with LTLAI curves for LT14
f02_plot <- ggplot(lt14_plot_data, aes(x = position)) +
  # LTLAI line (blue)
  geom_line(aes(y = LTLAI, color = "LAI"), linewidth = 2) +
  geom_point(aes(y = LTLAI, color = "LAI"), size = 3) +
  
  # LTLAI Difference line (orange)
  geom_line(aes(y = lai_difference, color = "LAI Difference"), linewidth = 2) +
  geom_point(aes(y = lai_difference, color = "LAI Difference"), size = 3) +
  
  # Add vertical lines for photo dates
  geom_vline(xintercept = plot_positions, color = "gray", linetype = "dotted", alpha = 0.7) +
  
  scale_color_manual(
    values = c("LAI" = "#31688E", "LAI Difference" = "#FDE725")
  ) +
  
  scale_x_continuous(
    breaks = 1:length(photo_dates),
    labels = format(photo_dates, "%d %B %Y"),
    expand = c(0.05, 0.05)
  ) +
  
  scale_y_continuous(
    limits = c(0, 6),
    breaks = c(0, 2, 4, 6)
  ) +
  
  labs(
    x = "Date",
    y = "LAI"
  ) +
  
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", size = 25),
    axis.title.x = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold"),
    axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.position = c(0.85, 0.8),
    legend.direction = "vertical",
    legend.background = element_rect(fill = "white", color = "grey"),
    legend.margin = margin(8, 8, 8, 8),
    legend.title = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )

# Arrange photos in a horizontal row
photos_row <- do.call(grid.arrange, c(photo_plots, list(ncol = 8, nrow = 1)))

# Combine photos and timeline plot with less spacing
f02_complete <- arrangeGrob(
  photos_row,
  f02_plot,
  ncol = 1,
  nrow = 2,
  heights = c(0.55, 0.45),  # Photos take 55% of height, timeline takes 45% (more squished)
  padding = unit(0.05, "line")  # Minimal spacing between photos and plot
)

# Save the complete f02 visualization
ggsave("results/f02_complete.png", plot = f02_complete, width = 20, height = 12, dpi = 300)
ggsave("results/f02_complete.pdf", plot = f02_complete, width = 20, height = 12)

# Also save individual components
ggsave("results/f02_timeline_plot.png", plot = f02_plot, width = 16, height = 8, dpi = 300)

cat("Complete F02 visualization with hemispherical photos saved to results/f02_complete.png\n")
