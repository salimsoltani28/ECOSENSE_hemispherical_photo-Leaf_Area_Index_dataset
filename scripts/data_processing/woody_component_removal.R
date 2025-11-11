library(dplyr)
library(readr)

# Set working directory
setwd("/mnt/gsdata/users/lotz/LittervsLens/")

# File paths
lai_results_path <- "results/LAI_results.csv"

# Read the data
lai_results <- read_csv(lai_results_path)

# Or if it's already character/numeric
lai_results$date <- as.Date(as.character(lai_results$date), format = "%Y%m%d")

lai_results$LAI_LXG1 <- lai_results$Le / lai_results$LXG1

# ===== WOODY COMPONENT REMOVAL =====
# Step 1: Extract baseline values (2024-12-16) for each ID+endVZA combination
baseline_values <- lai_results %>%
  filter(date == as.Date("2024-12-16")) %>%
  select(ID, endVZA, baseline_LAI = LAI_LXG1)

# Step 2: Join baseline with all data and subtract
lai_corrected <- lai_results %>%
  left_join(baseline_values, by = c("ID", "endVZA")) %>%
  mutate(
    # Subtract baseline from all measurements (wood-corrected LAI)
    LAI_LXG1_nowood = LAI_LXG1 - baseline_LAI,
    # Ensure no negative values
    LAI_LXG1_nowood = pmax(LAI_LXG1_nowood, 0, na.rm = TRUE)
  ) %>%
  select(ID, date, endVZA, LAI_LXG1, baseline_LAI, LAI_LXG1_nowood)

# Preview the result
cat("Results preview:\n")
print(head(lai_corrected, 10))

# Optionally save to CSV
write_csv(lai_corrected, "results/LAI_results_wood_removal.csv")

# Read both files
lai_results_wood_removal <- read_csv("results/LAI_results_wood_removal.csv")
lt_results <- read_csv("results/LT_results.csv")

# Check the structure
cat("LAI_results_wood_removal columns:", paste(colnames(lai_results_wood_removal), collapse = ", "), "\n")
cat("LT_results columns:", paste(colnames(lt_results), collapse = ", "), "\n")

# Since LTLAI is the same for each plot-date combination regardless of endVZA,
# we only need one LTLAI value per plot-date combination
lt_results_unique <- lt_results %>%
  select(plot, date, LTLAI) %>%
  distinct()  # Remove duplicates since LTLAI is same for all endVZA

# Join the datasets
# ID in lai_results_wood_removal corresponds to plot in lt_results
combined_data <- lai_results_wood_removal %>%
  left_join(lt_results_unique, by = c("ID" = "plot", "date" = "date"))

# Check the result
cat("Combined data preview:\n")
print(head(combined_data, 10))

cat("Combined data dimensions:", nrow(combined_data), "rows,", ncol(combined_data), "columns\n")

# Check if all records got LTLAI values
cat("Records without LTLAI:", sum(is.na(combined_data$LTLAI)), "\n")

# Save the combined dataset
write_csv(combined_data, "/mnt/gsdata/users/lotz/LittervsLens/results/LAI_results_wood_removal_full.csv")

cat("Combined dataset saved!\n")