#This is the script where we import the functions for image processing and apply it on our hemispherical images. Furthermore we import the resuts from the leaf area python script and combine both results into one csv file

library(dplyr)
library(parallel)

# Set working directory
setwd("/mnt/gsdata/users/lotz/LittervsLens/")

# Set number of cores to use
num_cores <- 4  # Change this based on your system capabilities

source("scripts/data_processing/functions.R") # Load helper functions

# Define paths - UPDATE THESE TO YOUR ACTUAL PATHS
input_folder <- "data/all/"
cropped_folder <- "data/hemi_photo_cropped/" 
csv_path <- "results/LAI_results.csv"

# Uncomment the following block to preview canopy parameters from a sample image
# # Step 1: Preview canopy parameters (run once to understand output structure)
# sample_images <- list.files(input_folder, pattern = "\\.(jpg|JPG)$", full.names = TRUE, recursive = TRUE)
# if (length(sample_images) > 0) {
#   cat("=== PREVIEW OF CANOPY PARAMETERS ===\n")
#   sample_canopy <- preview_canopy_params(sample_images[1])
# }

#Step 2: Crop images (if needed) - now with parallel processing
cat("\n=== CROPPING IMAGES ===\n")
crop_images(input_folder, cropped_folder, num_cores)


# Step 3: Get list of cropped images
image_files <- list.files(cropped_folder, pattern = "\\.(jpg|jpeg|JPG|JPEG)$", 
                         full.names = TRUE, recursive = TRUE)

cat("Found", length(image_files), "images to process\n")
cat("Using", num_cores, "parallel cores\n")

# Delete existing CSV to start fresh
if (file.exists(csv_path)) {
  file.remove(csv_path)
  cat("Removed existing CSV file\n")
}

# Step 4: Process all images with parallel processing
all_results_list <- mclapply(image_files, function(img) {
  tryCatch({
    result <- process_image_multi_vza(img, csv_path)
    return(list(success = TRUE, rows = ifelse(is.null(result), 0, nrow(result))))
  }, error = function(e) {
    cat("Error:", basename(img), "\n")
    return(list(success = FALSE))
  })
}, mc.cores = num_cores)

# Step 5: Final summary
successful_files <- sum(sapply(all_results_list, function(x) x$success))
failed_files <- length(image_files) - successful_files

cat("\nProcessing complete:\n")
cat("Success:", successful_files, "/ Failed:", failed_files, "\n")

if (file.exists(csv_path)) {
  final_data <- read.csv(csv_path, stringsAsFactors = FALSE)
  cat("Total rows in CSV:", nrow(final_data), "\n")
  cat("Results saved to:", csv_path, "\n")
} else {
  cat("No results file created\n")
}