#This Script contains helper functions to process fisheye images for LAI analysis using hemispheR package

# Set working directory
setwd("/mnt/gsdata/users/lotz/LittervsLens/")

# Set CRAN mirror to avoid selection prompt
options(repos = c(CRAN = "https://cloud.r-project.org"))

# First check if required packages are installed. If not they get installed.
required_packages <- c("hemispheR", "stringr", "dplyr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Function to crop all images from input_dir into output_dir
##### This is done so we only analyze the image itself and set the Zenith angles right. Before the images had some small black borders at the top and bottom
# Parallel version of crop_images function
crop_images <- function(input_dir, cropped_dir, num_cores = 10) {
  if (!require("imager")) install.packages("imager")
  library(imager)
  
  if (!dir.exists(cropped_dir)) dir.create(cropped_dir, recursive = TRUE)
  
  image_files <- list.files(input_dir, pattern = "\\.(jpg|jpeg|png|tif|bmp|JPG|JPEG|PNG|TIF|BMP)$", 
                            full.names = TRUE, recursive = TRUE)
  
  cat("Cropping", length(image_files), "images using", num_cores, "cores\n")
  
  # Process images in parallel
  mclapply(image_files, function(image_path) {
    tryCatch({
      image <- load.image(image_path)
      cropped_image <- crop.borders(image, ny = c(220, 280)) # set this to the size of the border to crop
      
      relative_path <- sub(input_dir, "", image_path)
      output_path <- file.path(cropped_dir, dirname(relative_path))
      if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
      
      output_file <- file.path(output_path, sub("\\.JPG$", ".jpg", basename(image_path)))
      
      save.image(cropped_image, output_file)
      cat("Processed:", basename(image_path), "\n")
    }, error = function(e) {
      cat("Failed to process:", basename(image_path), "Error:", e$message, "\n")
    })
  }, mc.cores = num_cores)
}

# Helper to extract date from filename, the naming needs to be as this: LT11_20240920.JPG
extract_date <- function(filename) {
  str_extract(filename, "\\d{8}")  # expects date like YYYYMMDD
}

# Helper to extract plot/ID from filename
extract_plot <- function(filename) {
  str_remove(filename, "_\\d{8}.*")  # remove the date and extension
}

# Function to preview canopy parameters from a sample image this helps for debugging and understanding the output structure of canopy_fisheye
preview_canopy_params <- function(sample_image_path) {
  cat("Analyzing sample image to understand canopy_fisheye output structure...\n")
  
  img <- import_fisheye(sample_image_path,
                        channel = 3,
                        circular = TRUE,
                        gamma = 2.2,
                        stretch = FALSE,
                        display = FALSE,
                        message = FALSE)
  
  img.bw <- binarize_fisheye(img,
                             method = 'Otsu',
                             zonal = FALSE,
                             display = FALSE,
                             export = FALSE)
  
  gap.frac <- gapfrac_fisheye(img.bw,
                              maxVZA = 90,
                              lens = "Sigma-4.5",
                              startVZA = 0,
                              endVZA = 20,
                              nrings = 5,
                              nseg = 8,
                              display = FALSE,
                              message = FALSE)
  
  canopy <- canopy_fisheye(gap.frac)
  cat("---\n")
  cat("Canopy parameters available:\n")
  if (is.list(canopy)) {
    cat("Parameter names:", paste(names(canopy), collapse = ", "), "\n")
    cat("Data structure:\n")
    str(canopy)
  } else {
    cat("Canopy object structure:\n")
    str(canopy)
  }
  
  return(canopy)
}
# From the list you can check which parameters you want to extract and include in the final CSV


# Simplified function that processes image and writes directly to CSV
process_image_multi_vza <- function(image_path, csv_path) {
  file_name <- basename(image_path)
  date <- extract_date(file_name)
  plot <- extract_plot(file_name)
  
  # Handle cases where extraction fails
  if (is.null(date) || is.na(date) || length(date) == 0) {
    cat("Warning: Could not extract date from", file_name, "\n")
    return(NULL)
  }
  if (is.null(plot) || is.na(plot) || length(plot) == 0) {
    cat("Warning: Could not extract plot from", file_name, "\n")
    return(NULL)
  }
  
  cat("Processing", file_name, "(Plot:", plot, ", Date:", date, ")\n")
  
  # Import and binarize image once
  tryCatch({
    img <- import_fisheye(image_path,
                          channel = 3,
                          circular = TRUE,
                          gamma = 2.2,
                          stretch = FALSE,
                          display = FALSE,
                          message = FALSE)
    
    img.bw <- binarize_fisheye(img,
                               method = 'Otsu',
                               zonal = FALSE,
                               display = FALSE,
                               export = FALSE)
  }, error = function(e) {
    cat("Error importing/binarizing", file_name, ":", e$message, "\n")
    return(NULL)
  })
  
  # Test different endVZA values from 10 to 90 in 10° steps
  vza_ranges <- seq(10, 90, by = 10)
  all_results <- data.frame()
  
  for (endVZA in vza_ranges) {
    tryCatch({
      gap.frac <- gapfrac_fisheye(img.bw,
                                  maxVZA = 90,
                                  lens = "Sigma-4.5",
                                  startVZA = 0,
                                  endVZA = endVZA,
                                  nrings = 5,
                                  nseg = 8,
                                  display = FALSE,
                                  message = FALSE)
      
      canopy <- canopy_fisheye(gap.frac)
      
      # Extract what we need from the results
      result_row <- data.frame(
        ID = plot,
        date = date,
        endVZA = endVZA,
        filename = file_name,
        Le = as.numeric(canopy$Le[1]),
        L = as.numeric(canopy$L[1]),
        LX = as.numeric(canopy$LX[1]),
        LXG1 = as.numeric(canopy$LXG1[1]),
        LXG2 = as.numeric(canopy$LXG2[1]),
        DIFN = as.numeric(canopy$DIFN[1]),
        MTA_ell = as.numeric(canopy$MTA.ell[1]),
        stringsAsFactors = FALSE
      )
      
      all_results <- rbind(all_results, result_row)
      
    }, error = function(e) {
      cat("Error processing", file_name, "with endVZA =", endVZA, ":", e$message, "\n")
    })
  }
  
  # Write directly to CSV (append mode)
  if (nrow(all_results) > 0) {
    # Check if CSV exists and has header
    if (file.exists(csv_path)) {
      # Append without header
      write.table(all_results, csv_path, 
                  sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    } else {
      # Create new with header
      write.csv(all_results, csv_path, row.names = FALSE)
    }
    cat("Wrote", nrow(all_results), "rows to CSV for", file_name, "\n")
  }
  
  return(all_results)
}

# Keep original function for backward compatibility
process_image <- function(image_path, endVZA) {
  file_name <- basename(image_path)
  date <- extract_date(file_name)
  plot <- extract_plot(file_name)

  img <- import_fisheye(image_path,
                        channel = 3, #only use 1 if you increased the contrast with 2BG; DEFAULT: 3 we run the analysis on the blue channel
                        circular = TRUE,
                        gamma = 2.2,
                        stretch = FALSE,
                        display = FALSE,
                        message = TRUE)

  img.bw <- binarize_fisheye(img,
                             method = 'Otsu',
                             zonal = FALSE,
                             manual = NULL,
                             display = FALSE,
                             export = FALSE)

  gap.frac <- gapfrac_fisheye(img.bw,
                              maxVZA = 90,
                              lens = "Sigma-4.5",
                              startVZA = 0,
                              endVZA = endVZA,
                              nrings = 5,
                              nseg = 8,
                              display = FALSE,
                              message = FALSE)

  canopy <- canopy_fisheye(gap.frac)

  results <- data.frame(
    ID = as.character(plot),
    date = as.character(date),
    LAI = as.numeric(canopy$Le / canopy$LXG2),
    lat = as.numeric(NA),
    lon = as.numeric(NA),
    stringsAsFactors = FALSE
  )

  return(results)
}

# Updated CSV function to handle expanded data structure
update_csv_multi <- function(new_data, csv_path) {
  # Ensure consistent data types for key columns
  new_data$ID <- as.character(new_data$ID)
  new_data$date <- as.character(new_data$date)
  new_data$endVZA <- as.numeric(new_data$endVZA)
  new_data$filename <- as.character(new_data$filename)
  
  # Load existing data if it exists
  if (file.exists(csv_path)) {
    existing <- read.csv(csv_path, stringsAsFactors = FALSE)
  } else {
    existing <- data.frame()
  }
  
  # Combine and write
  combined <- dplyr::bind_rows(existing, new_data)
  write.csv(combined, csv_path, row.names = FALSE)
  
  cat("Updated CSV with", nrow(new_data), "new rows. Total rows:", nrow(combined), "\n")
}

### Original update_csv function for backward compatibility
update_csv <- function(new_data, csv_path) {
  # Define expected structure
  enforce_types <- function(df) {
    df$ID   <- as.character(df$ID)
    df$date <- as.character(df$date)
    df$LAI  <- as.numeric(df$LAI_LXG1)
    df
  }

  # Apply structure to new data
  new_data <- new_data[, c("ID", "date", "LAI_LXG1", "lat", "lon"), drop = FALSE]
  new_data <- enforce_types(new_data)

  # Load or initialize existing data
  if (file.exists(csv_path)) {
    existing <- read.csv(csv_path, stringsAsFactors = FALSE)
    existing <- existing[, c("ID", "date", "LAI_LXG1", "lat", "lon"), drop = FALSE]
    existing <- enforce_types(existing)
  } else {
    existing <- data.frame(
      ID = character(),
      date = character(),
      LAI_LXG1 = numeric(),
      stringsAsFactors = FALSE
    )
  }

  # Combine and write
  combined <- dplyr::bind_rows(existing, new_data)
  write.csv(combined, csv_path, row.names = FALSE)
}