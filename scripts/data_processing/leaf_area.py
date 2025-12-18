#This script contains functions to process images of Litter Trap content for leaf area calculation.

# dependencies
import cv2
import numpy as np
import os
import csv
from scipy.interpolate import griddata

### FUNCTIONS ###

# Function to warp image with ArUco markers
def warp_image(image, width_mm, height_mm, file_in):
    """
    Warps an input image based on the positions of ArUco markers, transforming it into a
    specified size defined in millimeters (width and height).

    Parameters:
    -----------
        image (numpy array): The input image containing ArUco markers.
        width_mm (float): The desired width of the output image in millimeters.
        height_mm (float): The desired height of the output image in millimeters.
        filename (str): The name of the image file for logging purposes.

    Returns:
    --------
        numpy array: The warped image with the specified size, or None if not enough markers are detected.
    """
    aruco_dict = cv2.aruco.getPredefinedDictionary(cv2.aruco.DICT_4X4_50)
    parameters = cv2.aruco.DetectorParameters()
    detector = cv2.aruco.ArucoDetector(aruco_dict, parameters)

    # Detect ArUco markers in the image
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    corners, ids, _ = detector.detectMarkers(gray)

    # Define the expected marker IDs for the four corners
    expected_ids = [0, 1, 2, 3]  # Replace these with your actual expected IDs

    if ids is not None and len(ids) >= 4:
        # Calculate centroids of detected squares
        quadrilateral_centroids = [np.mean(quadrilateral[0], axis=0) for quadrilateral in corners]
        quadrilateral_centroids = np.array(quadrilateral_centroids)

        # Calculate overall centroid
        overall_centroid = np.mean(quadrilateral_centroids, axis=0)

        # Detect corners of large rectangle
        large_quad_corners = []
        for i, quadrilateral in enumerate(corners):
            points = quadrilateral[0]
            q_centroid = quadrilateral_centroids[i]
            direction_vector = q_centroid - overall_centroid
            direction_vector /= np.linalg.norm(direction_vector)  # Normalization
            # Project points on direction vector
            projections = np.dot(points - overall_centroid, direction_vector)
            max_index = np.argmax(projections)
            large_quad_corners.append(points[max_index])

        large_quad_corners = np.array(large_quad_corners)
        sorted_indices = np.argsort(ids.flatten())
        pts_src = large_quad_corners[sorted_indices]

        # Define target points for specified size (e.g., 210 mm x 297 mm)
        dpi = 72
        width, height = round(width_mm * dpi / 25.4), round(height_mm * dpi / 25.4)
        pts_dst = np.array([
            [0, 0],
            [width - 1, 0],
            [width - 1, height - 1],
            [0, height - 1]
        ], dtype="float32")

        # Perform perspective transformation
        M = cv2.getPerspectiveTransform(pts_src, pts_dst)
        warped = cv2.warpPerspective(image, M, (width, height))
        
        print('Warped successfully')
        return warped
    else:
        # Determine which markers are missing
        detected_ids = ids.flatten() if ids is not None else []
        missing_ids = set(expected_ids) - set(detected_ids)
        
        # Map marker IDs to corners
        corner_map = {0: "top-left", 1: "top-right", 2: "bottom-right", 3: "bottom-left"}
        
        # Print details when markers are missing
        if missing_ids:
            print(f"Missing markers for image '{file_in}':")
            for marker in missing_ids:
                corner = corner_map.get(marker, "Unknown corner")
                print(f" - Marker ID {marker} ({corner})")
        else:
            print(f"Unexpected error: no markers detected.")
        
        return None


# Variable image exposure correction function
def correct_image_exposure(image,
                           width_mm=980,
                           height_mm=1980,
                           square_size_mm=10, 
                           square_positions_fix=[],
                           square_positions_groups_brightest=[]):
    """
    Corrects the exposure of an image using interpolation based on reference points (squares) 
    across the image. The reference points are defined by their positions and their gray level
    values. The correction is applied by adjusting the gray values so that the selected squares 
    approach ideal white (255 in grayscale).

    Parameters:
    -----------
    image : numpy array (cv2.imread())
        Input image file (number of pixels height, number of pixels width, 3 channels BGR)
        
    width_mm : float, optional, default=980
        The physical width of the image in millimeters. Used to convert mm to pixels.
        
    height_mm : float, optional, default=1980
        The physical height of the image in millimeters. Used to convert mm to pixels.
        
    square_size_mm : float, optional, default=10
        The size of the reference squares in millimeters (1 cm x 1 cm by default). These squares 
        are used to determine correction values for exposure based on their average gray values.
        
    square_positions_fix : list of tuples, optional
        A list of fixed square center positions (in mm), where each tuple represents the x and y 
        coordinates in millimeters. These squares are directly used as reference points for the 
        interpolation.

    square_positions_groups_brightest : list of lists, optional
        A list of lists, where each inner list contains square center positions (in mm) that belong 
        to a group. For each group, only the brightest square (highest average gray value) will be 
        selected and used as a reference point for the interpolation. This is to avoid reference
        squares that are covered by leaves accidentially.

    Returns:
    --------
    corrected_image : numpy array
        Corrected image with exposure adjusted (number of pixels height, number of pixels width, 1 channel gray)
    """
    
    # Check if the image was loaded correctly
    if image is None:
        raise ValueError("The image could not be loaded.")
    
    # Convert physical dimensions from mm to pixels
    image_height, image_width = image.shape[:2]
    
    width_ratio = image_width / width_mm
    height_ratio = image_height / height_mm

    gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    
    # Convert mm into pixels for square size
    square_size_px = int(square_size_mm * width_ratio)
    
    # List of reference points (position and gray level difference from white)
    points = []
    values = []

    # Process the square_positions
    for (x_mm, y_mm) in square_positions_fix:
        # Convert the center positions from mm to pixel coordinates
        x_px = int(x_mm * width_ratio)
        y_px = int(y_mm * height_ratio)
        
        # Determine the area of the square around this center
        x_start = max(0, x_px - square_size_px // 2)
        y_start = max(0, y_px - square_size_px // 2)
        x_end = min(image_width, x_px + square_size_px // 2)
        y_end = min(image_height, y_px + square_size_px // 2)
        
        # Extract the area of the square
        square = gray_image[y_start:y_end, x_start:x_end]
        
        # Calculate the average gray value of the square
        avg_gray_value = np.mean(square)
        
        # Difference to white (255 is the ideal white value)
        correction_value = 255 - avg_gray_value
        
        # Add the reference point and the difference to the list
        points.append([x_px, y_px])
        values.append(correction_value)

    # Process the square_positions_groups_brightest
    for group in square_positions_groups_brightest:
        brightest_avg_value = -1  # Initial value to compare brightness
        best_point = None

        for (x_mm, y_mm) in group:
            x_px = int(x_mm * width_ratio)
            y_px = int(y_mm * height_ratio)
            
            # Determine the area of the square around this center
            x_start = max(0, x_px - square_size_px // 2)
            y_start = max(0, y_px - square_size_px // 2)
            x_end = min(image_width, x_px + square_size_px // 2)
            y_end = min(image_height, y_px + square_size_px // 2)
            
            # Extract the area of the square
            square = gray_image[y_start:y_end, x_start:x_end]
            
            avg_gray_value = np.mean(square)
            
            # Keep track of the brightest square in the group
            if avg_gray_value > brightest_avg_value:
                brightest_avg_value = avg_gray_value
                best_point = (x_px, y_px)
        
        if best_point:
            correction_value = 255 - brightest_avg_value
            
            points.append([best_point[0], best_point[1]])
            values.append(correction_value)

    # Add artificial edge points to ensure interpolation covers the borders
    edge_points = [
        [0, 0],  # Top-left corner
        [0, image_height-1],  # Bottom-left corner
        [image_width-1, 0],  # Top-right corner
        [image_width-1, image_height-1]  # Bottom-right corner
    ]
    
    # Function to find the nearest point to a given edge point
    def find_nearest_point(edge_point, points):
        distances = np.sqrt((np.array(points)[:, 0] - edge_point[0])**2 + (np.array(points)[:, 1] - edge_point[1])**2)
        nearest_index = np.argmin(distances)
        return values[nearest_index]

    # Use the correction value of the nearest reference point for each corner
    edge_values = [find_nearest_point(ep, points) for ep in edge_points]
    
    # Append the edge points and their corresponding values
    points.extend(edge_points)
    values.extend(edge_values)
    
    print(f"Square center coordinates: {points}")
    print(f"Brightness correction values: {values}")

    # Create a grid for the entire image
    grid_x, grid_y = np.meshgrid(np.arange(image_width), np.arange(image_height))
    
    # Interpolate the correction values over the entire image, including the edges
    correction_matrix = griddata(points, values, (grid_x, grid_y), method='linear', fill_value=0)
    
    # Apply the correction to the grayscale image
    corrected_image = np.clip(gray_image + correction_matrix, 0, 255).astype(np.uint8)
        
    return corrected_image


# Function to draw white corners on the image so that they are not detected as leaves
def draw_white_corners(image, corner_square_length_mm, width_image_mm, height_image_mm, frame_width_mm=5):
    """
    Draws white squares in the four corners of an image based on the specified dimensions 
    in millimeters. This is often used to mark corners for calibration or image transformation.

    Parameters:
    -----------
        image (numpy array): The image on which the white corner squares are to be drawn.
        corner_square_length_mm (float): The side length of the squares in millimeters to be drawn in the corners.
        width_image_mm (float): The width of the image in millimeters.
        height_image_mm (float): The height of the image in millimeters.
        frame_width_mm (float, optional): The frame width around the squares, default is 5 mm.

    Returns:
    --------
        numpy array: The image with white squares in the corners.
    """
    height_pixels, width_pixels = image.shape[:2]
    
    # Convert square length from millimeters to pixels
    square_length_pixels_w = int(np.ceil(((corner_square_length_mm + frame_width_mm) / width_image_mm) * width_pixels))
    square_length_pixels_h = int(np.ceil(((corner_square_length_mm + frame_width_mm) / height_image_mm) * height_pixels))
    
    # Ensure the largest square is chosen to account for rounding differences
    square_length_pixels = max(square_length_pixels_w, square_length_pixels_h)
    
    # Fill the four corners of the image with white squares
    # Top-left corner
    image[0:square_length_pixels, 0:square_length_pixels] = (255, 255, 255)
    # Top-right corner
    image[0:square_length_pixels, -square_length_pixels:] = (255, 255, 255)
    # Bottom-left corner
    image[-square_length_pixels:, 0:square_length_pixels] = (255, 255, 255)
    # Bottom-right corner
    image[-square_length_pixels:, -square_length_pixels:] = (255, 255, 255)
    
    print('Whitened markers successfully')

    return image

# Function to mask image based on threshold
def mask_image(image, mask):
    """
    Masks the input image based on the threshold value applied to the grayscale version.
    Pixels in the grayscale image greater than the threshold are masked (set to black) 
    in the original image.

    Parameters:
    -----------
    - image: np.ndarray
        The input color image (BGR format) as a numpy array, typically read using cv2.imread().
    - mask: np.ndarray
        The mask.
    
    Returns:
    --------
    - masked_image: np.ndarray
        The original image with masked regions (those above the threshold in the grayscale image) set to black.
    """
    
    # Copy the original image to avoid modifying the input
    masked_image = image.copy()
    
    # Set the pixels where mask is True to black [0, 0, 0] in the original image
    masked_image[mask] = [255, 255, 255]
    
    return masked_image

#Function to calculate leaf area in images from a folder
def calculate_leaf_area_in_images(folder_path, total_area_mm2, output_csv):
    """
    This function calculates the leaf area in cm² for all mask images in the specified folder.
    The result is saved into a CSV file with the plot, date, and corresponding leaf area for each image.

    Parameters:
    - folder_path: str
        The path to the folder containing mask images. Only images with filenames ending in 
        '_mask.png', '_mask.jpg', '_mask.jpeg', '_mask.PNG', '_mask.JPG', or '_mask.JPEG' will be processed.
    
    - total_area_mm2: float
        The total area of the image in millimeters squared. The leaf area in each image will be 
        calculated as a proportion of this total area and converted to cm².
    
    - output_csv: str
        The path to the output CSV file where the results (plot, date, and leaf area in cm²) will be saved.

    Returns:
    None
    """

    # Get all image files in the folder that match the mask naming convention
    image_files = [f for f in os.listdir(folder_path) if f.endswith(('_mask.png', '_mask.jpg', '_mask.jpeg', '_mask.PNG', '_mask.JPG', '_mask.JPEG'))]
    
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Plot', 'Date', 'Leaf Area (cm²)'])
        
        print(image_files)  # prints the list of image files for debugging

        for image_file in image_files:
            file_parts = image_file.split('_')
            plot = file_parts[0]  # LT plus plot number
            date = file_parts[1][:6]  # YYMMDD

            image_path = os.path.join(folder_path, image_file)
            
            image = cv2.imread(image_path)
            
            gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
            
            # Count the number of black pixels (assuming leaf area is represented by black pixels with value 255)
            leaf_pixels = np.sum(gray_image == 0)
            
            # Calculate the total number of pixels in the image (width * height)
            total_pixels = image.shape[0] * image.shape[1]
            
            # Calculate the leaf area in mm² and convert to cm²
            leaf_area_cm2 = (leaf_pixels / total_pixels) * total_area_mm2 / 100  # mm² to cm²
            
            # Write the plot, date, and calculated leaf area to the CSV file
            writer.writerow([plot, date, f"{leaf_area_cm2:.2f}"])


###-----------------------------------------------------------------------------------------------------###
### USER INPUT

# input image directory (takes all jpg/jpeg/png images of this directory)
path_dir_in = 'data/LT'

# output directory for processing results
path_dir_out = 'data/LT_processed'

# physical width of the image in mm. Used to convert mm to pixels.
width_mm = 980
# physical height of the image in mm. Used to convert mm to pixels.
height_mm = 1980

# aruco marker size in mm
corner_square_length_mm = 40

# frame of safety around aruco markers in mm
frame_width_mm = 5

# size of the reference squares in mm for exposure correction. 
# These squares are used to determine correction values for exposure based on their average gray values.
square_size_mm = 10

# Fixed squares for expoure correction. List of fixed square center positions (x,y in mm). 
# These squares are directly used as reference points for the interpolation.
square_positions_fix = [
    (frame_width_mm + square_size_mm/2, frame_width_mm + square_size_mm/2 + corner_square_length_mm),
    (frame_width_mm + square_size_mm/2 + corner_square_length_mm, frame_width_mm + square_size_mm/2),
    (width_mm - (frame_width_mm + square_size_mm/2), frame_width_mm + square_size_mm/2 + corner_square_length_mm),
    (width_mm - (frame_width_mm + square_size_mm/2 + corner_square_length_mm), frame_width_mm + square_size_mm/2),                
    (frame_width_mm + square_size_mm/2, height_mm - (frame_width_mm + square_size_mm/2 + corner_square_length_mm)),
    (frame_width_mm + square_size_mm/2 + corner_square_length_mm, height_mm - (frame_width_mm + square_size_mm/2)),
    (width_mm - (frame_width_mm + square_size_mm/2), height_mm - (frame_width_mm + square_size_mm/2 + corner_square_length_mm)),
    (width_mm - (frame_width_mm + square_size_mm/2 + corner_square_length_mm), height_mm - (frame_width_mm + square_size_mm/2))]

# Best-of squares for expoure correction. list of lists, where each inner list contains square center positions (x,y in mm) that belong to a group. 
# For each group, only the brightest square (highest average gray value) will be used as a reference point for the interpolation. 
# This is to avoid reference squares that are covered by leaves accidentially.
square_positions_groups_brightest = [
    # square at center-left
    [(frame_width_mm + square_size_mm/2, height_mm/2),
     (frame_width_mm + square_size_mm/2, height_mm/2 - 50),
     (frame_width_mm + square_size_mm/2, height_mm/2 - 100),
     (frame_width_mm + square_size_mm/2, height_mm/2 + 50),
     (frame_width_mm + square_size_mm/2, height_mm/2 + 100)],
    # square at center-right
    [(width_mm - (frame_width_mm + square_size_mm/2), height_mm/2),
     (width_mm - (frame_width_mm + square_size_mm/2), height_mm/2 -50),
     (width_mm - (frame_width_mm + square_size_mm/2), height_mm/2 -100),
     (width_mm - (frame_width_mm + square_size_mm/2), height_mm/2 + 50),
     (width_mm - (frame_width_mm + square_size_mm/2), height_mm/2 + 100)]
]

# gray threshold for leaf/non-leaf in grayscale image
threshold_gray = 240



### END OF USER INPUT
###--------------------------------------------------------------------###

###Apply on all images in folder
# Execute for all images in folder and subfolders
for root, _, files in os.walk(path_dir_in):  # root gives the current directory path, _ ignores subdirectories, and files is the list of files
    for file_in in files:
        if file_in.split('.')[-1] in ['jpg', 'JPG', 'png', 'PNG', 'jpeg', 'JPEG']:
            print(file_in)
            full_path_in = os.path.join(root, file_in)
            image = cv2.imread(full_path_in)
            image_processed = warp_image(image, width_mm, height_mm, file_in)
            # Whiten ArUco markers in the original image
            image_processed = draw_white_corners(image_processed,
                                                 corner_square_length_mm=corner_square_length_mm,
                                                 width_image_mm=width_mm,
                                                 height_image_mm=height_mm,
                                                 frame_width_mm=frame_width_mm-1) # Avoid overlaps in whitening and future exposure correction squares
            # Correct image exposure
            image_gray = correct_image_exposure(image_processed, width_mm, height_mm, square_size_mm, square_positions_fix, square_positions_groups_brightest)

            file_out_name, file_out_ext = os.path.splitext(file_in)  # os.path.splitext handles file extensions more robustly
            output_gray_path = os.path.join(path_dir_out, f"{file_out_name}_gray{file_out_ext}")
            cv2.imwrite(output_gray_path, image_gray)

            # Create a binary mask where pixels in image_gray greater than the threshold are set to 1 (True)
            mask = image_gray > threshold_gray
            image_mask = np.where(mask, 255, 0).astype(np.uint8)
            output_mask_path = os.path.join(path_dir_out, f"{file_out_name}_mask{file_out_ext}")
            cv2.imwrite(output_mask_path, image_mask)

            # Mask the background
            image_processed = mask_image(image_processed, mask)
            output_processed_path = os.path.join(path_dir_out, f"{file_out_name}_processed{file_out_ext}")
            cv2.imwrite(output_processed_path, image_processed)

# calculate leaf area for each image and save to csv file
calculate_leaf_area_in_images(folder_path=path_dir_out,
                              total_area_mm2=width_mm*height_mm, 
                              output_csv= os.path.join(path_dir_out, 'leaf_area.csv'))
