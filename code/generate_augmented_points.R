# A Novel Approach to Field Data Augmentation with Remote Sensing and Machine Learning in Rangelands
# R Code by Javier Osorio Leyton & Hailey E. Schmidt
# Part II: Augmenting Field Data
# -------------------------------------------------------

# load libraries
library(sf)
library(dplyr)

# load the ground truth shapefile
ground_truth_data <- st_read("path_to/ground_truth_data.shp")

# create offsets for neighboring pixels
offsets <- data.frame(
  x_offset = c(3, -3, 0, 0),  # right, left, top, bottom
  y_offset = c(0, 0, 3, -3)   # right, left, top, bottom
)

# function to generate neighbors for a single point
generate_neighbors <- function(point, attributes, offsets) {
  # extract coordinates
  coords <- st_coordinates(point)
  # create shifted points using offsets
  neighbor_points <- lapply(1:nrow(offsets), function(i) {
    st_point(c(coords[1] + offsets$x_offset[i], coords[2] + offsets$y_offset[i]))
  })
  
  # create a dataframe with the same attributes as the original point
  neighbor_attributes <- do.call(rbind, replicate(length(neighbor_points), attributes, simplify = FALSE))
  
  # return as an sf object
  st_sf(
    neighbor_attributes,
    geometry = st_sfc(neighbor_points, crs = st_crs(point))
  )
}

# apply the function to all points
neighbor_data <- do.call(rbind, lapply(1:nrow(ground_truth_data), function(i) {
  generate_neighbors(
    point = st_geometry(ground_truth_data)[i],
    attributes = as.data.frame(ground_truth_data[i, , drop = FALSE]),
    offsets = offsets
  )
}))

# combine original points with neighbors
augmented_data <- rbind(ground_truth_data, neighbor_data) %>%
  distinct(geometry, .keep_all = TRUE)

# plot the original shapefile with a base raster
plot(b2data$blue, main = "Shapefile Overlay")
plot(st_geometry(augmented_data), 
     pch = 19,
     cex = 0.5,
     col = "red",
     add = TRUE)

# check if they are identical
identical(augmented_data$class, augmented_data$Class)

# drop one column (e.g., 'Class')
augmented_data <- augmented_data %>% select(-Class)

# save shapefile - simple (writes four files with extension: .dbf; .prj; .shx; .shp)
st_write(obj = augmented_data, # sf_object,
         dsn = "augmented_field_data.shp",
         quiet = FALSE,
         factorsAsCharacter = TRUE,
         driver = "ESRI Shapefile",
         overwrite = TRUE)

# important note: this just generates neighbor points for each field point
# the shapefile should be loaded into GIS for inspection/quality checks 
# as many points will likely need to be deleted to avoid misclassification 

# also note: once the cleaned shapefile is ready, it can be loaded into Part I 
# the same as the other points datasets for sample extraction



