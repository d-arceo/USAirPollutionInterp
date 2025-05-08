# AUTHOR: DANIEL ARCEO
# KRIGING AND IDW INTERPOLATION

# Load required libraries
library(tidyr)
library(tidyverse)
library(dplyr)
library(gstat)
library(sp)
library(geoR)
library(rgdal)
library(raster)


# Change from scientific notation
options(scipen = 999)

# Read in the cleaned data
pollutant_data <- read.csv("GGRA3_cleaned.csv")

# Create subsets for each of the 3 pollutants

# NITROGEN DIOXIDE (NO2)
nitrogen_dioxide <- filter(pollutant_data, `Parameter.Name` == "Nitrogen dioxide (NO2)")

nitrogen_dioxide <- nitrogen_dioxide %>%
  arrange(Address, 'Event Type') %>%
  filter(duplicated(Address) == FALSE)

# Turn into SpatialPointsDataFrame
coordinates(nitrogen_dioxide) <- ~Longitude+Latitude
proj4string(nitrogen_dioxide) <- CRS("+proj=longlat +datum=WGS84 +no_defs")


# OZONE (O4)
ozone <- filter(pollutant_data, `Parameter.Name` == "Ozone")

ozone <- ozone %>%
  arrange(Address, 'Event Type') %>%
  filter(duplicated(Address) == FALSE)

# Turn into SpatialPointsDataFrame
coordinates(ozone) <- ~Longitude+Latitude
proj4string(ozone) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# PM2.5
pm25 <- filter(pollutant_data, `Parameter.Name` == "PM2.5 - Local Conditions")

pm25 <- pm25 %>%
  arrange(Address, 'Event Type') %>%
  filter(duplicated(Address) == FALSE)

# Turn into SpatialPointsDataFrame
coordinates(pm25) <- ~Longitude+Latitude
proj4string(pm25) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Read in geoJSON for each state as a SpatialPolygonsDataFrame
arizona_spdf <- readOGR("https://raw.githubusercontent.com/georgique/world-geojson/master/states/usa/arizona.json")
colorado_spdf <- readOGR("https://raw.githubusercontent.com/georgique/world-geojson/master/states/usa/colorado.json")
idaho_spdf <- readOGR("https://raw.githubusercontent.com/georgique/world-geojson/master/states/usa/idaho.json")
nevada_spdf <- readOGR("https://raw.githubusercontent.com/georgique/world-geojson/master/states/usa/nevada.json")
utah_spdf <- readOGR("https://raw.githubusercontent.com/georgique/world-geojson/master/states/usa/utah.json")
combined_spdf <- rbind(arizona_spdf, colorado_spdf, idaho_spdf, nevada_spdf, utah_spdf)

# Create the testing grid based on the joined spatial polygons
testing_grid <- as.data.frame(spsample(combined_spdf, "regular", n = 50000))
names(testing_grid) <- c("Longitude", "Latitude")
coordinates(testing_grid) <- c("Longitude", "Latitude")
gridded(testing_grid) <- TRUE
proj4string(testing_grid) <- proj4string(combined_spdf)

# Note: we will use the arithmetic mean as values of interest
nitrogen_dioxide_trans <- nitrogen_dioxide@data %>%
  mutate(logArithmeticMean = log(Arithmetic.Mean+1))
coordinates(nitrogen_dioxide_trans) <- coordinates(nitrogen_dioxide)
proj4string(nitrogen_dioxide_trans) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

ozone_trans <- ozone@data %>%
  mutate(logArithmeticMean = 1/(Arithmetic.Mean+1))
coordinates(ozone_trans) <- coordinates(ozone)
proj4string(ozone_trans) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

pm25_trans <- pm25@data %>%
  mutate(logArithmeticMean = log(Arithmetic.Mean+1))
coordinates(pm25_trans) <- coordinates(pm25)
proj4string(pm25_trans) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Plot histograms of the transformed data -> Ordinary Kriging does not require ND data - but we can check anyway.
hist(nitrogen_dioxide_trans$logArithmeticMean)
hist(ozone_trans$logArithmeticMean)
hist(pm25_trans$logArithmeticMean)

# ASSESS STATIONARITY

# Find extent values
area <- extent(testing_grid@bbox)

# Create a raster of 5x5 covering the area
raster_area <- raster(nrows = 5, ncols = 5, ext = area, vals = -999)

# Convert Raster to Polygons
polygons_area <- rasterToPolygons(raster_area)

# Remove all columns in each subset except logArithmeticMean
nitrogen_dioxide2 <- nitrogen_dioxide_trans[, -(1:54)]
ozone_2 <- ozone_trans[, -(1:54)]
pm25_2 <- pm25_trans[, -(1:54)]

# Number of observations per polygon
assessment_NO2 <- over(polygons_area, nitrogen_dioxide2)
assessment_O4 <- over(polygons_area, ozone_2)
assessment_pm25 <- over(polygons_area, pm25_2)

# Mean of values
assessment_NO2$mean <- over(polygons_area, nitrogen_dioxide2, fn = mean, na.rm = T)[,1]
assessment_O4$mean <- over(polygons_area, ozone_2, fn = mean, na.rm = T)[,1]
assessment_pm25$mean <- over(polygons_area, pm25_2, fn = mean, na.rm = T)[,1]

# Variance of Values
assessment_NO2$var <- over(polygons_area, nitrogen_dioxide2, fn = var, na.rm = T)[,1]
assessment_O4$var <- over(polygons_area, ozone_2, fn = var, na.rm = T)[,1]
assessment_pm25$var <- over(polygons_area, pm25_2, fn = var, na.rm = T)[,1]

# Change Names
names(assessment_NO2)[1] <- c("counts")
names(assessment_O4)[1] <- c("counts")
names(assessment_pm25)[1] <- c("counts")

# Replace Polygon Data
polygons_area@data <- assessment_NO2

# Summary
summary(polygons_area@data) # for NO2
summary(assessment_O4)
summary(assessment_pm25)



# CREATE VARIOGRAMS FOR EACH POLLUTANT -> use transformed values

# nitrogen dioxide variogram
nitrogen_dioxide_vgm <- variogram(logArithmeticMean ~ 1, data = nitrogen_dioxide_trans)

# ozone variogram
ozone_vgm <- variogram(logArithmeticMean ~ 1, data = ozone_trans)

# pm25 variogram
pm25_vgm <- variogram(logArithmeticMean ~ 1, data = pm25_trans)

# FITTING VARIOGRAMS

# fitting nitrogen dioxide variogram using spherical variogram model
nitrogen_dioxide_fit <- fit.variogram(nitrogen_dioxide_vgm, model = vgm("Sph"))
plot(nitrogen_dioxide_vgm, nitrogen_dioxide_fit)

# fitting ozone dioxide variogram using 
ozone_fit <- fit.variogram(ozone_vgm, model = vgm(psill = 0.000025, "Sph"))
plot(ozone_vgm, ozone_fit)

# fitting pm25 variogram using spherical variogram model
pm25_fit <- fit.variogram(pm25_vgm, model = vgm("Sph"))
plot(pm25_vgm, pm25_fit)

# Kriging Setup
nitrogen_dioxide_krige <- krige(logArithmeticMean ~ 1, nitrogen_dioxide_trans, testing_grid, model = nitrogen_dioxide_fit)
ozone_krige <- krige(logArithmeticMean ~ 1, ozone_trans, testing_grid, model = ozone_fit)
pm25_krige <- krige(logArithmeticMean ~ 1, pm25_trans, testing_grid, model = pm25_fit)

# Nitrogen dioxide krig plot
nitrogen_dioxide_krige %>% as.data.frame %>%
  ggplot(aes(x=Longitude, y=Latitude)) + geom_tile(aes(fill = var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_x_continuous() + 
  scale_y_continuous() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("Nitrogen Dioxide (NO2) Concentration Interpolation in Study Area", subtitle = "Using Ordinary Kriging")

# Ozone krig plot
ozone_krige %>% as.data.frame %>%
  ggplot(aes(x = Longitude, y = Latitude)) + geom_tile(aes(fill = var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_x_continuous() +
  scale_y_continuous() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("Ozone (04) Concentration Interpolation in Study Area", subtitle = "Using Ordinary Kriging")

#PM2.5 krig plot
pm25_krige %>% as.data.frame %>%
  ggplot(aes(x = Longitude, y = Latitude)) + geom_tile(aes(fill = var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_x_continuous() +
  scale_y_continuous() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("PM2.5 Concentration Interpolation in Study Area", subtitle = "Using Ordinary Kriging")

#PM2.5 IDW plot
idwems_nit <- idw(Arithmetic.Mean ~ 1, nitrogen_dioxide, testing_grid, idp = 2)
idwems_nit %>% as.data.frame %>%
  ggplot(aes(x=Longitude, y=Latitude)) + geom_tile(aes(fill = var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_x_continuous() + 
  scale_y_continuous() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("Nitrogen Dioxide (NO2) Concentration Interpolation in Study Area", subtitle = "Using IDW Interpolation")

#PM2.5 IDW plot
idwems_oz <- idw(Arithmetic.Mean ~ 1, ozone, testing_grid, idp = 2)
idwems_oz %>% as.data.frame %>%
  ggplot(aes(x=Longitude, y=Latitude)) + geom_tile(aes(fill = var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_x_continuous() + 
  scale_y_continuous() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("Ozone Concentration Interpolation in Study Area", subtitle = "Using IDW Interpolation")


#PM2.5 IDW plot
idwems_pm <- idw(Arithmetic.Mean ~ 1, pm25, testing_grid, idp = 2)
idwems_pm %>% as.data.frame %>%
  ggplot(aes(x=Longitude, y=Latitude)) + geom_tile(aes(fill = var1.pred)) + coord_equal() +
  scale_fill_gradient(low = "yellow", high = "red") +
  scale_x_continuous() + 
  scale_y_continuous() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("particulate matter2.5 Concentration Interpolation in Study Area", subtitle = "Using IDW Interpolation")
