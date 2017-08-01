## Evaluating spatial autocorrelation and distances among locations in the vampire data
## 
## Mark A. Hayes
## 8/1/2017

## Install "geoR", "sp", "rgeos", "geosphere"  packages

install.packages("geoR") # if the package is not yet installed

library(geoR)
library(sp)
library(rgeos)
library(geosphere)

## A note on the data used:

# Here, we use the merged dataset (mds) file created by SAHM, which is the "Hayes and Piaggio data.csv".  
# The latitude and longitude for each occurrence location is shown in  X and Y columns
# of the csv. These coordinates are in decimal degrees. We first convert the latitude  
# and longitude coordinates into distances in kilometers from an origin for each location. 
# These distances will be labeled "x_km" and "y_km", and will represent the approximate 
# distance from the origin (the minimum latitude and longitude) in km of x and y for each 
# location. We use the following functions to calculate these values, then add them
# as new columns to the dataframe using the "cbind" command:
# x_km: =111*(X-(-111))*COS(Y/(360/(2*PI())))
# y_km: =111*(Y-14)

## Housekeeping:

# Remove any objects as needed.

ls()
remove()

rm(list=ls())

# Bring in the data file ("Hayes and Piaggio data.csv"):

data = read.table(file.choose(), header=TRUE, sep=",")

summary(data) # Summarize the data

## Converting decimal degrees to km, and adding "x_km" and "y_km" to dataframe

x_km <- 111*(data$X-(-111))*cos(data$Y/(360/(2*pi)))
y_km <- 111*(data$Y-14)

# Bind these new vectors to the original dataframe

df <- cbind(data, x_km)
df <- cbind(df, y_km)

# Now the "df" dataframe includes columns for "x_km" and y_km".

# In this dataframe x_km and y_km are the x and y coordinates.
# Plotting these:

plot(df$x_km, df$y_km)

# Printing the plot, as desired, for example, as a tiff file:

dev.print(tiff, "locs.tiff", height=4, width=6, units='in', res=300)

## Evaluating distances among locations, in km:

# The "dist" command calculates distanced between all pairs of points.
# The [,9:10] command indicates that the spatial locations are in column 9 & 10.
# The "hist" command creates a histogram for visualization of these distances.

dists = dist(df[,9:10])
summary(dists)
hist(dists)

# How many occurrence records in the final data set?

length(df$bio6)

# How many items in the "dists" vector? This should be the same length as the columns in "df".

length(dists)

# How many items in "dists" are less than 5 km?

sum(dists<5.0) 

# Note that here "sum" adds the number of items in the vector, but doesn't add 
# the vector values together. Crawley (pp 39) discusses how to do this, if needed. 

# Percent of dists that are < 5.0 km?

(sum(dists<5.0)/length(dists))*100

# Percent of dists that are less than 100 km? 

(sum(dists<100)/length(dists))*100 # < 5%, so spatial autocorrelation should not a problem. 

## Adding a column to the dataframe with minimum distance to nearest point, for each location:

# Here, we use the raw decimal degree coordinates to calculate distance in kilometers,
# use package "geosphere".

# Transform the dataframe to spatial objects:
  
sp.df <- df
coordinates(sp.df) <- ~X+Y

class(sp.df) 
attr(,"package") # Checking that class is a spatial dataframe, and the package used

# Now calculate pairwise distances between points

d <- distm(sp.df)

# Finding the second shortest distance (closest distance is of point to itself, 
# therefore use second shortest)

min_d <- apply(d, 1, function(x) order(x, decreasing=F)[2])

# Bind the "min_d" vector to the "df" dataframe, which now will have 11 columns.
# But first checking to be sure the vectors are the same length.
length(df$X)
length(min_d)
df <- cbind(df,min_d)

# Sumarize the "min_d" vector, then visualize as a histogram.

summary(min_d)
hist(min_d)

# How many occurrence records in the final data set?

length(df$X)

# How many items in the "min_d" vector?

length(min_d)

# How many items in "min_d" are less than 5 km?

sum(min_d<5.0) 

# Percent of min_d that are < 5.0?

(sum(min_d<5.0)/length(min_d))*100

# What % less than 50 km?

(sum(min_d<50)/length(min_d))*100 # < 5%, so spatial autocorrelation probably not a problem. 


## Evaluating patterns in spatial autocorrelation.  
# Plotting variograms for each variable, bio2, bio6, bio12, bio15, bio19

bio2_variog <- variog(coords = df[,9:10], data = df$bio2)
plot(bio2_variog, main="bio2 variogram")

bio6_variog <- variog(coords = df[,9:10], data = df$bio6)
plot(bio6_variog, main="bio6 variogram")

bio12_variog <- variog(coords = df[,9:10], data = df$bio12)
plot(bio12_variog, main="bio12 variogram")

bio15_variog <- variog(coords = df[,9:10], data = df$bio15)
plot(bio15_variog, main="bio15 variogram")

bio19_variog <- variog(coords = df[,9:10], data = df$bio19)
plot(bio19_variog, main="bio19 variogram")

# Printing the variograms, as desired, for example:

dev.print(tiff, "figure1_bio6.tiff", height=4, width=6, units='in', res=300)

## End
