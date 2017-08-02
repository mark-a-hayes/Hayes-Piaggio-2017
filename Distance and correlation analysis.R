## Evaluating spatial autocorrelation and distances among locations in the vampire data
## 
## Mark A. Hayes
## 8/2/2017

## Install "spatstat" package

install.packages("spatstat") # if the package is not yet installed

library(spatstat)

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

## Looking at nearest neighbor distances (Crawley pp. 830+)

# How many occurrence records in the final data set?

length(df$bio6)

# Setting up a loop to calculate distances

x <- df$x_km
y <- df$y_km

distance <- function(x1, y1, x2, y2) sqrt ((x2 - x1)^2 + (y2 - y1)^2)

r <- numeric(1145)
nn <- numeric(1145)
d <- numeric(1145)
for (i in 1:1145) {
  for (k in 1:1145) d[k] <- distance(x[i],y[i],x[k],y[k])
  r[i] <- min(d[-i])
  nn[i] <- which (d==min(d[-i]))
}

for (i in 1:1145) lines(c(x[i], x[nn[i]]), c(y[i], y[nn[i]]), col = "green")

# Printing the plot, as desired, for example, as a tiff file:

dev.print(tiff, "locs_nn.tiff", height=4, width=6, units='in', res=300)

## Creating a pair correlation function (Crawley pp. 844+)

# Convert coordinate data to a point pattern data set in a 2-dimentional plane:

summary(df)

vampires <- ppp(df$x_km, df$y_km, c(0,2500), c(0,1550), marks = df$response)

plot(vampires)

summary(vampires)

pc <- pcf(vampires)
plot(pc, main = "Pair correlation function")
plot(pc, main = "Pair correlation function", xlim = c(0,100)) # looking closer
plot(pc, main = "Pair correlation function", xlim = c(0,50)) # and closer still
plot(pc, main = "Pair correlation function", xlim = c(0,20)) # and still closer

# Printing the plot, as desired, for example, as a tiff file:

dev.print(tiff, "pcf.tiff", height=4, width=6, units='in', res=300)

# The function 'distmap' shows the distance map around individual vampire occurrences

Z <- distmap(vampires) 
plot(Z, main = "")
dev.print(tiff, "distmap.tiff", height=4, width=6, units='in', res=300)

neighbors <- nndist(vampires)
neighbors
summary(neighbors)
hist(neighbors)
length(neighbors)
df <- cbind(df, neighbors)
sum(neighbors < 5.0)
(sum(neighbors<5.0)/length(neighbors))*100

## End


