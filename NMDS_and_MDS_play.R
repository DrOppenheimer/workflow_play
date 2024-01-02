
# SET THE STAGE -----------------------------------------------------------


# set the workign directory
setwd("~/Documents/GitHub/workflow_play/")

# simple function to import data
import_data <- function(file_name)
{
  data.matrix(read.table(file_name, 
                         row.names=1, 
                         header=TRUE, 
                         sep="\t", 
                         comment.char="", 
                         quote="", check.names=FALSE))
}

# simple function to import metadata
import_metadata <- function(group_table){ #, group_column, sample_names){
  metadata_matrix <- as.matrix( # Load the metadata table (same if you use one or all columns)
    read.table(
      file=group_table,row.names=1,header=TRUE,sep="\t",
      colClasses = "character", check.names=FALSE,
      comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE
    )
  )
}





# NMDS --------------------------------------------------------------------
my_data <- import_data("filtered_counts.txt")
my_metadata <- import_metadata("filtered_counts.metadata.txt")

# data expected as a sample per ROW so have to rotate it
colnames(my_data)
library(matlab)
my_data <- rot90(my_data)
rownames(my_data)

# Load the required packages 
library(vegan) 
library(tidyverse)
# Create a sample data matrix 
# data <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3) 
# Compute the distance matrix 
dist_matrix <- vegdist(my_data) 
# method # bray by default	
# Dissimilarity index, partial match to 
###"manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", 
###"jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", 
###"binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", 
###"aitchison", or "robust.aitchison".
# Perform NMDS on the distance matrix 
nmds_result <- metaMDS(dist_matrix, k=2) # use k=3 for 3d 
# Print the NMDS result 
nmds_result
nmds_result$points
# take a quick look
plot(x =nmds_result$points[,"MDS1"] , y =nmds_result$points[,"MDS2"])
# prep tidy data for ggplot
# make sure the data and metadata are sorted by rowname
# Sort the rows based on row names
sorted_metadata <- my_metadata[order(rownames(my_metadata)), ]
results <- cbind(nmds_result$points[,"MDS1"], nmds_result$points[,"MDS2"])
sorted_results <- results[order(rownames(results)), ]
# make the output tidy for ggplot
my_tidy_data <- cbind(sorted_results, sorted_metadata[,"env_package.data.body_site"])
colnames(my_tidy_data) <- c("MDS1","MDS2", "env_package.data.body_site" )
my_tidy_data[1:3,]
my_tidy_data <- as.data.frame(my_tidy_data)
# plot the results with ggplot
ggplot(
  data = my_tidy_data,
  mapping = aes(x=MDS1, y=MDS2, color=env_package.data.body_site)) +
    geom_point()



# Do a 3d NMDS -- requires functions below to generate color from metadata

dist_matrix <- vegdist(my_data) 
nmds_result <- metaMDS(dist_matrix, k=3) # use k=3 for 3d 
results <- cbind(nmds_result$points[,"MDS1"], nmds_result$points[,"MDS2"],nmds_result$points[,"MDS3"])
sorted_results <- results[order(rownames(results)), ]
# make the output tidy for ggplot
my_tidy_data <- cbind(sorted_results, sorted_metadata[,"env_package.data.body_site"])
colnames(my_tidy_data) <- c("MDS1","MDS2","MDS3", "env_package.data.body_site" )
my_tidy_data[1:3,]
my_tidy_data <- as.data.frame(my_tidy_data)

# plot it WITHOUT GGPLOT
# Create colors from the metadata
# generate metadata colors here
metadata_column = my_tidy_data[,'env_package.data.body_site', drop=FALSE]
metadata_column <- metadata_column[order(rownames(metadata_column)),,drop=FALSE]
# generate colors for selected column
metadata_colors <- create_colors(metadata_column)


# render 3d plot
plot_ly( x=my_tidy_data$MDS1,
         y=my_tidy_data$MDS2,
         z=my_tidy_data$MDS3,
         type = "scatter3d",
         mode = "markers",
         #marker = list(color="red")  # list(color = metadata_colors[,1]))
         marker = list(color = metadata_colors[,1])
)








# MDS ---------------------------------------------------------------------

# import data and metadata and make sure that both are similarly sorted by mgid
my_data <- import_data("filtered_counts.txt")
my_metadata <- import_metadata("filtered_counts.metadata.txt")
identical(colnames(my_data), rownames(my_metadata))
my_data <- my_data[,order(colnames(my_data))]
my_metadata <- my_metadata[order(rownames(my_metadata)),]
identical(colnames(my_data), rownames(my_metadata))





# data expected as a sample per ROW so have to rotate it
colnames(my_data)
library(matlab)
my_data <- rot90(my_data)
rownames(my_data)
colnames(my_data)

# create dist matrix 
dist_matrix <- vegdist(my_data) 

# Perform Multidimensional Scaling (MDS)
mds_result <- cmdscale(dist_matrix, k = 3)  # k=3for 3D visualization

# Extract coordinates for plotting
mds_coordinates <- as.data.frame(mds_result)
mds_coordinates[1:3,]

# Quick Plot of the MDS result
plot(mds_coordinates$V1, mds_coordinates$V2, type = "n", xlab = "Dimension 1", ylab = "Dimension 2", main = "MDS Plot")
text(mds_coordinates$V1, mds_coordinates$V2, labels = rownames(mds_coordinates))

# make tidy and prep for scatterplot 3d and plotly
# first make sure data and metadata are sorted the same
my_data <- mds_coordinates[order(rownames(mds_coordinates)),]
my_metadata <- my_metadata[order(rownames(my_metadata)),]
# combine data and metadata to make it tidy
my_tidy_data <- cbind(my_data,my_metadata[,"env_package.data.body_site"])
colnames(my_tidy_data) <- c("V1","V2","V3", "env_package.data.body_site" )

# make a 2 d plot with ggplot
ggplot(
  data = my_tidy_data,
  mapping = aes(x=V1, y=V2, color=env_package.data.body_site)) + 
  geom_point()


# Now do an interactive 3d plot
library(plotly)
library(scatterplot3d)


# FUNCTIONS TO GENERATE COLOR FROM METADATA -------------------------------

# Functions to create colors from metadata
## ######################
## # SUB(3): FUNCTION TO GENERATE COLORS FOR THE PLOT 
## ######################
create_colors <- function(metadata_column, color_mode = "auto"){ # function to     
  my_data.color <- data.frame(metadata_column)
  column_factors <- as.factor(metadata_column[,1])
  column_levels <- levels(as.factor(metadata_column[,1]))
  num_levels <- length(column_levels)
  color_levels <- col.wheel(num_levels)
  levels(column_factors) <- color_levels
  my_data.color[,1]<-as.character(column_factors)
  return(my_data.color)
}

######################
# SUB(8): Create optimal contrast color selection using a color wheel
# adapted from https://stat.ethz.ch/pipermail/r-help/2002-May/022037.html 
######################
col.wheel <- function(num_col, my_cex=0.75) {
  cols <- rainbow(num_col)
  col_names <- vector(mode="list", length=num_col)
  for (i in 1:num_col){
    col_names[i] <- getColorTable(cols[i])
  }
  cols
}
######################
######################
######################
# SUB(10): Convert all colors into format "#rrggbb"
# adapted from https://stat.ethz.ch/pipermail/r-help/2002-May/022037.html
######################
getColorTable <- function(col) {
  rgb <- col2rgb(col);
  col <- rgb2col(rgb);
  sort(unique(col))
}
######################
######################
######################
# SUB(9): The inverse function to col2rgb()
# adapted from https://stat.ethz.ch/pipermail/r-help/2002-May/022037.html
######################
rgb2col <- function(rgb) {
  rgb <- as.integer(rgb)
  class(rgb) <- "hexmode"
  rgb <- as.character(rgb)
  rgb <- matrix(rgb, nrow=3)
  paste("#", apply(rgb, MARGIN=2, FUN=paste, collapse=""), sep="")
}
######################
######################


# Create colors from the metadata
# generate metadata colors here
metadata_column = my_tidy_data[,'env_package.data.body_site', drop=FALSE]
metadata_column <- metadata_column[order(rownames(metadata_column)),,drop=FALSE]
# generate colors for selected column
metadata_colors <- create_colors(metadata_column)

  
# render 3d plot
plot_ly( x=my_tidy_data$V1,
         y=my_tidy_data$V2,
         z=my_tidy_data$V3,
         type = "scatter3d",
         mode = "markers",
         #marker = list(color="red")  # list(color = metadata_colors[,1]))
         marker = list(color = metadata_colors[,1])
)


