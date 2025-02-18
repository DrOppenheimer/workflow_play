# WORKFLOW ----------------------------------------------------------------

# load necessary packages 
library(plotly)
library(scatterplot3d)

# FUNCTIONS ---------------------------------------------------------------
# FUNCTIONS ---------------------------------------------------------------
# FUNCTIONS ---------------------------------------------------------------

# Import PCoA
## ######################
## # SUB(1): FUNCTION TO LOAD A PRECALCULATED *.PCoA
## ######################
load_pcoa_data <- function(PCoA_in){
  # create two connections to the file
  con_1 <- file(PCoA_in)
  con_2 <- file(PCoA_in)
  # read through the first time to get the number of samples
  open(con_1);
  num_values <- 0
  data_type = "NA"
  while ( length(my_line <- readLines(con_1,n = 1, warn = FALSE)) > 0) {
    if ( length( grep("PCO", my_line) ) == 1  ){
      num_values <- num_values + 1
    }
  }
  close(con_1)
  # create object for values
  eigen_values <- matrix("", num_values, 1)
  dimnames(eigen_values)[[1]] <- 1:num_values
  # create object for vectors
  eigen_vectors <- matrix("", num_values, num_values)
  dimnames(eigen_vectors)[[1]] <- 1:num_values
  # read through the input file a second time to populate the R objects
  value_index <- 1
  vector_index <- 1
  open(con_2)
  current.line <- 1
  data_type = "NA"
  
  while ( length(my_line <- readLines(con_2,n = 1, warn = FALSE)) > 0) {
    if ( length( grep("#", my_line) ) == 1  ){
      if ( length( grep("EIGEN VALUES", my_line) ) == 1  ){
        data_type="eigen_values"
      } else if ( length( grep("EIGEN VECTORS", my_line) ) == 1 ){
        data_type="eigen_vectors"
      }
    }else{
      split_line <- noquote(strsplit(my_line, split="\t"))
      if ( identical(data_type, "eigen_values")==TRUE ){
        dimnames(eigen_values)[[1]][value_index] <- noquote(split_line[[1]][1])
        eigen_values[value_index,1] <- noquote(split_line[[1]][2])       
        value_index <- value_index + 1
      }
      if ( identical(data_type, "eigen_vectors")==TRUE ){
        dimnames(eigen_vectors)[[1]][vector_index] <- noquote(split_line[[1]][1])
        for (i in 2:(num_values+1)){
          eigen_vectors[vector_index, (i-1)] <- as.numeric(noquote(split_line[[1]][i]))
        }
        vector_index <- vector_index + 1
      }
    }
  }
  close(con_2)
  
  # finish labeling of data objects and return them in a single list object
  dimnames(eigen_values)[[2]] <- "EigenValues"
  dimnames(eigen_vectors)[[2]] <- dimnames(eigen_values)[[1]]
  class(eigen_values) <- "numeric"
  class(eigen_vectors) <- "numeric"
  return(list(eigen_values=eigen_values, eigen_vectors=eigen_vectors))
  
}



# Import metadata
## ######################
## # SUB(2): FUNCTION TO LOAD METADATA FILE
## ######################
load_metadata <- function(metadata_file){
  metadata_matrix <- as.matrix( # Load the metadata table (same if you use one or all columns)
    read.table(
      file=metadata_file,row.names=1,header=TRUE,sep="\t",
      colClasses = "character", check.names=FALSE,
      comment.char = "",quote="",fill=TRUE,blank.lines.skip=FALSE
    )
  ) 
  return(metadata_matrix)
}




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


## ######################
## # SUB(4): FUNCTION TO GENERATE THE INTERACTIVE 3D PLOT 
## ######################
# Create an interactive 3D scatter plot using the selected_eigen_vectors and metadata_colors
plot_interactive_colored_3d_pcoa <- function(
    pcoa_data_file = "data.PCoA",
    selected_eigen_vectors = c(1,2,3),
    pcoa_metadata_file = "metadata.txt",
    metadata_column = 1
){
  
  #import PCoA data file here
  pcoa_data <- load_pcoa_data(pcoa_data_file)
  
  #import metadata file here
  my_metadata <- load_metadata(pcoa_metadata_file)
  
  # generate metadata colors here
  metadata_column = my_metadata[,metadata_column, drop=FALSE]
  
  # get metadata column name
  metadata_column_name <- colnames(metadata_column) 

  # create a data frame with a category axis that contains the PCoAs and the metadata
  my_df <- cbind(
    metadata_column,
    pcoa_data$eigen_vectors[,selected_eigen_vectors[1]],
    pcoa_data$eigen_vectors[,selected_eigen_vectors[2]],
    pcoa_data$eigen_vectors[,selected_eigen_vectors[3]]
  )

  # add column names to the df
  colnames(my_df) <- c(metadata_column_name, 
                     paste("PCoA_",selected_eigen_vectors[1],"_(", round(pcoa_data$eigen_values[selected_eigen_vectors[1]], digits = 4),")", sep=""),
                     paste("PCoA_",selected_eigen_vectors[2],"_(", round(pcoa_data$eigen_values[selected_eigen_vectors[2]], digits = 4),")", sep=""),
                     paste("PCoA_",selected_eigen_vectors[3],"_(", round(pcoa_data$eigen_values[selected_eigen_vectors[3]], digits = 4),")", sep="")
                     )

  # Extract the unique categories
  categories <- unique(my_df[,1])  # Extract unique categories

  # Create colors for the plot
  metadata_colors <- create_colors(matrix(categories, ncol = 1))

  # Create an empty plot
  fig <- plot_ly()

  # Add each category as a separate trace with custom legend name and color
  for (i in seq_along(categories)) {
    cat_data <- my_df[my_df[,1] == categories[i], ]  # Filter category data
    fig <- fig %>%
        add_trace(
        x = cat_data[,2], 
        y = cat_data[,3], 
        z = cat_data[,4], 
        type = "scatter3d", 
        mode = "markers",
        marker = list(color = metadata_colors[i,1], size = 10),  # Custom color
        name = legend_labels[i]  # Custom legend text
        )
    }

  # Customize layout
  fig <- fig %>%
    layout(
        title = paste("3d PCoA colored by", metadata_column_name),
        scene = list(#xaxis = list(title = colnames(my_df)[2], tickformat = ".2f"),
            xaxis = list(title = colnames(my_df)[2]), 
            yaxis = list(title = colnames(my_df)[3]), 
            zaxis = list(title = colnames(my_df)[4])),
            legend = list(title = list(text = 'Legend'))
        )
    

 fig

}
######################
######################


## ######################
## # SUB(5): FUNCTION TO GENERATE STATIC 3D PLOTS FROM ALL OF THE METADATA 
## ######################
######################
# SUB(1): Workhorse function that creates the 3d plot
######################
plot_static_colored_3d_pcoas <- function(
    pcoa_filename = my_test_pcoa,
    selected_eigen_vectors = c(1,2,3),
    metadata_filename = my_test_metadata,
    debug = TRUE
){
  # import the calculated PCoA
  pcoa_data <- load_pcoa_data(pcoa_filename)
  
  # import metadata
  my_metadata <- load_metadata(metadata_filename)
  
  # iterate through the metadata
  # create plot and color it automatically
  for (i in 1:ncol(my_metadata)){
    
    metadata_column_name <- colnames(my_metadata)[i]
    metadata_column = my_metadata[,metadata_column_name, drop=FALSE]
    my_output_png_filename <- paste(pcoa_filename, ".", metadata_column_name, ".png", sep = "")
    #print(metadata_column_name)
    column_color = create_colors(metadata_column)
    if(debug==TRUE){ print(paste(dim(column_color))); testing <<- column_color }
    
    png(
      filename=my_output_png_filename, 
      width = 8, 
      height = 8, 
      units = "in", 
      res = 300) # 300 dpi
    scatterplot3d(
      x = pcoa_data[["eigen_vectors"]][,selected_eigen_vectors[1]], 
      y = pcoa_data[["eigen_vectors"]][,selected_eigen_vectors[2]], 
      z = pcoa_data[["eigen_vectors"]][,selected_eigen_vectors[3]], 
      color = column_color[,1], 
      pch = 19, 
      main = metadata_column_name,
      xlab = paste("xPC_",selected_eigen_vectors[1],"_", pcoa_data[["eigen_values"]][selected_eigen_vectors[1]]),
      ylab = paste("yPC_",selected_eigen_vectors[2],"_", pcoa_data[["eigen_values"]][selected_eigen_vectors[2]]), 
      zlab = paste("zPC_",selected_eigen_vectors[3],"_", pcoa_data[["eigen_values"]][selected_eigen_vectors[3]]))
    dev.off()
    
  }
  
}

## import the calculated PCoA
#my_test_pcoa <- load_pcoa_data("HMP.Jumpstart.DESeq_normed.euclidean.PCoA")

## import metadata
#my_metadata <- load_metadata("HMP_jumpstart_metadata.txt")

## select a single column of metadata for color generation
#column_name = "env_package.data.body_site"
##column_name = "mixs.seq_method"
##column_name = "library.data.seq_center"
#metadata_column = metadata_matrix[,column_name, drop=FALSE]
## generate colors for selected column
#metadata_colors <- create_colors(metadata_column)

## generate the interactive 3d plot
#plot_interactive_colored_3d_pcoa(
#  pcoa_data_file = "HMP.Jumpstart.DESeq_normed.euclidean.PCoA", #my_test_pcoa,
#  selected_eigen_vectors = c(1,2,3),
#  pcoa_metadata_file = "HMP_jumpstart_metadata.txt",
#  metadata_column = "env_package.data.body_site"#metadata_colors
#)

## iterate through the metadata creating a static 3d plot for each metadata column
#plot_static_colored_3d_pcoas(
#  pcoa_filename = "HMP.Jumpstart.DESeq_normed.euclidean.PCoA",
#  selected_eigen_vectors = c(1,2,3),
#  metadata_filename = "HMP_jumpstart_metadata.txt",
#  debug = TRUE
#)






























############ JUST CRAP BELOW HERE
############ JUST CRAP BELOW HERE
############ JUST CRAP BELOW HERE





