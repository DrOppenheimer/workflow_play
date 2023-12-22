# for norm with quantiles
BiocManager::install("preprocessCore")
library(preprocessCore)


# get the metadata subset that corresponds to the small HMP sample
# did this just once - not meant for repeating
metadata_subset <- matrix( nrow=30, ncol=76 )
colnames(metadata_subset) <- colnames(my_metadata)
rownames <- list()
for ( i in 1:ncol(my_data)){ # my data is the subset of HMP
  #print(colnames(my_data)[i]) 
  metadata_subset[i,] <- my_metadata[ colnames(my_data)[i], ]
  rownames[i] <- colnames(my_data)[i]
  }
rownames(metadata_subset) <- rownames
# export data as a text file
source("~/Documents/GitHub/Kevin_R_scripts/export_data.r")
export_data(metadata_subset, "filtered_counts.metadata.txt")
##################################################################

# Start here --------------------------------------------------------------

# set the working directory
setwd("~/Documents/GitHub/workflow_play")

# Preprocess the data
source("preprocessing_tool.r")
preprocessing_tool("filtered_counts.txt", )

# calculate PCoA on the preprocessed data
source("calculate_pco.r")
calculate_pco(file_in="filtered_counts.txt.standardize.PREPROCESSED.txt")

# render interactive 3d PCoA plot from the PCoA and the corresponding metadata
source("render_calculated_pcoa.r")
plot_interactive_colored_3d_pcoa(
  pcoa_data_file = "filtered_counts.txt.standardize.PREPROCESSED.txt.euclidean.PCoA", #my_test_pcoa,
  selected_eigen_vectors = c(1,2,3),
  pcoa_metadata_file = "filtered_counts.metadata.txt",
  metadata_column = "env_package.data.body_site"#metadata_colors
)

# iterate through the metadata creating a static 3d plot for each metadata column
plot_static_colored_3d_pcoas(
  pcoa_filename = "filtered_counts.txt.standardize.PREPROCESSED.txt.euclidean.PCoA",
  selected_eigen_vectors = c(1,2,3),
  metadata_filename = "filtered_counts.metadata.txt",
  debug = TRUE
)
system("open filtered_counts.txt.standardize.PREPROCESSED.txt.euclidean.PCoA.*.png")

# perform a stat test on the data/metadata (selected metadata column is used to create groups for chosen stat test)
source("sigtest.R")
sigtest(data_file="filtered_counts.txt", 
                    metadata_file="filtered_counts.metadata.txt",  
                    metadata_column="env_package.data.body_site", 
                    stat_test="ANOVA-one-way",
                    p_adjust_method = "BH"
  )

# Load the stat test results and subselect data based on the stat results
library(tidyverse)
my_stat_results <- import_data("filtered_counts.ANOVA-one-way.env_package.data.body_site.STAT_RESULTS.txt")
my_stat_results <- as_tibble(my_stat_results,rownames=NA)
my_stat_results_subselected <- my_stat_results |> filter( bonferroni_p < 0.0001 ) # woot, I used a pipe
# export results with stats
export_data(data_object = my_stat_results_subselected, file_name = "my_stat_results_subselected.txt")
# remove the columns that contain the stats(this is hacky)
source("remove_last_n_columns.R")
my_stat_results_subselected_and_cleaned <- remove_last_n_columns(my_stat_results_subselected, n=7)
# export the data ready to create a HD
export_data(data_object = my_stat_results_subselected_and_cleaned, file_name = "my_stat_results_subselected_and_cleaned.txt")


# create heatmap dendrograms of original and stat subselected data
source("heatmap_dendrogram.r")
heatmap_dendrogram(file_in = "filtered_counts.txt", # should really be using normalized data here
                   metadata_table = "filtered_counts.metadata.txt",
                   metadata_column="env_package.data.body_site"
  )
system("open filtered_counts.txt.HD.png")
heatmap_dendrogram(file_in = "my_stat_results_subselected_and_cleaned.txt",
                   metadata_table = "filtered_counts.metadata.txt",
                   metadata_column="env_package.data.body_site"
)
system("open my_stat_results_subselected_and_cleaned.txt.HD.png")






