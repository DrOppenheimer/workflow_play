# set the working directory
setwd("~/Documents/GitHub/workflow_play")


# Install BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# Use bioconductor to install DESeq2 package
BiocManager::install("DESeq2")

# Load package
library(DESeq2)

# Simple way from ChatGPT
# Assuming you have count data in a matrix or data frame named 'counts'
# where rows are genes and columns are samples
# Replace 'counts' with your actual count data


source("import_data.r")
my_data = import_data("filtered_counts.txt")
# contains function/gene per row, sample per column counts

source("import_metadata.r")
my_metadata <- import_metadata("filtered_counts.metadata.txt")
selected_column = "env_package.data.body_site" # or # colnames(metadata)[1]
# contains sample per row, metadata catagory per column
condition_variable <- metadata[,selected_column]
## levels(as.factor(metadata[,selected_column]))
# samples names
names(condition_variable)
# sample metadata values
as.listr(condition_variable)




# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = my_data,
                              colData = my_metadata,
                              design = ~ env_package.data.body_site)

# Replace counts with your actual count matrix or data frame. 
# your_metadata_dataframe is your metadata containing information 
# about each sample, and condition_variable is the variable in your
# metadata indicating the different conditions or groups.

# Perform differential expression analysis
dds_processed <- DESeq(dds)

# Get differential expression results
results <- results(dds_processed)

# Summarize differential expression results
summary(results)

# From # https://bookdown.org/ggiaever/RNAseq_analysis/count-normalization.html
normalized_counts <- counts(dds_processed, normalized=TRUE)

# compare raw and normed distributions
png(
  filename = "my_boxplots.png",
  height = 8.5,
  width = 11,
  res = 300,
  units = 'in'
)
plot.new()
split.screen(c(2,1))
screen(1)
graphics::boxplot(my_data, las=2, cex.axis=0.5)
screen(2)
graphics::boxplot(normalized_counts,las=2, cex.axis=0.5)
dev.off()
#boxplot_message <- paste("output boxplot:        ", boxplots_file, "\n", sep="", collapse="")



# try it again with violin plots
# Loading the required library
library(tidyr)

# Using pivot_longer to collapse multiple columns into two
collapsed_data <- pivot_longer(data, cols = -ID, names_to = "Variable", values_to = "Value")

# reshape data to make it ready to become "tidy"
normalized_counts <- as.data.frame(normalized_counts)
old_matrix <- normalized_counts
new_column <- rownames(normalized_counts)
new_matrix <- cbind(new_column,old_matrix)
new_matrix <- as.data.frame(new_matrix)
# make data "tidy". by pivoting 
collapsed_normalized_counts <- pivot_longer(new_matrix, col=-new_column, names_to = "Sample_id", values_to = "Values" )

# Creating a violin plot
ggplot(collapsed_normalized_counts, aes(x = Sample_id, y = Values)) +
  geom_violin() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  labs(
    title = "Violin Plot Example",
    x = "Group",
    y = "Value"
  )

# Creating a boxplot
ggplot(collapsed_normalized_counts, aes(x = Sample_id, y = Values)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Violin Plot Example",
    x = "Group",
    y = "Value"
  )











# simple violin plot example
# Generating sample data
set.seed(42)
group1 <- rnorm(100, mean = 0, sd = 1)
group2 <- rnorm(100, mean = 1, sd = 1)
group3 <- rnorm(100, mean = 2, sd = 1)

# Creating a data frame
data <- data.frame(
  group = rep(c("Group 1", "Group 2", "Group 3"), each = 100),
  value = c(group1, group2, group3)
)

# Loading ggplot2 library
library(ggplot2)














screen(2)
#graphics::boxplot(normalized_counts,las=2, cex.axis=0.5)


dev.off()


