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







