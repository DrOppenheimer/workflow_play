source("~/Documents/GitHub/Kevin_R_scripts/import_data.r")
# my_data <- import_data("~/Downloads/GPL1708_family.soft") # DOES NOT WORK
library(tidyverse)
my_data <- read_csv("~/Downloads/GPL1708_family.soft") # Try this

source("~/Documents/GitHub/Kevin_R_scripts/import_metadata.r")
my_metadata <- import_metadata("~/Downloads/GPL1708-20418.txt")
# Fails like this
# Error in vroom_(file, delim = delim %||% col_types$delim, col_names = col_names,  :                                          
#                  R character strings are limited to 2^31-1 bytes

BiocManager::install("Seurat")
library(Seurat)
browseVignettes(package = "Seurat")