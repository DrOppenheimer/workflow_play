# Install BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# Use bioconductor to install DESeq2 package
BiocManager::install("DESeq")
# FAIL

# download the source for the package from here
# https://www.bioconductor.org/packages//2.12/bioc/html/DESeq.html
# install.packages(path_to_file, repos = NULL, type="source")
# FAILED
# got newer version from here # WORKED
# https://bioconductor.riken.jp/packages/3.8/bioc/html/DESeq.html
install.packages("DESeq/", repos = NULL, type="source")

# then as usual
library(DESeq)
# THEN FOUND THIS

# From # https://support.bioconductor.org/p/129152/
#Yes, DESeq is depreacated since at least five years, Maybe the message that we 
#put there five years ago 
#("Welcome to 'DESeq'. For improved performance, usability and functionality, 
#  please consider migrating to 'DESeq2'.") is too mild.