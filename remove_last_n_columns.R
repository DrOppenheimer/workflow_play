remove_last_n_columns <- function(data,n) {

  # # Store row names
  row_names <- rownames(data)
  
  # Removing the last n columns
  data <- data[, -seq(ncol(data) - n + 1, ncol(data))]
  
  # # Restore row names
  rownames(data) <- row_names
  
  # View the modified matrix
  return(data)
  
}








# # Check if the input is a tibble
# if (!inherits(data, "tbl_df")) {
#   stop("Input should be a tibble.")
# }
# 
# # Store row names
# row_names <- rownames(data)
# 
# # Remove the last 7 columns
# result <- data %>%
#   select(-tail(seq_along(.), n))
# 
# # Restore row names
# rownames(result) <- row_names
# 
# return(result)