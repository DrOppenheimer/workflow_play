remove_last_n_columns <- function(data,n) {
  # Check if the input is a tibble
  if (!inherits(data, "tbl_df")) {
    stop("Input should be a tibble.")
  }
  
  # Store row names
  row_names <- rownames(data)
  
  # Remove the last 7 columns
  result <- data %>%
    select(-tail(seq_along(.), n))
  
  # Restore row names
  rownames(result) <- row_names
  
  return(result)
}