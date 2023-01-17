# Function to convert an upper triangular matrix to a full correlation matrix
upper_triangular_to_cor_matrix <- function(triangular_mat) {
  # Get the dimension of the matrix
  dim_mat <- dim(triangular_mat)[1]
  # Create an empty matrix
  full_mat <- matrix(0, nrow = dim_mat, ncol = dim_mat)
  # Fill in the upper triangular values
  full_mat[upper.tri(full_mat, diag = TRUE)] <- triangular_mat
  # Copy the upper triangular values to the lower triangular part
  full_mat <- t(full_mat)
  full_mat[upper.tri(full_mat, diag = TRUE)] <- triangular_mat
  return(full_mat)
}

# Test the function
mat_data <- c(1,2,3,NA,5,6,NA,NA,9)
mat_data
triangular_mat <- matrix(mat_data,nrow=3,ncol=3,byrow=TRUE)
triangular_mat
dim(triangular_mat)[1]
upper_triangular_to_cor_matrix(triangular_mat)

# Function to normal score transform columns of a dataframe
# will do select columns pass in as 'columns'
normal_score_transform <- function(df, columns) {
  for (col in columns) {
    col_data <- df[,col]
    # Compute the normal scores
    normal_scores <- qnorm(rank(col_data)/(length(col_data)+1))
    # Replace the original column data with the normal scores
    df[,col] <- normal_scores
  }
  return(df)
}

# Test the function
df <- data.frame(x1 = c(1, 2, 3, 4, 5), x2 = c(6, 7, 8, 9, 10))
normal_score_transform(df, c("x1", "x2"))
normal_score_transform(df, 1:2) # can take a vector of column numbers
