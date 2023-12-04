#' CoefMat_Create_Spike_Slab
#'@description
#' Uses spike and slab prior to get the corresponding coefficient matrix Bij corresponding to \eqn{Y = B_{ij}X_{i} + e_{i}}. The coefficients represents the connections from the jth entry of Y to the ith entry of Y.
#' @param Y The data the causal discovery adjacency matrix wants to show the connections with.
#' @param A The adjacency matrix that shows us the structure of the DAG. The input should be a lower triangular matrix where the possible matrix entry can only be 0
#' @returns A matrix that contains all the coefficients corresponding to specific relationships between each entry of the Y data
#'
#' @examples
#' CoefMat_Create_Spike_Slab(Y_vector,A_matrix)
#' @export


CoefMat_Create_Spike_Slab = function(Y, A){
  # Does all the dimension checks where the Y should have number of rows and columns the same as number of data points of Y
  if(length(Y) != nrow(A)){
    stop('The Adjacency matrix should have number of rows that are equal to number of data in Y')
  }

  if(length(Y) != ncol(A)){
    stop('The Adjacency matrix should have number of columns that are equal to number of data in Y')
  }

  # Checks if the matrix is upper triangular
  if(sum(A[upper.tri(A)]==0)!=length(A[upper.tri(A)])){
    stop('The Adjacency matrix should be lower triangular')
  }


  A_edit = t(t(A)*Y)

  unwinded_adj_mat <- c(t(A_edit))

  # Makes sure the value corresponding to each of the 1 entries in the adjacency matrix has the corresponding value of Yj where j is the column number (This will help create the design matrix)
  design_matrix <- matrix(0, nrow = length(Y), ncol = length(unwinded_adj_mat))

  index_A = 1:nrow(design_matrix)
  index_design = 1:length(unwinded_adj_mat)

  # Creates the design matrix row that will be used to predict Yi
  for(i in 1:nrow(design_matrix)){

    test = unwinded_adj_mat
    test[index_design != index_A] = 0 # Makes sure the edges that do not connect to Yi are not used to predict Yi
    design_matrix[i,] = test

    index_A = index_A+ncol(A)
  }

  # Uses spikeslab prior on the coefficient matrix Bij to use sampling to get values of Bij
  coefficient_matrix = matrix(spikeslab(Y~design_matrix)$bma,nrow=nrow(A),ncol=ncol(A),byrow=TRUE)
  coefficient_matrix[abs(coefficient_matrix)<10^-6]=0

  return(coefficient_matrix)
}
