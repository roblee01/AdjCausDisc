#' CoefMat_Create_SSVS
#' @description
#' Uses ssvs bayesian model selection to get the corresponding coefficient matrix Bij corresponding to \eqn{Y = B_{ij}X_{i} + e_{i}}. The coefficients represents the connections from the jth entry of Y to the ith entry of Y.
#' @param Y The data the causal discovery adjacency matrix wants to show the connections with. Y1 = c(y1, ..., yn)
#' @param A The adjacency matrix that shows us the structure of the DAG. The input should be a lower triangular matrix where the possible matrix entry and have the number of rows and number of columns same as the number of data points in Y
#' @returns A matrix that contains all the coefficients corresponding to specific relationships between each entry of the Y data using Stochastic Search Variable Search to get the values of the coeffecients.
#'
#' @examples
#' CoefMat_Create_SSVS(Y_vector,A_matrix)
#'
#'
#' @export


CoefMat_Create_SSVS = function(Y, A){
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


  # Makes sure the value corresponding to each of the 1 entries in the adjacency matrix has the corresponding value of Yj where j is the column number (This will help create the design matrix)
  A_edit = t(t(A)*Y)

  unwinded_adj_mat <- c(t(A_edit))

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

  col_names_design<-paste(as.character(1),1:length(Y),sep="")
  for(i in 2:length(Y)){
    col_names_design = append(col_names_design,paste(as.character(i),1:length(Y),sep=""))
  }

  colnames(design_matrix) = col_names_design

  X<-design_matrix[,colSums(design_matrix)!=0] #SSVS requires the input design matrix to be invertible, so we remove the columns that are of 0's

  big_dat = as.data.frame(cbind(Y,X))
  corresponding_Y_names = colnames(X)
  results_ssvs = setNames(summary(ssvs(big_dat,colnames(big_dat)[1],colnames(big_dat)[2:length(colnames(big_dat))]))$`Avg Beta`,corresponding_Y_names)

  indices = strsplit(corresponding_Y_names,"")
  coefficient_matrix = matrix(0, nrow = nrow(A), ncol = ncol(A))
  # Creates coefficient matrix and makes sure corresponding coefficients match the edge it corresponds to.
  for(i in 1:length(indices)){
    indices_2 = as.numeric(indices[[i]])
    coefficient_matrix[indices_2[1],indices_2[2]] = results_ssvs[i]
  }

  return(coefficient_matrix)
}
