#' CoefMat_Create_SSVS
#'
#' @param Y The data the causal discovery adjacency matrix wants to show the connections with.
#' @param A The adjacency matrix that shows us the structure of the DAG. The input should be a lower triangular matrix where the possible matrix entry can only be 0
#' @returns A matrix that contains all the coefficients corresponding to specific relationships between each entry of the Y data using Stochastic Search Variable Search to get the values of the coeffecients.
#'
#' @examples
#' CoefMat_Create_SSVS(Y_vector,A_matrix)


CoefMat_Create_SSVS = function(Y, A){
  # given in the Yi = ∑ on (where the index j of Yj corresponding to parent nodes of Yi) of (bij)*(Yj) + ei. The coeffecient matrix will have all the possible bij results that is calculated using spike and slab prior.

  if(length(Y) != nrow(A)){
    stop('The Adjacency matrix should have rows that are equal to number of data in Y')
  }

  if(sum(A[upper.tri(A)]==0)!=length(A[upper.tri(A)])){
    stop('The Adjacency matrix should be lower triangular')
  }


  A_edit = t(t(A)*Y)

  unwinded_adj_mat <- c(t(A_edit))

  design_matrix <- matrix(0, nrow = length(Y), ncol = length(unwinded_adj_mat))

  index_A = 1:nrow(design_matrix)
  index_design = 1:length(unwinded_adj_mat)
  for(i in 1:nrow(design_matrix)){

    test = unwinded_adj_mat
    test[index_design != index_A] = 0
    design_matrix[i,] = test

    index_A = index_A+ncol(A)
  }

  col_names_design<-paste(as.character(1),1:length(Y),sep="")
  for(i in 2:length(Y)){
    col_names_design = append(col_names_design,paste(as.character(i),1:length(Y),sep=""))
  }

  colnames(design_matrix) = col_names_design

  X<-design_matrix[,colSums(design_matrix)!=0]

  big_dat = as.data.frame(cbind(Y,X))
  corresponding_Y_names = colnames(X)
  results_ssvs = setNames(summary(ssvs(big_dat,colnames(big_dat)[1],colnames(big_dat)[2:length(colnames(big_dat))]))$`Avg Beta`,corresponding_Y_names)

  indices = strsplit(corresponding_Y_names,"")
  coefficient_matrix = matrix(0, nrow = nrow(A), ncol = ncol(A))
  for(i in 1:length(indices)){
    indices_2 = as.numeric(indices[[i]])
    coefficient_matrix[indices_2[1],indices_2[2]] = results_ssvs[i]
  }

  return(coeffecient_matrix)
}
