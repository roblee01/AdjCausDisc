#' CoefMat_Create_Spike_Slab
#'
#' @param Y The data the causal discovery adjacency matrix wants to show the connections with.
#' @param A The adjacency matrix that shows us the structure of the DAG. The input should be a lower triangular matrix where the possible matrix entry can only be 0
#' @returns A matrix that contains all the coefficients corresponding to specific relationships between each entry of the Y data
#'
#' @examples
#' CoefMat_Create_Spike_Slab(Y_vector,A_matrix)


CoefMat_Create_Spike_Slab = function(Y, A){
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

  coeffecient_matrix = matrix(spikeslab(Y~design_matrix)$bma,nrow=nrow(A),ncol=ncol(A),byrow=TRUE)
  coeffecient_matrix[abs(coeffecient_matrix)<10^-6]=0

  return(coeffecient_matrix)
}
