#' AdjDAG
#'
#' @param A The adjacency matrix that shows us the structure of the DAG. The input should be a lower triangular matrix where the possible matrix entry can only be 0
#' @return The plot of the Directed Acyclic Graph that is described within the adjacency matrix
#' @examples
#' Adjacency_matrix = matrix(nrow=3,ncol=3,0)
#' Adjacency_matrix[,1]=c(0,1,1)
#' Adjacency_matrix[,2]=c(0,0,1)
#' Adjacency_matrix[,3]=c(0,0,0)
#' AdjDAG(Adjacency_matrix)
#' @export


AdjDAG <- function(A){

  names <- paste('Y',1:nrow(A),sep="")

  # Creating all the DAG strings that will contain all the information and different edges present in the DAG
  DAG_Strings <- c(as.formula(paste(names[2],paste(names[A[2,]==1],collapse="+"),sep="~")))
  for(i in 3:nrow(A)){
    DAG_Strings = append(DAG_Strings,c(as.formula(paste(names[i],paste(names[A[i,]==1],collapse="+"),sep="~"))))
  }

  # Prints out the plot
  return(ggdag(do.call(dagify, DAG_Strings)))
}
