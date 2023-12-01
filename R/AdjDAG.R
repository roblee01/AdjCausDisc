#' AdjDAG
#'
#' @param A The adjacency matrix that shows us the structure of the DAG. The input should be a lower triangular matrix where the possible matrix entry can only be 0
#' @return The plot of the Directed Acyclic Graph that is described within the adjacency matrix
#' @examples
#' AdjDAG(A_matrix)


AdjDAG <- function(A){

  names <- paste('Y',1:nrow(A),sep="")

  DAG_Strings <- c(as.formula(paste(names[2],paste(names[A[2,]==1],collapse="+"),sep="~")))

  DAG_Strings <- c(as.formula(paste(names[2],paste(names[A[2,]==1],collapse="+"),sep="~")))
  for(i in 3:nrow(A)){
    DAG_Strings = append(DAG_Strings,c(as.formula(paste(names[i],paste(names[A[i,]==1],collapse="+"),sep="~"))))
  }

  pdf("DAG_result.pdf")


  return(ggdag(do.call(dagify, DAG_Strings)))

  dev.off()
}
