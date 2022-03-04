#' Make an Adjacency List out of a Logical Matrix
#' 
#' Produces an adjacency list of the type that can be passed to various
#' functions of package \link{musette}
#' 
#' 
#' @param matrix A logical matrix, with mandatory row and column names ! Row
#' names must correspond to alterations, and column names must correspond to
#' samples.
#' @return An named list of character vectors. Its names are the alterations'
#' names, and each vector contains the names of the samples affected by the
#' corresponding alteration.
#' @examples
#' #Toy example
#' toy.matrix <- matrix(c(1,0,1,0,0,0,1,1,0,0,1,0,0,0,1,1,0,1,0,1,0,0,0,1),4,6,byrow=TRUE)
#' rownames(toy.matrix) <- LETTERS[1:4]
#' colnames(toy.matrix) <- letters[1:6]
#' toy.graph <- matrix2graph(toy.matrix)
#' 
#' #Real data
#' data("tcga_bladder",package="musette")
#' graph.dele=matrix2graph(matrices$dele)
#' 
#' 
#' @export matrix2graph
matrix2graph<- function(matrix)
  sapply(rownames(matrix),function(x) names(which(matrix[x,]==1)),simplify=FALSE )





#' Turn  Vectors in a Dataframe into Strings for Better Rendering in a CSV File
#'
#' Cells of the solution dataframe containing vectors are uneasy to read. This function solves the problem when exporting the dataframe to a csv file.
#'
#' @param d A dataframe
#' @return The same dataframe, where cells that are vectors have been replaced by semicolon-punctuated strings.
#' @examples
#' dataframe=data.frame(person=c("sam","ben","alice"))
#' dataframe$food=list(c("cheese","egg","bacon"),c("tomato"),c("bread","ham"))
#' \dontrun{write.csv(dataframe,"file.csv")} # causes an error ("unimplemented type")
#' d2=musette2csv(dataframe)
#' \dontrun{write.csv(d2,"file.csv")} # no error
#' 
#' @export
musette2csv=function(d){
  d2=data.frame(lapply(d,function(c) sapply(c,paste,collapse=";") ) )
  rownames(d2)=rownames(d)
  d2
}




