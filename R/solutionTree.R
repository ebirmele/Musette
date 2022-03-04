
#' Initialize a Solution Tree
#'
#' Use the provided arguments to initialize data structures that
#' will contain, and help to build, a solution tree.
#'
#' (See \link{musette-algorithm} to learn what "solution tree" means)
#' This is an internal function meant for use by other functions in the package. 
#' An expert user wanting to be in full control could use it, followed 
#' by calls to \link{growTree}, \link{fullyGrowTree}, or \link{exploreTree}.
#' @param graph A named list of character vectors. Each of these vectors 
#' should have the name of one alteration, and contain the vector of the 
#' names of the samples affected by this alteration. Such a list can be 
#' obtained from a logical matrix using \code{\link{matrix2graph}}
#' @param reds A named logical vector whose names are the samples' names. The
#' values are \code{TRUE} for red samples, i.e. those which our solutions
#' should affect.
#' @param blues A similar vector specifying the blue samples, those that our
#' solutions should not affect. It is forbidden for a sample to be red and blue
#' at the same time !
#' @param threshold The maximum step-score for creating new solutions (see
#' details about the step-score in \link{musette-algorithm}).
#' @param bound The program will stop generating new solutions when this number
#' of solutions is reached. Note : if there are many \dQuote{equally good}
#' solutions they are added at the same time, possibly exceeding this
#' \code{bound}, then the algorithm stops.
#' @param stepmode A string to set how the step-score is calculated. Valid values
#' are \code{"original"}, \code{"unique"}, \code{"bestFirst"}, \code{"strict"}.
#'
#' @export
initSolutionTree=function(graph,reds,blues=!reds,threshold=1,bound=2000,stepmode="strict")
{
  if (is.null(names(reds))) stop ("names(red) must contain the names of red and blue nodes.")
  if(any(reds & blues) ) stop ( "It s forbidden for a sample to be both red and blue.")
  indivindex=seq_along(reds)-1
  names(indivindex)<-names(reds)
  indexedgraph=lapply(graph,function(line) indivindex[line]) 
  
  takeStepScore(list(original=0,bestFirst=1,strict=2)[[stepmode]])
  takeGraph(indexedgraph)
  setBound(bound)
  setThreshold(threshold)
  takeGroups(reds,blues)
  resetSolutionTree()
}

#' Monitor the Growth of the Solution Tree with Respect to the Threshold
#' 
#' Grow the solution tree completely, only to see how its size evolves
#' in relation to the threshold used.
#' 
#' See \link{musette-algorithm} to learn about the algorithms and the its parameters.
#' 
#' @inheritParams initSolutionTree
#' @return A dataframe with 3 columns showing the evolution of the threshold,
#' the number of nodes and the number of leaf nodes, during the growth of the
#' solution tree. Its study can give an insight as to what could be a good
#' threshold value to use.
#' @examples 
#' #Toy example
#' toy.matrix <- matrix(c(1,0,1,0,0,0,1,1,0,0,1,0,0,0,1,1,0,1,0,1,0,0,0,1),4,6,byrow=TRUE)
#' rownames(toy.matrix) <- LETTERS[1:4]
#' colnames(toy.matrix) <- letters[1:6]
#' toy.graph <- matrix2graph(toy.matrix)
#' reds <- c(rep(TRUE,3),rep(FALSE,3))
#' names(reds) <- colnames(toy.matrix)
#' treeGrowth(toy.graph,reds)
#' 
#' #Real bladder data, restricted on the 100 first solution sets including only mutations
#' data("tcga_bladder",package="musette")   
#' reds= (groups == 'basal')
#' names(reds)=names(groups)
#' graph.muta <- matrix2graph(matrices$muta)
#' tg <- treeGrowth(graph.muta,reds,bound=100)
#' 
#' @export
treeGrowth=function(graph,reds,blues=!reds,threshold=1,bound=2000,stepmode="strict")
{
  initSolutionTree(graph,reds,blues,threshold,bound,stepmode)
  return(data.frame(exploreTree()))
}






#' Enumerate Sets of Alterations Specific to a Given Group of Samples
#' 
#' See \link{musette-algorithm} for a detailed explanation 
#' of what is meant by \dQuote{linked}, and the algorithm used.
#' 
#' Two halting criteria are used : \code{bound} and \code{threshold}. Their
#' default values are high enough, so providing a sensible value for only one
#' of them effectively results in ignoring the other.
#' 
#' @inheritParams initSolutionTree
#' @return A data.frame with a row for each generated solution, in order of
#' decreasing score.
#' Columns include :
#' \itemize{
#' \item alterations of each solution (a character vector)
#' \item size : the number of alterations
#' \item number of red samples hit
#' \item number of blue samples hit
#' \item sensitivity
#' \item specificity
#' \item total number of red samples
#' \item total number of blue samples
#' \item step-score of adding the last element of this solution
#' \item threshold : maximum step-score of all the steps in constructing this 
#' solution
#' \item leaf : TRUE indicates that this solution has not been extended into a larger solution; FALSE indicates that the dataframe contains at least an extended solution.
#' \item children-threshold : minimum step-score of extending this solution into
#'  a new solution
#' \item parent : row number of the solution that was extended into this one
#' \item score : score of this solution.
#' }
#' @examples
#' #Toy example
#' toy.matrix <- matrix(c(1,0,1,0,0,0,1,1,0,0,1,0,0,0,1,1,0,1,0,1,0,0,0,1),4,6,byrow=TRUE)
#' rownames(toy.matrix) <- LETTERS[1:4]
#' colnames(toy.matrix) <- letters[1:6]
#' toy.graph <- matrix2graph(toy.matrix)
#' reds <- c(rep(TRUE,3),rep(FALSE,3))
#' names(reds) <- colnames(toy.matrix)
#' sol <- solutions(toy.graph,reds)
#' View(sol)
#' 
#' #Real bladder data, restricted on the 100 first solution sets including only mutations
#' data("tcga_bladder",package="musette")   
#' reds= (groups == 'basal')
#' names(reds)=names(groups)
#' graph.muta <- matrix2graph(matrices$muta)
#' sol <- solutions(graph.muta,reds,bound=100)
#' View(sol)
#' 
#' @export
solutions=function(graph,reds,blues=!reds,threshold=1,bound=2000,stepmode="strict")
{
  initSolutionTree(graph,reds,blues,threshold,bound,stepmode)
  fullyGrowTree()
  res=encodeSolutionTree()
  res=append(res,list(size=sapply(res$alterations,length),red.total=sum(reds),blue.total=sum(blues),
                      sensitivity=res$red/sum(reds),specificity=1-res$blue/sum(blues)),6)
  
  sol.df=data.frame(res[2:length(res)],stringsAsFactors = FALSE)
  sol.df$alterations=sapply(res[[1]], function(indexset) names(graph)[indexset+1])
  sol.df=sol.df[c(length(sol.df),1:(length(sol.df)-1))]
  rownames(sol.df)=as.character(1:nrow(sol.df)) # the solution's name (number) is fixed already, this value is used in column 'parent'.
  sol.df=sol.df[order(sol.df$score,decreasing = TRUE),]
  return(sol.df)
}

