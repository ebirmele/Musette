#' Compute the Domination Relationship in a Graph, its Strongly Connected Components and Normal Forms
#'
#' (See details below about what "domination" means) 
#' Compute for every alteration A :
#'   \itemize{
#'     \item the list of alterations dominated by A "loosely" (i.e. including those dominating A reciprocally)
#'     \item the list of alterations dominating A "loosely" (i.e. including those dominated by A reciprocally)
#'     \item the list of alterations dominated by A "strictly" (i.e. excluding those dominating A)
#'     \item the list of alterations dominating A "strictly" (i.e. excluding those dominated by A)
#'     \item the strongly connected component to which A belongs 
#'     \item one "normal form", i.e. an alteration without strict dominators, which dominates A (directly or indirectly) 
#'   }
#' 
#' 
#' Given a network connecting samples and alterations, and given two groups of samples, red and blue, alteration A
#' will dominate alteration B if : 
#' \itemize{
#'  \item At least <redPercent>\% of the red samples affected by B are also affected by A, 
#'  \item At least <bluePercent>\% of the blue samples affected by A are also affected by B,
#'  \item The loci of A and B lie on the same chromosome, at most <distance> base pairs apart (this part is for alterations which are genomic, and can simply be ignored for other applications)
#' }
#' This is a somewhat more flexible version of the more natual criteria according to which A would dominate B if : 
#' \itemize{
#'  \item The red samples affected by B are also affected by A
#'  \item The blue samples affected by A are also affected by B
#' } 
#' To compute this simpler and more restrictive relation, use \code{redPercent=100} and \code{bluePercent=100}.
#' 
#' @inheritParams initSolutionTree
#' @param redPercent The required percentage of red samples
#' @param bluePercent The required percentage of blue samples
#' @param chromosome (for genomic alterations) An integer vector with the same length as \code{graph}, containing
#' the chromosome of each alteration. If omitted, this parameter has no effect.
#' @param locus (for genomic alterations) An integer vector with the same length as \code{graph}, containing
#' the position of each alteration (e.g. in base pairs).  If omitted, this parameter has no effect.
#' @param distance The maximum distance accepted between two alterations' loci for a
#' domination to occur. If omitted, this parameter has no effect.
#' @return A list \code{l} of 6 named lists of character vectors. If \code{A} is the name of an alteration, then :
#' \itemize{
#' \item \code{l$loose_dominators[[A]]} contains the names of alterations "loosely" dominating \code{A}.
#' \item \code{l$loose_dominated[[A]]} contains the names of alterations "loosely" dominated by \code{A}.
#' \item \code{l$strict_dominators[[A]]} contains the names of alterations "strictly" dominating \code{A}.
#' \item \code{l$strict_dominated[[A]]} contains the names of alterations "strictly" dominated by \code{A}.
#' \item \code{l$strongroot[[A]]} contains the name of a "reference" alteration in \code{A}'s strongly 
#' connected component (the set of alterations directly or indirectly dominated by \code{A} which also 
#' directly or indirectly dominate \code{A}). The list \code{l$strongroot} can be used to check whether
#' two alterations are in the same strongly connected component.
#' \item \code{l$strongnormal[[A]]} contains the name of a "normal form", a gene which dominates
#' \code{A} and is not dominated strictly by any other alterations.
#' } 
#' 
#' @examples
#' 
#' #Real data
#' #bipartite graph data
#' data("tcga_bladder",package="musette")
#' graph.dele=matrix2graph(matrices$dele)
#' reds= (groups == 'basal')
#' names(reds)=names(groups)
#' 
#' # chromosome and position data
#' genes = unlist(rownames(matrices$dele))
#' chromosome=chromosome[genes]
#' position=position[genes]
#' dg = domination.graph(graph.dele,reds,chromosome=chromosome,locus=position)
#' 
#' dg$strict_dominated$CDKN2A
#' dg$strict_dominators$CDKN2A
#' 
#' @export domination.graph
domination.graph=function(graph,reds,blues=!reds,chromosome=rep(1,length(graph)),locus=rep(1,length(graph)),distance=500000,redPercent=50,bluePercent=50)
{
  
  if (is.null(names(reds))) return ("names(red) must contain the names of red and blue nodes.")
  indivindex=seq_along(reds)-1
  names(indivindex)<-names(reds)
  indexedGraph=lapply(graph,function(line) indivindex[line]) 
  
  takeGraph(indexedGraph)
  takeGroups(reds,blues)
  l=rawDominationGraphFineTuned(chromosome,locus,distance,redPercent,bluePercent)
  strongcomp=scc(l$loose_dominated)
  
  l$loose_dominators=sapply(l$loose_dominators,function(line) names(graph)[line+1])
  l$loose_dominated=sapply(l$loose_dominated,function(line) names(graph)[line+1])
  l$strict_dominators=sapply(l$strict_dominators,function(line) names(graph)[line+1])
  l$strict_dominated=sapply(l$strict_dominated,function(line) names(graph)[line+1])
  names(l$loose_dominators)=names(graph)
  names(l$loose_dominated)=names(graph)
  names(l$strict_dominators)=names(graph)
  names(l$strict_dominated)=names(graph)
  
  l$strongroot=names(graph)[strongcomp[[1]]+1] # the root of the component is its defining vertex
  l$strongnormal=names(graph)[strongcomp[[2]]+1]
  
  return(l)
}
