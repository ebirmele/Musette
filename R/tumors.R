#' Run the Whole Musette Pipeline in One Function
#'
#' For a set of alteration matrices (for example amplifications, deletions and mutations) on the same samples partitioned into reds and blues, 
#' \itemize{
#'  \item the blind-domination and color-domination relatioships are computed using the function \link{domination.graph} 
#'  \item for each relationship type (blind or color), «leader alterations» are chosen, which dominate a nuber of «follower alterations»
#'  This information will be returned in the \code{alterations} dataframe. The procedure and goes on with the leader alterations only.
#'  \item the solution tree is constructed using the function \link{solutions}. Solutions are returned in the \code{solutions} dataframe.
#' }
#' 
#' @inheritParams initSolutionTree
#' @param matrices The alteration matrices, which columns correspond to the n samples and rows are named by the corresponding alterations 
#' @param reds A named logical vector whose names are the samples' names. The
#' values are \code{TRUE} for red samples, i.e. those which our solutions
#' should affect.
#' @param blues A similar vector specifying the blue samples, those that our
#' solutions should not affect. It is forbidden for a sample to be red and blue
#' at the same time !
#' @param blind_domination_step A character array containing the alteration types for which blind domination must be used
#' @param color_domination_step A character array containing the alteration types for which color domination must be used
#' @param blind_percent The required percentage of coverage for the blind domination step (see \link{domination.graph})
#' @param red_percent The required percentage of red samples for the domination step (see \link{domination.graph})
#' @param blue_percent The required percentage of blue samples for the domination step (see \link{domination.graph})
#' @param blind_distance The maximum distance accepted between two alterations' loci for a blind
#' domination to occur. If omitted, this parameter has no effect.
#' @param color_distance The maximum distance accepted between two alterations' loci for a blind
#' domination to occur. If omitted, this parameter has no effect.
#' @param chromosome (for genomic alterations) An integer vector with the same length as \code{graph}, containing
#' the chromosome of each alteration. If omitted, this parameter has no effect.
#' @param position (for genomic alterations) An integer vector with the same length as \code{graph}, containing
#' the 'locus' of each alteration (e.g. in base pairs).  If omitted, this parameter has no effect.
#' @param longname Character vector containing the long name or description of each gene
#' @param pathways Array containing the pathways each gene belongs to.  NOT USED (maybe in the future ?)
#' @param stepmode A string to set how the step-score is calculated. Valid values
#' are \code{"original"}, \code{"unique"}, \code{"bestFirst"}, \code{"strict"}.
#' @param threshold The maximum step-score for creating new solutions (see
#' details about the step-score in \link{musette-algorithm}).
#' @param bound The program will stop generating new solutions when this number
#' of solutions is reached. Note : \dQuote{equally good}
#' solutions in terms of the step-score are added at the same time, so that the final number of solutions possibly exceeds this
#' \code{bound}. 
#' @return A list of two objects. The \code{alterations} dataframe contains the list of considered alterations.
#'  Each row corresponds to an alteration and columns include, for each alteration:
#' \itemize{
#' \item its \code{name}, \code{type} and the concerned \code{gene}
#' \item its location coded by its  \code{chromosome} and \code{position} in terms of base-pairs
#' \item \code{longname}: its extended name
#' \item the \code{pathways} it belongs to
#' \item \code{neighbours} the samples it is present in
#' \item \code{redneighbours, blueneighbours, nredneighbours, nblueneighbours}: the set of red (resp. blue) samples it is present in, and the cardinal of those sets
#' \item the hypergeometric score \code{hyper} and the musette score \code{score} corresponding to an alteration set reduced to this single alteration
#' \item \code{blind_dominators}, \code{blind_dominated}, \code{color_dominators}, \code{color_dominated}: sets of alteration that dominate or are dominated by this alteration in the blind or colored way
#' \item \code{blind_leader}, \code{color_leader}: alteration chosen to lead the group with dominations (it may be not be a direct dominator of the considered alteration but dominate a dominator) 
#' \item \code{blind_repr} (\code{color_repr}): a representative chosen in this alteration's strongly connected component of the blind domination (color domination) directed graph. All the alterations in a given strongly connected component will have the same representative
#' \item \code{followers}: the set of all alterations whose color_leader is this particular alteration
#' \item \code{merged_in}: pseudo-alteration in which a real alteration is potentially merged 
#' \item \code{merged_from}: real alterations merged to form a pseudo-alteration
#' \item \code{stepscore}: stepscore of this alteration seen as an evolution of the empty set.
#' }
#
#'  The \code{solutions} data.frame contains a row for each generated solution, in order of decreasing score.
#' Columns include :
#' \itemize{
#' \item alterations of each solution (a character vector)
#' \item size : the number of alterations
#' \item number of red samples hits
#' \item number of blue samples hits
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
#' \item all_genes : list of all the alterations belonging to this solution or
#'  dominated, directly or indirectly by such an alteration
#' \item all_loci : list of all chromosomes of the genes in this solution
#' } 
#' @examples 
#' #Real bladder data, restricted on the 20 first solution sets including only mutations
#' data("tcga_bladder",package="musette")   
#' reds= (groups == 'basal')
#' names(reds)=names(groups)
#' 
#' ll=do.musette(matrices=matrices, reds=reds, blind_domination_step = c("ampli","dele"), 
#' color_domination_step = c("ampli","dele"), blind_distance=5000000, 
#' color_distance=5000000, blind_percent=90, red_percent=80, 
#' blue_percent=80, chromosome = chromosome, longname=longname,
#' pathways=pathways,position=position, bound=20)
#' 
#' View(ll$solutions)
#' View(ll$alterations)
#' 
#' @export do.musette

do.musette=function(matrices,reds, blues=!reds, blind_domination_step=c(),color_domination_step=c(),
                    blind_distance=5e5,blind_percent=80,color_distance=5e5,red_percent=80,blue_percent=80,
                    chromosome=NULL,position=NULL, longname=NULL,pathways=NULL,bound=2000, stepmode="strict",threshold=1){
  # préparation des ingrédients du dataframe des altérations
  type=array(unlist( sapply(names(matrices),function(name) rep(name,nrow(matrices[[name]])))))
  gene=array(unlist(lapply(matrices,function(m) rownames(m))))
  name=array(paste(type,"_",gene,sep=""))
  fullmatrix=do.call(rbind,matrices)
  rownames(fullmatrix)=name
  
  print("computing full graph...")
  neighbours=sapply(name,function(x) names(which(fullmatrix[x,]==1)),simplify=FALSE )
  print("done.")
  
  #création du dataframe des altérations
  alterations=data.frame(name=name,type=type,gene=gene,
                         chromosome=chromosome[gene], position=position[gene], longname=longname[gene], stringsAsFactors = FALSE)
  alterations$pathways=pathways[gene]
  rownames(alterations)=name
  
  alterations$neighbours=neighbours # was previously : sapply(neighbours,unlist)
  
  alterations$redneighbours=sapply(alterations$neighbours, function(l) l[reds[l]])
  alterations$blueneighbours=sapply(alterations$neighbours, function(l) l[blues[l]])
  
  alterations$nbredneighbours=lengths(alterations$redneighbours)
  alterations$nbblueneighbours=lengths(alterations$blueneighbours)
  
  alterations$hyper=-phyper(alterations$nbblueneighbours,sum(blues),sum(reds)
			    ,alterations$nbredneighbours+alterations$nbblueneighbours,log.p=TRUE)
  alterations$score=alterations$hyper*alterations$nbredneighbours/sum(reds)
  alterations$stepscore=sapply(seq(nrow(alterations)),function(i) 
			       singleStepScore(0,0,alterations$nbredneighbours[i],alterations$nbblueneighbours[i])
			      )
  
  # Domination considering only the number of neighbours, ignoring red/blue
  alterations$blind_dominators=sapply(rownames(alterations), function(x) character())
  alterations$blind_dominated=sapply(rownames(alterations), function(x) character())
  alterations$blind_leader=alterations$name
  alterations$blind_repr=alterations$name
  names(alterations$blind_leader)=alterations$name
  
  for(alt_type in blind_domination_step){
    partial_graph=alterations$neighbours[alterations$type==alt_type]
    names(partial_graph)=alterations$name[alterations$type==alt_type]
    partial_chromosome=alterations$chromosome[alterations$type==alt_type]
    partial_position=alterations$position[alterations$type==alt_type]
    reds_all=sapply(reds,function(x) TRUE,USE.NAMES = TRUE)
    
    print(paste("computing blind domination for ", alt_type))
    dom_result=domination.graph(partial_graph,reds_all,chromosome=partial_chromosome,locus=partial_position,
                                redPercent=blind_percent,bluePercent=0,distance=blind_distance)
    print("done.")
    # get the list of members of each scc
    l_scc=split(names(partial_graph),dom_result$strongroot)[unique(dom_result$strongroot)]
    # N.B. :l_scc has names, which will be passed on to scc_leaders and used when filling columns "blind_leader" and "blind_repr"
    
    
    alterations$blind_dominators[names(dom_result$loose_dominators)]=dom_result$loose_dominators
    alterations$blind_dominated[names(dom_result$loose_dominated)]=dom_result$loose_dominated
    
    print(paste("computing blind leaders for",alt_type))
    scc_leader=lapply(l_scc,function(scc){
      len=lengths(alterations[scc,"neighbours"])
      champions=which(len== max(len)) # champions' indices 
      positions=alterations[scc[champions],"position"]
      scc[[champions[[which.min(abs(positions-mean(positions)))]]]] # winner = closest to middle
    })
    print("done.")
    alterations[names(partial_graph),"blind_leader"]=as.character(scc_leader[dom_result$strongnormal])
    alterations[names(partial_graph),"blind_repr"]=as.character(scc_leader[dom_result$strongroot])
  }
  
  # Now domination taking red/blue samples into account
  alterations$color_dominated=sapply(rownames(alterations), function(x) character())
  alterations$color_dominators=sapply(rownames(alterations), function(x) character())
  alterations$color_leader=alterations$name
  alterations$color_repr=alterations$name
  
  for(alt_type in color_domination_step){
    partial_graph=alterations$neighbours[alterations$type==alt_type & alterations$blind_leader==alterations$name]
    names(partial_graph)=alterations$name[alterations$type==alt_type & alterations$blind_leader==alterations$name]
    partial_chromosome=alterations$chromosome[alterations$type==alt_type & alterations$blind_leader==alterations$name]
    partial_position=alterations$position[alterations$type==alt_type & alterations$blind_leader==alterations$name]
    
    print(paste("computing color-aware domination for ", alt_type))
    dom_result=domination.graph(partial_graph,reds,blues,chromosome=partial_chromosome,locus=partial_position,redPercent=red_percent,bluePercent=blue_percent,distance=color_distance)
    print("done.")
    l_scc=split(names(partial_graph),dom_result$strongroot)[unique(dom_result$strongroot)]
    
    alterations$color_dominators[names(dom_result$loose_dominators)]=dom_result$loose_dominators
    alterations$color_dominated[names(dom_result$loose_dominated)]=dom_result$loose_dominated
    
    print(paste("computing color-aware leaders for",alt_type))
    scc_leader=lapply(l_scc,function(scc){
      len=alterations[scc,"score"]
      champions=which(len== max(len)) # champions : indices dans len
      positions=alterations[scc[champions],"position"]
      scc[[champions[[which.min(abs(positions-mean(positions)))]]]]
    })
    print("done.")
    alterations[names(partial_graph),"color_leader"]=as.character(scc_leader[dom_result$strongnormal])
    alterations[names(partial_graph),"color_repr"]=as.character(scc_leader[dom_result$strongroot])
  }
  alterations[,"color_leader"]=alterations[as.character(alterations$blind_leader),"color_leader"]
  alterations[,"color_repr"]=alterations[as.character(alterations$blind_leader),"color_repr"]
  
  alterations$followers=sapply(rownames(alterations), function(x) character())
  followers=split(alterations$name,alterations$color_leader)
  alterations$followers[names(followers)]=followers
  
  
  # Are there alterations which affect exactly the same samples ? Merge them into pseudo-alterations.
  # The final decision was to not split them afterwards.
  partial_graph=alterations$neighbours[lengths(alterations$followers)>0]
  canonic=match(partial_graph,partial_graph)# or use S4Vectors::selfmatch
  family_indices<-split(seq_along(partial_graph),canonic) # reverse : whose canonic is this ? 
  family_indices=family_indices[lengths(family_indices)>1]
  
  pseudo_alterations=data.frame(name=as.character(paste("pseudo",names(family_indices),sep="_")),stringsAsFactors = FALSE)
  rownames(pseudo_alterations)=pseudo_alterations$name
  
  print("computing attributes for pseudo_alterations...")
  for(field in setdiff(names(alterations),"name")){
    pseudo_alterations[[field]]=sapply(names(family_indices),function(j){
      t=unlist(alterations[names(partial_graph)[family_indices[[j]]],field])
      if(is.null(t)) character(0) else sort(unique(t)) # there WERE cases when t was null, but which ones ?
    })
  }
  print("done.")
  
  pseudo_alterations$merged_from=sapply(family_indices, function(l) names(partial_graph)[l])
  pseudo_alterations$merged_in=sapply(family_indices,function(x) "")
  alterations$merged_from=sapply(alterations$name,function(x) character(0))
  alterations$merged_in=sapply(alterations$name,function(x) "")
  for (i in seq_along(family_indices)){
    alterations[pseudo_alterations$merged_from[[i]],"merged_in"]=pseudo_alterations$name[[i]]
  }
  pseudo_alterations$type="pseudo"
  rownames(pseudo_alterations)=pseudo_alterations$name
  
  alterations=rbind(alterations,pseudo_alterations)
  
  eligible=alterations$name[(lengths(alterations$followers)>0 & alterations$merged_in=="") |
                              alterations$type=="pseudo"]
  graph=alterations[eligible,"neighbours"]
  names(graph)=eligible
  
  sol=solutions(graph,reds,blues=blues,bound=bound,stepmode=stepmode,threshold=threshold)
  sol$all_genes=sapply(sol$alterations,function(l){
    unique(do.call("c",rbind(l,as.list(alterations[l,"followers"])))) # to "print" each leader before its followers
  })
  sol$all_loci=sapply(sol$all_genes,function(l){
    unique(do.call("c",as.list(alterations[l,"chromosome"])))
  })
  sol$allred=sum(reds)
  sol$allblue=sum(blues)
  
  
  
  ### final result
  list(alterations=alterations,solutions=sol)

}


#' Analysis of the solutions in terms of pathways
#' 
#' Given the solutions and alterations 
#' dataframes returned by do.musette, it returns a new dataframe with solutions
#' reorganized by pathways shared by several genes in a common solution
#' 
#' @param solutions  a solution dataframe returned by do.musette 
#' @param alterations an alteration dataframe returned by do.musette
#' 
#' @return The returned data.frame contains a row for each couple solution/pathway such
#' that at least two alterations do concern genes belonging to the pathway.
#' Columns include :
#' \itemize{
#' \item sol: identification number of the solution in the \code{solutions} dataframe
#' \item alterations describing of the solution 
#' \item pathway
#' \item gene: alterations of the solution (including those dominated by the descriptors) that belong to the pathway
#' \item leader: leaders of the alterations in the gene list
#' \item red neighbors, blue neighbors of the alterations in the gene list, and the cardinal of those sets
#' } 
#' @examples 
#' #Real bladder data, restricted on the 20 first solution sets including only mutations
#' data("tcga_bladder",package="musette")   
#' reds= (groups == 'basal')
#' names(reds)=names(groups)
#' 
#' ll=do.musette(matrices=matrices, reds=reds, blind_domination_step = c("ampli","dele"), 
#' color_domination_step = c("ampli","dele"), blind_distance=5000000, 
#' color_distance=5000000, blind_percent=90, red_percent=80, 
#' blue_percent=80, chromosome = chromosome, longname=longname,
#' pathways=pathways,position=position, bound=20)
#' 
#' spdf = sharedPathways(ll$solutions,ll$alterations)
#' View(spdf)
#' 
#' 
#' @export sharedPathways
sharedPathways=function(sol, alterations){
  
  alterations$follower_pathways=sapply(rownames(alterations),function(x) character(0))
  solgene=as.array(unique(unlist(sol$alterations)))
  alterations$follower_pathways[solgene]=lapply(solgene,
                                                function(x) {unique(do.call(c,alterations[alterations[[x,"followers"]],"pathways"]))})
  
  print("computing shared pathway informations...")
  sol$shared_pathways=sapply(sol$alterations,function(x){ 
    t=unlist(alterations[x,"follower_pathways"]) ; unique( t[duplicated(t)])
  })
  
  
  spdf=data.frame(sol=do.call(c, sapply(seq(nrow(sol)),function(i) rep(rownames(sol)[[i]],length(sol$shared_pathways[[i]])))),
                pathway=unlist(sol$shared_pathways),stringsAsFactors = FALSE)
  spdf$alterations=sol[spdf$sol,"alterations"]
  spdf=data.frame(spdf[1],spdf[3],spdf[2]) # just reordering columns

  spdf$gene=lapply(seq(nrow(spdf)), function(i) {x=sol[[spdf$sol[[i]],"all_genes"]]; y=sapply(x,function(g) spdf$pathway[[i]] %in% alterations[[g,"pathways"]]) ; x[y] })
  spdf$leader=lapply(seq(nrow(spdf)), function(i) {x=sol[[spdf$sol[[i]],"alterations"]]; y=sapply(x,function(g) spdf$pathway[[i]] %in% alterations[[g,"follower_pathways"]]) ; x[y] })
  spdf$redneighbours=lapply(spdf$gene, function(x) unique(unlist(alterations[x,"redneighbours"])))
  spdf$blueneighbours=lapply(spdf$gene, function(x) unique(unlist(alterations[x,"blueneighbours"])))
  spdf$nbredneighbours=lengths(spdf$redneighbours)
  spdf$nbblueneighbours=lengths(spdf$blueneighbours)

  return(spdf)
}


