
#' 
#' Find sets of alterations linked to a user-defined tumor subgroup. 
#' 
#' The user chooses among the data a subset of tumors which are called 'red'
#' and another, disjoint subset of 'blue' tumors. A good solution S is a set
#' of alterations (A1,A2,...An) such that : 
#' \enumerate{
#' \item  many of the red tumors are affected by at least one of A1,A2,An. 
#' \item  few blue tumors are affected by any of A1,A2, An. 
#' \item  the number n of alterations in S is kept small.}
#' A score is defined for any set S of alterations and the best solutions are
#' generated.
#' 
#' 
#' \bold{Data}
#' A data frame is constructed containing information about the alterations. 
#' Some are mandatory (name, number of blue/red neighbours), others are optional (alteration type,
#' gene, chromosome location, etc).
#'
#' \bold{Score definition}
#' The default score of a solution is the opposite logarithm of the hypergeometric
#' p-value, times the fraction f of red tumors (0<=f<=1) affected by any
#' alteration in S. This allows to give a preference to the alteration sets yielding 
#' a good coverage of the red tumors, compared to the hypergeometric score alone.
#'
#' \bold{Preprocessing}
#' For some alteration types (e.g. deletions and amplifications)
#' which often act on a whole chromosome segment, we define two notions of
#' domination in order to concentrate the whole segment into one alteration.
#' 
#' Consider two alterations A and B which are close enough on the genome, and A such that concerns more
#' tumors than B, including a sufficient percentage (blind_percent) of the
#' tumors touched by B. Then A is said to 'dominate' B in the color-blind
#' fashion.
#'
#' In a preprocessing phase, blind domination is computed for alterations of
#' the same type, for the types listed in 'blind_domination_step'. The set of considered alterations 
#' is reduced to a set of blind-leader alterations, that is a set such each alteration A is dominated by a 
#' some alteration, which is dominated by some alteration, ... in a chain 
#' leading to some blind-leader alteration. Only the leaders are kept to run the main algorithm.
#' 
#' If A and B are close enough, A has a better score than B, touching a
#' sufficient fraction of red tumors touched by B, and if a suffcient fraction
#' of blue tumors touched by A are also touched by B, then A dominates B in the
#' color-aware fashion.
#' The number of considered alterations is again reduced, in a similar manner to the blind domination,
#' to a set of color_leaders. Only those are kept for the
#' remaining steps of the analysis.
#' 
#' In order to be able to get back to the original alterations, a list is created for every color-leader gene A, containing the list
#' of the alterations which color_leader is A.
#' 
#'
#' \bold{Enumeration algorithm}
#' A tree of solutions (sets of alterations) is generated from a root which is
#' the empty set. Every other solution S is the child of a solution S'
#' having one less alteration. The child S is constructed and generated only if
#' its score score(S) is significantly better than its parents'. This is
#' evaluated by computing a 'step score', which can be defined as : -
#' 'original' mode : the probability of getting a better score than score(S) by
#' adding a random alteration to S'. - 'best-first' mode : the probability of
#' getting a better score than score(S) by replacing the 'worst' alteration in
#' S by a random alteration (here 'worst' is in terms of the score of the
#' individual alteration) - 'strict' mode : the highest of the p_i, where p_i
#' is the probability of getting a better score than score(S) by replacing
#' alteration number i in S by a random alteration. The solution S is added to
#' the tree only if this step-score is below a certain threshold.  The tree of solutions is grown by gradually
#' raising this threshold.
#' 
#' \bold{Outcome}
#' The solutions are displayed in a data frame which columns include the alterations included in the solutions,
#' its specificity and sensitivity in terms of covering of the red tumors or the fact that it can be extended or not to a better solution.
#' 
#' Function \link{do.musette} does all the above and returns both the 'alterations'
#' data frame (with informations on indiviual alterations and their
#' domination), and the 'solutions' data frame.  Function \link{musette2csv} returns a data frame containg the same info
#' as d but items which were vectors are converted into semicolon-delimited
#' strings.
#' 
#' @name musette-package
#' @aliases musette-package musette
#' @docType package
#' @import Rcpp 
#' @importFrom Rcpp evalCpp
#' @useDynLib musette
#' @author Thomas Picchetti, Jennifer Wong, Etienne Birmelé.
#' 
#' Maintainer: Thomas Picchetti ans Etienne Birmelé <etienne.birmele@@parisdescartes.fr>
#' @references This optional section can contain literature or other references
#' for background information.
#' @keywords package
#' @examples cf vignette
#' 
#' 
NULL

