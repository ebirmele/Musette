% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rawDominationGraphFineTuned}
\alias{rawDominationGraphFineTuned}
\title{Compute the Domination Relationship in the Current Graph}
\usage{
rawDominationGraphFineTuned(
  chromosome,
  locus,
  distance,
  redPercent,
  bluePercent
)
}
\arguments{
\item{chromosome}{(for genomic alterations) An integer vector containing the chromosome of each
alteration in the current graph (see \code{\link{takeGraph}}). If omitted, this parameter has no effect.}

\item{locus}{(for genomic alterations) An integer vector with the same length as \code{chromosome}, containing
the 'locus' of each alteration (e.g. in base pairs).  If omitted, this parameter has no effect.}

\item{distance}{The maximum distance accepted between two alterations' loci for
domination to take place. If omitted, this parameter has no effect.}

\item{redPercent}{The required percentage of red samples}

\item{bluePercent}{The required percentage of blue samples}
}
\value{
A list \code{l} of 6 named lists of character vectors. If \code{A} is the name of an alteration, then :
\itemize{
\item \code{l$loose_dominators[[A]]} contains the names of alterations "loosely" dominating \code{A}.
\item \code{l$loose_dominated[[A]]} contains the names of alterations "loosely" dominated by \code{A}.
\item \code{l$strict_dominators[[A]]} contains the names of alterations "strictly" dominating \code{A}.
\item \code{l$strict_dominated[[A]]} contains the names of alterations "strictly" dominated by \code{A}.
}
}
\description{
(See details below about what "domination" means) 
Compute for every alteration A :
  \itemize{
    \item the list of alterations dominated by A "loosely" (i.e. including those dominating A)
    \item the list of alterations dominated by A "strictly" (i.e. excluding those dominating A)
    \item the list of alterations dominating A "loosely" (i.e. including those dominated by A)
    \item the list of alterations dominating A "strictly" (i.e. excluding those dominated by A)
  }
}
\details{
This is an inner function, called by other R functions in this
package, but it is left available for the expert user or developer.

Given a network connecting samples and alterations, and given two groups of samples, red and blue, alteration A
will dominate alteration B if : 
\itemize{
 \item At least <redPercent>\% of the red samples affected by B are also affected by A, 
 \item At least <bluePercent>\% of the blue samples affected by A are also affected by B,
 \item The loci of A and B lie on the same chromosome, at most <distance> base pairs apart (this part is for alterations which are genomic, and can simply be ignored for other applications)
}
This is a somewhat more flexible version of the more natual criteria according to which A would dominate B if : 
\itemize{
 \item The red samples affected by B are also affected by A
 \item The blue samples affected by A are also affected by B
} 
To compute this simpler and more restrictive relation, use \code{redPercent=100} and \code{bluePercent=100}.
}
\seealso{
\code{\link{takeGraph}}, \code{\link{domination.graph}}

Other inner functions: 
\code{\link{connectedComponents}()},
\code{\link{encodeSolutionTree}()},
\code{\link{exploreTree}()},
\code{\link{fullyGrowTree}()},
\code{\link{growTree}()},
\code{\link{resetSolutionTree}()},
\code{\link{scc}()},
\code{\link{setBound}()},
\code{\link{setThreshold}()},
\code{\link{takeGraph}()},
\code{\link{takeGroups}()},
\code{\link{takeStepScore}()}
}
\concept{inner functions}
