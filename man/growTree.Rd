% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{growTree}
\alias{growTree}
\title{Perform One Growth Step of the Solution Tree}
\usage{
growTree()
}
\value{
nothing.
}
\description{
This is an inner function, called by other R functions in this
package, but it is left available for the expert user or developer.
}
\details{
Increase the step-core threshold to the next possible value (which is the
step-score of the first node in the waiting queue). This allows one or 
many nodes to be added to the solution tree.
}
\seealso{
Other inner functions: 
\code{\link{connectedComponents}()},
\code{\link{encodeSolutionTree}()},
\code{\link{exploreTree}()},
\code{\link{fullyGrowTree}()},
\code{\link{rawDominationGraphFineTuned}()},
\code{\link{resetSolutionTree}()},
\code{\link{scc}()},
\code{\link{setBound}()},
\code{\link{setThreshold}()},
\code{\link{takeGraph}()},
\code{\link{takeGroups}()},
\code{\link{takeStepScore}()}
}
\concept{inner functions}
