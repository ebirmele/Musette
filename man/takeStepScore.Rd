% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{takeStepScore}
\alias{takeStepScore}
\title{Pass a Step-score Mode to the C++ Code}
\usage{
takeStepScore(stepMode)
}
\arguments{
\item{stepMode}{: an integer. 0 for "original", 1 for "best-first", 2 for "strict"}
}
\value{
nothing.
}
\description{
This is an inner function, called by other R functions in this
package, but it is left available for the expert user or developer.
}
\seealso{
Other inner functions: 
\code{\link{connectedComponents}()},
\code{\link{encodeSolutionTree}()},
\code{\link{exploreTree}()},
\code{\link{fullyGrowTree}()},
\code{\link{growTree}()},
\code{\link{rawDominationGraphFineTuned}()},
\code{\link{resetSolutionTree}()},
\code{\link{scc}()},
\code{\link{setBound}()},
\code{\link{setThreshold}()},
\code{\link{takeGraph}()},
\code{\link{takeGroups}()}
}
\concept{inner functions}
