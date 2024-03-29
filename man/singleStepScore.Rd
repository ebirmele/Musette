% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{singleStepScore}
\alias{singleStepScore}
\title{Original Step-score}
\usage{
singleStepScore(r0, b0, r1, b1)
}
\arguments{
\item{r0}{The number of red samples to start from}

\item{b0}{The number of blue samples to start from}

\item{r1}{The number of red samples to get to}

\item{b1}{The number of blue samples to get to}
}
\value{
the step-score value. If the parameters have meaningless values 
(negative or greater than the number of red/blue samples) the result is -1.
}
\description{
Computes the original step-score from a solution to another.
}
\details{
This function can only be called after \link{solutions} has been called
or at least \link{takeGraph}, to create the data structure which 
contain the score distributions. \cr
Only the number of red and blue samples altered in the start and finish
solutions are needed. The corresponding solutions do not have to 
actually exist.
}
