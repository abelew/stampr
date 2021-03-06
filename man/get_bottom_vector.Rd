% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtisan.r
\name{get_bottom_vector}
\alias{get_bottom_vector}
\title{Given a vector and set of steps, use rmvhyper to make a new vector.}
\usage{
get_bottom_vector(vector, steps)
}
\arguments{
\item{vector}{input vector}

\item{steps}{Set of steps over which to iterate.}
}
\value{
New vector containing the number of times the rhvhyper
 distribution was not zero for each element of the original vector.
}
\description{
This takes an input vector and for every element of it provides a
new element which is calculated by doing the following:
1.  Run rmvhyper on the original vector in order to generate a set
of random numbers using the multivariate hypergeometric
distribution.
2.  Finding the number of values from #1 which are not 0.
3.  Putting the number from #2 into the new element.
4.  Finally, returns the new vector.
}
\details{
It seems to me that there must be a much more efficient way to get
this information?
}
