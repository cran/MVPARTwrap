\name{code.tree} 
\alias{code.tree} 
\title{Multivariate regression tree analysis dummy coding} 
\description{This function dummy codes an mvpart tree} 
 
\usage{ 
code.tree(obj,LWHERE,RWHERE,LABELS)
} 

\arguments{ 
\item{obj}{obj of class mvpart} 
\item{LWHERE}{list of left objects in every node}
\item{RWHERE}{list of right objects in every node}
\item{LABELS}{\code{mat_lables} object} 
} 

\value{
\item{code}{Matrix representing the coded tree}}

\details{This function is used internally by the \code{plot.MRT} generic function and is not intended for external use. Either way, each coded variable in matrix 'code' represents one of the splits of the MRT. The objects on the left and right of the split have weights corresponding to the inverse of the number of observations in each group. The left objects get positive weights and the right objects get negative weights. The objects not concerned with a split have a weight of zero. As a consequence, each column sums to 0. Coded by Marie-Helene Ouellette, August 2009, modified from Pierre Legendre, August 2009}

\keyword{Multivariate regression tree}
\keyword{dummy coding}