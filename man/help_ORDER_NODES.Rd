\name{ORDER_NODES} 
\alias{ORDER_NODES} 
\title{Order mvpart nodes in terms of R^2} 
\description{ Order mvpart nodes in terms of R^2 for better tree output} 
 
\usage{ 
ORDER_NODES(obj,R2,splits)
} 
\arguments{ 
\item{obj}{mvpart object} 
\item{R2}{matrix of explained variation (objects x splits)}
\item{splits}{contains the line numbers of table 'frame' of the mvpart object that are nodes, not leafs}
} 

\value{Vector of modes ordering in terms of R^2}

\seealso{\code{\link{mvpart}} \code{\link{MRT}}}

\keyword{Multivariate regression tree}