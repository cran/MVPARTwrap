\name{plot_tree} 
\alias{plot_tree} 
\title{Plot Multivariate Regression Tree} 
\description{ Plot mvpart tree according to the relative error } 
 
\usage{ 
plot_tree(obj, NodeMarking, R2, Cex, widthtree, heighttree, 
LABELS, R2A=FALSE, X=NULL, T=NULL, \dots)
} 
\arguments{ 
\item{obj}{mvpart object} 
\item{NodeMarking}{TRUE or FALSE, should the level of \eqn{{R}^2}{R^2} be printed}
\item{R2}{Matrix of contribution to \eqn{{R}^2}{R^2} (object x splits)}
\item{Cex}{see par(cex)} 
\item{widthtree}{width of plot window}
\item{heighttree}{height of plot window}
\item{LABELS}{Given labels for plotting. Used in partial analysis.}
\item{R2A}{logical, should the adjusted \eqn{{R}^2}{R^2} (see references) be plotted under the tree}
\item{X}{explanatory matrix, must be provided if R2A=TRUE}
\item{T}{Number of runs for GDF estimation, must be provided if R2A=TRUE}
\item{\dots}{Further arguments to be passed to the plotting devices or functions. See \code{Details} for more details.}
} 

\value{plotted tree (no value returned)}

\keyword{Multivariate regression tree}