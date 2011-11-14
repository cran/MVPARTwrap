
\name{summary.MRT}
\alias{summary.MRT}
\title{Returns a detailed listing of a fitted MRT object.}
\description{Gives a summary of the results contained in an MRT object : output shows an overview of the information available for each node.}
\usage{\method{summary}{MRT}(object,IndvalPART=TRUE,IndvalNODE=TRUE,\dots) 
} 
\arguments{

  \item{object }{fitted model object of class \code{MRT}}
  \item{IndvalPART}{Logical. Print IndVal species for the final partition}
  \item{IndvalNODE}{Logical. Print IndVal species for each node}
  \item{... }{arguments to be passed to or from other methods.}
}

\details{This function is a method for the generic function summary for class \code{MRT}. It can be invoked by calling summary for an object of the appropriate class. Dicriminant species (who contribute the most to the deviance), their mean at each subgroup formed by the split (left and right), Indval species at the node level and at the full partition level, and the complexity (contribution to \eqn{R^{2}}{R^2}) is available.}

\seealso{summary, MRT} 

\examples{ 

  ## MRT of the spider data available in the mvpart package. ##
# ------------------------------------------------------------- #

require(mvpart)
data(spider)

result<-mvpart(data.matrix(spider[,1:12]) ~ herbs + reft + moss + 
sand + twigs + water, spider, xv='1se', xval=10, xvmult=100)

## Format the object of class mvpart into MRT ##
res<-MRT(result,10)

## Summary of this MRT object ##
summary(res)

## This summary has the following structure ##
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#       --- Node # --- : # 
# corresponds to the number given on the tree plot
# Complexity (r2) ##
# left child  right child
# ~ Discriminant species : 
#                      		 sp1         sp2        ...
# % of expl. deviance of##   46.2782706 11.5209795  ...
# Mean value on the left     2.0758065  0.0000000   ...
# Mean value on the right    0.5257732  0.7371134   ...

# ~ INDVAL species for this node: : left is 1, right is 2
#    cluster indicator_value probability
# SPA       1          0.8094       0.001
# SPB       1          0.7193       0.001
# SPC       1          0.6496       0.005
# SPD       2          1.0000       0.001
# ...       ...        ...          ...
# Note : the output will only show the species which have a higher 
#'percent' value than the one specified by the user when the MRT object
# was created. See ?MRT for more details
} 

\keyword{multivariate}
\keyword{tree}