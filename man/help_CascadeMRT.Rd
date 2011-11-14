\name{CascadeMRT}
\alias{CascadeMRT}
\title{Cascade of Multivariate Regression Trees Analysis} 
\description{Function CascadeMRT performs a cascade multivariate regression tree analysis using \code{mvpart} as a basic function to persue the multivariate tree analyses. The first model is at the top of the cascade and is called the first drop. All models built on the leafs of the first drop are called subsequent drops and consist the second wave. See Ouellette, Legendre & Borcard (accepted) for further informations and details.}

\usage{

CascadeMRT(Y, VSG, VA, xv1='min', xvSUB='min', minbucket1=5,
minbucketSUB=2, cp1=0.1, cpSUB=0.01, ...)}

\arguments{
  \item{Y}{Response matrix}
  \item{VSG}{Explanatory variables used for second wave} 
  \item{VA}{Explanatory variables used for first wave}
  \item{xv1}{\code{xv} argument for the first wave}
  \item{xvSUB}{\code{xv} argument for second wave}
  \item{minbucket1}{\code{minbucket} argument for the first wave}
  \item{minbucketSUB}{\code{minbucket} argument for second wave}
  \item{cp1}{\code{cp} argument for the first wave}
  \item{cpSUB}{\code{cp} argument for second wave}
  \item{... }{further arguments given to \code{mvpart}. See \code{rpart.control} and \code{mvpart} for further details. These controls will be used for both the first wave and second wave of models.}
}

\value{An object of class \code{CascadeMRT}, basically a list including the following elements

  \item{drops}{List of mvpart models built, the first drop follow by the subsequent drops given in the order of \code{where1}.}
  \item{part1}{Unique leaf numbers of the first drop}
  \item{where1}{Partition membership of the first drop}
  \item{Y}{Original response matrix. Useful for use of \code{CasMRTR2} function.}}
  
\details{\code{cp1} and \code{cpSUB} arguments are used to specify the complexity parameter of each wave (\code{cp1} for the principal, and \code{cpSUB} for the subsequent wave). The complexity parameter of an MRT model is the minimum contribution to the \eqn{{R}^2}{R^2} of the tree for a split to be considered. As the main explanatory set is set to explain the largest response composition variation, the complexity parameter \code{cp1} could be set to be higher than the usual 0.01 value (by default 0.1).}
 

\seealso{\code{\link{mvpart}},\code{\link{CasMRTR2}}}

\references{
Dea'th, G. 2002. Multivariate regression tree : a new technique for modeling species-environment relationships. Ecology 83:1105-1117.

Ouellette, M.-H. & P. Legendre. An adjusted \eqn{{R}^2}{R^2} statistic for multivariate regression tree analysis. Manuscript.

Ouellette, M.-H., P. Legendre and D. Borcard. 2012. Cascade Multivariate Regression Tree: a novel approach for modelling nested explanatory sets. Methods in Ecology & Evolution 3: 234-244. 
} 

\examples{ 
# Initialisation

require(vegan)
require(mvpart)
require(MVPARTwrap)
require(ade4)

                   ### -- Doubs example --- ###
# ------------------------------------------------------------- #
# Note that we modify the data found in ade4 to match the data  
# set provided with Numerical Ecology in R (Borcard et al 2011),  
# the data set used in the paper Ouellette, Legendre & Borcard 
# (2012).

data(doubs)

# Data importation #
spe<-doubs$fish
envmod<-doubs$env
envmod[,1]<-envmod[,1]/10
envmod[7,1]<-36.8
envmod[,3]<-round((exp(envmod[,3])-1)/10,1)
envmod[,4]<-envmod[,4]/100
envmod[,5]<-envmod[,5]/10
envmod[,7]<-envmod[,7]/100
envmod[,8]<-envmod[,8]/100
envmod[,9]<-envmod[,9]/100
envmod[,10]<-envmod[,10]/10
envmod[,11]<-envmod[,11]/10

envmod_nodas<-envmod[,-1]
spe.hel <- decostand(spe, "hellinger")

POISSON<-CascadeMRT(Y=spe.hel, VSG=envmod_nodas, 
VA=as.data.frame(envmod[,1]), xvmult=100)
CasMRTR2(POISSON,NodeADMIT=3)

               ### -- Oribatid mite example --- ###
# ------------------------------------------------------------- #
data(mite)
data(mite.env)

mite.hel <- decostand(mite, "hellinger")
W<-mite.env[,4]

MITE<-CascadeMRT(Y=mite.hel, VSG=mite.env[,-4], VA=data.frame(W), 
xvmult=500,xvSUB='pick',cp1=0.1)
CasMRTR2(MITE,NodeADMIT=c(3,2))
} 
\keyword{Cascade}
\keyword{multivariate}
\keyword{tree}

