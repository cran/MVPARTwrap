\name{CasMRTR2}
\alias{CasMRTR2}
\title{Diagram of explained variation for a cascade MRT analysis} 
\description{Diagram of explained variation for a cascade MRT analysis. See details for more precisions.}
 
\usage{CasMRTR2(obj,NodeADMIT=obj$where1)} 

\arguments{
  \item{obj}{Object of class \code{CascadeMRT}}
  \item{NodeADMIT}{Vector with numbers corresponding to the group membership of the first drop to be admitted in the diagram. Default is all are included.}
}

\value{Graphical output of explained variation.}


\details{This diagram represents the percentage of variation explained by each drop. Explained variance of the original response is in greyed area, unexplained in is white. The rectangle representing the explained variation of the first drop MRT model is annoted VA and is placed to the left. The area is proportional to the explained percentage which is written under the VA annotation. The other greyed areas each represent a model from the second wave. Each area has a height equal to the explained variance of the model build on the node objects and a width equal to the variance unexplained in that node by the first drop. This provides a proportional area to the variation explained of the original response.}

\seealso{\link{mvpart},\link{CascadeMRT}}

\references{

Dea'th, G. 2002. Multivariate regression tree : a new technique for modeling species-environment relationships. Ecology 83:1105-1117.

Ouellette, M.-H. & P. Legendre. An adjusted R^2 statistic for multivariate regression tree analysis (Manuscript).

Ouellette, M.-H., Legendre, P. & B. Borcard. Cascade Multivariate Regression Tree: a novel approach for modelling nested explanatory sets (accepted for publication in Methods in Ecology and Evolution).

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
# (accepted).

data(doubs)

# Data importation #
spe<-doubs$poi
envmod<-doubs$mil
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

