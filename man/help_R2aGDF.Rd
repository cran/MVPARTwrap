\name{R2aGDF} 
\alias{R2aGDF} 
\title{Adjustment for \eqn{{R}^2}{R^2} of Multivariate regression tree analysis (MRT) based on general degrees of freedom (GDF) estimation} 
\description{ This function returns the adjusted \eqn{{R}^2}{R^2} for an MRT model} 
 
\usage{ 
R2aGDF(TREE,T,tau_const,X=NULL,...)
} 
\arguments{ 
\item{TREE}{object of class \code{MRT}.} 
\item{T}{The number of times the algorithm should run. }
\item{tau_const}{constant for tuner parameter multiplicative value.}
\item{X}{explanatory matrix (must be provided when TREE is an \code{MRT} object)}
\item{...}{further options to provided to mvpart}
} 


\value{
\item{R2a}{Adjusted \eqn{{R}^2}{R^2} for the provided TREE}
}
\details{ This function performs the adjustment of the coefficient of determination for an MRT model. The adjustment is based on Ye's algorithm 1 used to estimate the number of general degrees of freedom. The adjustement itself is defined in Ouellette & Legendre (Manuscript)}

\seealso{\code{\link{mvpart}}, \code{\link{MRT}}}
\references{

	Ye, J. (1998). "On Measuring and Correcting the Effects of Data Mining and Model Selection." Journal of the American Statistical Association 93(441): 120-131.
	
	Dea'th, G. 2002. Multivariate regression tree : a new technique for modeling species-environment relationships. Ecology 83:1105-1117.

	Ouellette, M.-H. & P. Legendre. An adjusted \eqn{{R}^2}{R^2} statistic for multivariate regression tree analysis. (Manuscript)
	}
 
\examples{ 
                     # -------------------- #
                     # -   MRT analysis   - #
                     # -------------------- #
# ------------------------------------------------------------- #
require(vegan)
require(labdsv)
data(mite)
data(mite.env)
mite.hell<-decostand(mite,method='h')


# Polynomial water content
WaterCont.poly = poly(mite.env[,2], 3)
# All polynomial of water content are retained 
# (according to forward selection)

mite.tree<-mvpart(data.matrix(mite.hell) ~ data.matrix(mite.env[,-2])
+ WaterCont.poly, data=as.data.frame(cbind(mite.env,WaterCont.poly)),
xv='min',xvmult=1000)

# Get graphic window ready #
quartz(width=10, height=7.5)
par(cex=0.4)
# Get prettier tree #
res<-MRT(mite.tree,10)
# Make prettier lables
LABELS<-res$mat_labels
LABELS[1,2]<-'Water.p1>=-0.02159'
LABELS[1,3]<-'Water.p1< -0.02159'

LABELS[2,2]<-'SubsDens>=51.8'
LABELS[2,3]<-'SubsDens<51.8'

LABELS[3,2]<-'Water.p1>=-0.07365'
LABELS[3,3]<-'Water.p1< -0.07365'

LABELS[4,2]<-'SubsDens < 36.38'
LABELS[4,3]<-'SubsDens >=36.38'

LABELS[5,2]<-'Topo< 1.5'
LABELS[5,3]<-'Topo>= 1.5'

LABELS[6,2]<-'Shrub< 1.5'
LABELS[6,3]<-'Shrub>= 1.5'

MRT.mite.tree<-MRT(mite.tree,10,LABELS=LABELS)
summary(MRT.mite.tree)
plot(MRT.mite.tree,type='tree')

# Show summary of results
summary(MRT.mite.tree)

# Calculate R2a
R2aGDF(MRT.mite.tree, T=100, tau_const=0.6, X=cbind(mite.env[,-2],
WaterCont.poly))

} 

\keyword{Multivariate regression tree}