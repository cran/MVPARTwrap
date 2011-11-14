# This function estimates the generalized degrees of freedom according to : YE, Jianming. 1998. On Measuring and Correcting the Effects of Data Mining and Model Selection. Journal of the American Statistical Association, Vol. 93, No. 441, pp. 120-131.

# Definition of the arguments :
# TREE : object of class 'MRT'
# T : number of times the algorithm should run
# Tau_const : tuner parameter value
# error_type : type of errior to add as perturbation
# Explanatory matrix (for MRT objects)

# Output : Adjusted R2

R2aGDF<-function(TREE,T,tau_const,X=NULL,...)
{
	# Initialisation
	require(mvpart)
	
	# Defining the response (Y) and explanatory matrices (X), number of objects n, number of species s, variance of each species sigmas, coefficient of determination of the tree R2
	
	
	
	if(class(TREE)=='MRT')
	{
		if(is.null(TREE$obj$y)) stop ('Use argument y=TRUE when building the mvpart object so the response matrix is readily available')
		if(is.null(X)) stop ('Must provide argument X as non null for the explanatory matrix to be readily available')
		Y<-TREE$obj$y
		X<-X
		n<-nrow(Y)
		# Number of columns (species) in Y
		s<-ncol(Y)
		# Variance of columns of Y : sigmas
		sigmas<-diag(var(Y))
		R2<-sum(TREE$R2)
		size<-length(unique(TREE$obj$where))
	}
	
	# Tuning parameter : tau
	tau<-tau_const*sqrt(sigmas)
	# Simulation results
	DELTA_T<-mat.or.vec(n,s)
	FIT<-mat.or.vec(n,s)
	# Sensitivities by species
	diagB<-mat.or.vec(n,s)
	
	# Estimate of D(M) using Monte Carlo : multivariate algorithm 1 adapted from Ye 1998.
	for(i in 1:T)
	{
		# Generate the delta_t for each column of Y, a matrix of object by species, each column a random standart normal variable ponderated by the variance of each species.
		delta_t<-mat.or.vec(n,s)
		for(j in 1:s)
		{
			delta_t[,j]<-rnorm(n)*tau[j]
		}

		
		# Evaluate the predicted values of Y+delta_t using the modelling technique of choice  : mvpart.
		invisible(capture.output(fit<-predict(mvpart(data.matrix(Y+delta_t)~data.matrix(X),data=as.data.frame(cbind(Y+delta_t,X)),size=size,minbucket=1,plot.add = FALSE),type='matrix')))
		
		#,...),type='matrix')))
		
		
		# Stock fit and delta_t for future use : stocking by rows, thus each permutation result is stacked on top of one another.
		DELTA_T<-rbind(DELTA_T,delta_t)
		FIT<-rbind(FIT,fit)
		
	}	
	
	# Get rid of empty values previously added
	DELTA_T<-DELTA_T[-c(1:n),]
	FIT<-FIT[-c(1:n),]
	
	# *** Use linear modelling to assess the sensitivity of each species for each object, thus each regression slope. We we're going to use the rda formulation, until we realized thatt(DELTA_Ti)%*%DELTA_Ti) is not always inversible. (--We use the rda formulation to calculate B, all slopes of linear regression, and take the diagonal to get the slope , for each species seperatly, of the corresponding columns of FIT and DELTA_T.--)
	
	for(i in 1:n)
	{
		# Seperating the results for each object as required by the algorithm. lines_obji states the lines of both DELTA_T adn FIT that correspond to the same object.
		lines_obji<-seq(from = i, to = n*T, by=n)
		DELTA_Ti<-DELTA_T[lines_obji,]
		FITi<-FIT[lines_obji,]
		# Need to transpose the matrix so we have the species as colomuns and the perturbations as lines
		#DELTA_Ti<-t(DELTA_Ti)
		#FITi<-t(FITi)
		for(j in 1:s)
		{
			diagB[i,j]=coefficients(lm(FITi~DELTA_Ti))[2]
			#diag(inv(t(DELTA_Ti)%*%DELTA_Ti)%*%t(DELTA_Ti)%*%FITi)
		}
	}
		
	# Sum all lines of this matrix to get sensitivity of each species
	sensivit<-colSums(diagB)
	
	
	# Calculate MS(full)
	DF<-sensivit[1] # Degrees of freedom for the first species : the same for all.
	R2a<-1-(1-R2)*((n-1)/(n-DF))
	rm(n,s,sigmas,tau,DELTA_T,FIT,diagB)
	gc(verbose=FALSE)
	
	
	
	# Calculate MS(reduced)
	
	return(R2a)	
}
