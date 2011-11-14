# Arguments
# Y : response (transformed if needed)
# VSG : eVSGplanatory matriVSG
# VA : covariates
# PART : initial partition
# ... : further arguments to be passed to mvpart.

CascadeMRT<-function(Y,VSG,VA,xv1='min',xvSUB='min',minbucket1=5,minbucketSUB=2,cp1=0.1,cpSUB=0.01,...)
{
	# Format of matrices : all have to be data.frame
	if(!is.data.frame(Y)) stop('Y should be a data.frame')
	if(!is.data.frame(VSG)) stop('VSG should be a data.frame')
	if(!is.data.frame(VA)) stop('VA should be a data.frame')
	
	# If not previously loaded
	#if(!any(search()=="package:mvpart")) require(mvpart)
	
	# Initialisation
	drops<-list()
	
	# First drop
	drops[[1]]<-mvpart(data.matrix(Y)~data.matrix(VA),data=cbind(Y,VA),xv=xv1,minbucket=minbucket1,cp=cp1,...)
	
	# Subsequent drops : one for each previous leafs
	part1<-drops[[1]]$where
	where1<-unique(drops[[1]]$where)
	
	# Get residuals from first drop
	Yrw<-Y-predict(drops[[1]])
	
	for(i in 1: length(where1))
	{
		# Add interaction for subsequent drops here
		
		# Selection of the objects in this leaf-drop to be
		Yrwi<-Yrw[which(part1==where1[i]),]
		VSGi<-VSG[which(part1==where1[i]),]
		drops[[i+1]]<-mvpart(data.matrix(Yrwi)~data.matrix(VSGi),data=cbind(Yrwi,VSGi),xv=xvSUB,minbucket=minbucketSUB,cp=cpSUB,...)
	}
	
	# Add information for R2 here : folloVA diagram
	res<-list(drops,part1,where1,Y)
	names(res)<-c('drops','part1','where1','Y')
	class(res)<-'CascadeMRT'
	return(res)
}

# This function builds the R2 
CasMRTR2<-function(obj,NodeADMIT=obj$where1)
{

	if(class(obj)!='CascadeMRT') stop : 'The object provided was not of class CascadeMRT'
	# Initilisation
	PropRE_noeud_drops<-list()
	dropsR2<-list()
	SST<-sum(apply(obj$Y,2,var))*(nrow(obj$Y)-1)
	
	# The R2 of the first drop
	drops1<-MRT(obj$drops[[1]],10)
	drops1R2<-sum(drops1$R2)
	
	
	
	# Get residuals from first drop
	Yrw<-obj$Y-predict(obj$drops[[1]])
	SSE_TOT<-SST*(1-drops1R2)
	
	# The RE attributed to each node from the first drop
	for(i in 1:length(obj$where1))
	{
		Yrwi<-Yrw[which(obj$part1==obj$where1[i]),]
		MEANSi<-matrix(rep(apply(Yrwi,2,mean),each=nrow(Yrwi)),nrow(Yrwi),ncol(Yrwi),byrow=FALSE)
	
		PropRE_noeud_drops[[i]]<-sum((Yrwi-MEANSi)^2)/SSE_TOT
		# 
		# The R2 of subsequent drops
		dropsR2[[i]]<-sum(MRT(obj$drops[[i+1]],10)$R2)
	}
	# - # - # Building diagrame VAith polygons # - # - #

	# Initial diagram is 1 VSG 1.
	par(mar=c(1,1,1,1))
	plot(c(0,1),c(0,1),type='n',xlab='',ylab='',xaxt='n',yaxt='n',bty='n',asp=1)
	# Square borders
	lines(c(0,0),c(0,1))
	lines(c(0,1),c(0,0))
	lines(c(1,0),c(1,1))
	lines(c(1,1),c(1,0))

	# R2 portion explained by VA
	
	polygon(c(0,drops1R2,drops1R2,0), y = c(0,0,1,1),col='grey')
	# Lable 'VA' and portion of eVSGplained variance
	text(drops1R2/2,0.5,labels='VA')
	text(drops1R2/2,0.5,labels=paste(round(drops1R2*100,2),'%'),pos=1,cex=0.75,col='grey30')
	

	# Prop RE by node of the first drop x total RE of first drop + R2 of first drop => the horizontal width of the node bar
	HORIZ_BREAK<-c(0,cumsum(unlist(PropRE_noeud_drops)*(1-drops1R2)))+drops1R2
	
	
	# R2 of each subsequent drop => gives the limit vertical position in the node boxes.
	VERT_BREAK<-unlist(dropsR2)
	
	
	# Build the per node R2 portion boxes and calculate total R2
	R2_TOTAL=drops1R2
	
	for(j in 1:length(obj$where1))
	{
		if(any(NodeADMIT==obj$where1[j]))
		{
			polygon(c(HORIZ_BREAK[j],HORIZ_BREAK[j+1],HORIZ_BREAK[j+1],HORIZ_BREAK[j]), y = c(0,0,VERT_BREAK[j],VERT_BREAK[j]),col='grey')
		
			# Node label vertically in the middle, at the right of the R2 line + R2 in %
			middle_horiz<-(HORIZ_BREAK[j]+HORIZ_BREAK[j+1])/2
			text(middle_horiz,VERT_BREAK[j],labels=paste('#',obj$where1[j],sep=''),pos=3,crt=90,offset=0.1,cex=0.75)
			text(middle_horiz,VERT_BREAK[j],labels=paste(round((HORIZ_BREAK[j+1]-HORIZ_BREAK[j])*VERT_BREAK[j]*100,2),'%',sep=''),pos=1,col='grey30',crt=90,offset=0.2,cex=0.60)
			
			# Contribution of each subsequent drop : width*height
			R2_TOTAL=R2_TOTAL+(HORIZ_BREAK[j+1]-HORIZ_BREAK[j])*VERT_BREAK[j]
		}
	}
	
	# Add unexplained variation in corner : 
	text(0.9,0.9,labels=paste('Unexplained variation : ',round(100-R2_TOTAL*100,2),'%'),pos=2,col='grey30',offset=0.5,cex=0.75)
	
	# Plot GLobal R2 and subsequent drops contribution outside the square plot
	#text(0.9,0.9,labels=paste('Global : ',round(R2_TOTAL*100,2),'%'),pos=2,col='grey30',offset=0.5,cex=0.75)
	#text(0.9,0.85,labels=paste('Total SG : ',round((R2_TOTAL-drops1R2)*100,2),'%'),pos=2,col='grey30',offset=0.5,cex=0.75)
	

}