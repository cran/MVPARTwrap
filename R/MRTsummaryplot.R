# ---------------------------------------------------------- #
# Get determinant species : this has to be done at each node #
# ---------------------------------------------------------- #

MRT<-function(obj,percent=10,species=NULL,LABELS=FALSE,...)
# This function recalculates the complexity of each node, giving the species' contribution to the R2 at each node, thus the output is the table 1 in Dea'th (2002). It contains the total species variance partitionned by species, by the  tree, and by the splits of the tree.

{
	
	#require(Hmisc)

# obj is the mvpart object
# percent : percentage considered discriminant species
# typeplot : a vector with states tree the usual tree output is needed, and triordi if the triplot ordination representation of the tree is of interest
# ... : further options for mvpart. See rpart.option, rpart and mvpart functions for more details
# species : A vector of species names used to build the table of partitioned variance. If it is NULL, the colnames of obj$y will be used.

	# Check if obj is of proper class
	if(class(obj)!='rpart') stop('obj is not of class rpart')
	
	# splits contains the line numbers of table 'frame' of the mvpart object that are nodes, not leafs
	splits<-which(obj$frame[,1]!="<leaf>")
	
	# node numbers NOT in order of their contribution (their is no need yet)
	nodes<-as.numeric(row.names(obj$frame)[splits])  

	# Initilizing the matrix to calculate the contributions to R2 of each species, at each node
	SSE<-mat.or.vec(length(nodes)*2+1,ncol(obj$y))
	MOYs<-mat.or.vec(length(nodes)*2+1,ncol(obj$y))
	moy<-apply(obj$y,MARGIN=2,FUN=mean)	
	MOY<-rep(moy,times=nrow(obj$y))
	MOY<-matrix(MOY,nrow(obj$y),length(moy),byrow=TRUE)
	SSE[1,]<-colSums((obj$y-MOY)^2)
	# List of the objects in each node
	LWHERE<-list()
	RWHERE<-list()
	
	compt_SSE<-2 # Already got the SST... so we start at the second slot
	
	# Loop that computes the SSEs for all nodes (2 per nodes, 1 per child), and create LWHERE and RWHERE
	
	for(i in 1:length(nodes)) 
	{

		# Initialisation of :
		
		# left child a leaf, than this will change to nodes[i]*2
		Lleaf=0
		# rigth child a leaf, than this will change to nodes[i]*2+1
		Rleaf=0
		# left child a node, than this will change to nodes[i]*2
		Lnode=0
		# rigth child a node , than this will change to nodes[i]*2+1
		Rnode=0
		
		# What are the two children ? Are they leafs or nodes ?
		# Left child ?
		ifelse(any(nodes==nodes[i]*2),Lnode<-nodes[i]*2,Lleaf<-nodes[i]*2) # is it a node or a leaf ?
		# Rigth child ?
		ifelse(any(nodes==nodes[i]*2+1),Rnode<-nodes[i]*2+1,Rleaf<-nodes[i]*2+1) # is it a node or a leaf ?
		
		
		# For each child seperatly (the bipartition) , 
		
		# -- # Left child # -- #
		
		if(Lleaf==0)
		{
			
			# Initialisation
			node=0
			leaf=0
			ifelse(any(nodes==Lnode*2),node<-Lnode*2,leaf<-Lnode*2)
			ifelse(any(nodes==Lnode*2+1),node<-c(node,Lnode*2+1),leaf<-c(leaf,Lnode*2+1))
			
			# We know at this point that 'node' is not empty and Lnode was a node (at least the right one)

			for(j in 1:length(nodes)) # Since left child is a node : check younger relatives for leafs
			{
				# Get this node's children if it is one decendant of the left node
				
				if(any(node==nodes[j])) # if the node we're workin on is in the node vector (left child of interest) : look at it's children and classify them
				{
					# Check if these children are nodes or leafs and update consequently 'node' and 'leaf' ('leaf' is really the important one, but we need 'node' for 'leaf')
					ifelse(any(nodes==nodes[j]*2),node<-c(node,nodes[j]*2),leaf<-c(leaf,nodes[j]*2))
					ifelse(any(nodes==nodes[j]*2+1),node<-c(node,nodes[j]*2+1),leaf<-c(leaf,nodes[j]*2+1))
				}
			}
			
			
	
		# Than create the Lwhere : which object are in the left node
		
		Lwhere<-0
		for(k in 1:length(leaf)) Lwhere<-c(Lwhere,which(obj$where==which(as.numeric(row.names(obj$frame))==leaf[k])))
		Lwhere<-Lwhere[-1]
		
		}
		
		if(Lnode==0)
		{
			
			Lwhere<-which(obj$where==which(as.numeric(row.names(obj$frame))==nodes[i]*2))
			
		}
		
		
		
		# -- # Rigth child # -- #
	
		if(Rleaf==0)
		{
			
			# Initialisation
			node=0
			leaf=0
			ifelse(any(nodes==Rnode*2),node<-Rnode*2,leaf<-Rnode*2)
			ifelse(any(nodes==Rnode*2+1),node<-c(node,Rnode*2+1),leaf<-c(leaf,Rnode*2+1))
			
			# We know at this point that 'node' is not empty and Rnode was a node (at least the right one)
			# Since right child is a node : check younger relatives for leafs
			for(j in 1:length(nodes)) 
			{
				# Get this node's children if it is one decendant of the rigth node
				
				if(any(node==nodes[j])) # if the node we're workin on is in the node vector (rigth child of interest) : look at it's children and classify them
				{
					# Check if these children are nodes or leafs and update consequently 'node' and 'leaf' ('leaf' is really the important one, but we need 'node' for 'leaf')
					ifelse(any(nodes==nodes[j]*2),node<-c(node,nodes[j]*2),leaf<-c(leaf,nodes[j]*2))
					ifelse(any(nodes==nodes[j]*2+1),node<-c(node,nodes[j]*2+1),leaf<-c(leaf,nodes[j]*2+1))
					
				}
			}
			
		
			# Than create the Rwhere : the objects in the rigth child
			Rwhere<-0
			for(k in 1:length(leaf)) Rwhere<-c(Rwhere,which(obj$where==which(as.numeric(row.names(obj$frame))==leaf[k])))
			Rwhere<-Rwhere[-1]
			
		}
	
		if(Rnode==0)
		{
			
			Rwhere<-which(obj$where==which(as.numeric(row.names(obj$frame))==nodes[i]*2+1))
			
		}
		
		
		
		
		# We now have all elements needed to calculate the SSEs and means
		if(length(Lwhere)==1)
		{
			moy<-obj$y[Lwhere,]
			MOY<-moy
			SSE[compt_SSE,]<-(moy-MOY)^2
			MOYs[compt_SSE,]<-moy
			compt_SSE<-compt_SSE+1
		}
		if(length(Lwhere)>1)
		{
			moy<-apply(obj$y[Lwhere,],MARGIN=2,FUN=mean)
			MOY<-rep(moy,times=nrow(obj$y[Lwhere,]))
			MOY<-matrix(MOY,nrow(obj$y[Lwhere,]),length(moy),byrow=TRUE)
			SSE[compt_SSE,]<-colSums((obj$y[Lwhere,]-MOY)^2)
			MOYs[compt_SSE,]<-moy
			compt_SSE<-compt_SSE+1
		}
		
		
		if(length(Rwhere)==1)
		{
			moy<-obj$y[Rwhere,]
			MOY<-moy
			SSE[compt_SSE,]<-(moy-MOY)^2
			MOYs[compt_SSE,]<-moy
			compt_SSE<-compt_SSE+1
		}
		if(length(Rwhere)>1)
		{
			moy<-apply(obj$y[Rwhere,],MARGIN=2,FUN=mean)
			MOY<-rep(moy,times=nrow(obj$y[Rwhere,]))
			MOY<-matrix(MOY,nrow(obj$y[Rwhere,]),length(moy),byrow=TRUE)
			SSE[compt_SSE,]<-colSums((obj$y[Rwhere,]-MOY)^2)
			MOYs[compt_SSE,]<-moy
			compt_SSE<-compt_SSE+1
		}
		
		# Keep the L and R where for futur use (i.e. objects present at each node/child)
		LWHERE[[i]]<-Lwhere
		RWHERE[[i]]<-Rwhere
		
		
	}	
	
	# Now that we have all SSEs... let's calculate the contribution to R2 of each species, for each node.
	
	# Build the node/leaf for each level matrix
	
	Tree<-mat.or.vec(length(nodes),length(as.numeric(row.names(obj$frame))))
	rownames(Tree)<-as.character(nodes)
	colnames(Tree)<-row.names(obj$frame)
	
	for(i in 1:length(nodes))
	{
		# For every node, which ones are nodes and leafs #
		
		# Repeat previous line (nodes r nodes, leafs r leaf)
		if(i!=1) Tree[i,]=Tree[i-1,]
		
		# Update : This node is a node
		Tree[i,which(as.numeric(row.names(obj$frame))==nodes[i])]<-1
		
		# Update : It's leafs are leafs
		Tree[i,which(as.numeric(row.names(obj$frame))==nodes[i]*2)]<-2
		Tree[i,which(as.numeric(row.names(obj$frame))==nodes[i]*2+1)]<-2
		
		# 
		
	}
	
	# Calculation of contribution #
	
	SST<-sum(SSE[1,])
	
	R2<-mat.or.vec(length(nodes),ncol(obj$y))
	
	for(i in 1:length(nodes))
	{
		
		## leafs before ##
		if(i!=1)
		{
			bef_pos<-which(Tree[i-1,]==2)
			leafs<-1
			for(j in 1:length(bef_pos))
			{
				
				leaf<-as.numeric(colnames(Tree)[bef_pos[j]])
			
				if(leaf%%2==0)
				{
					# if it is a left node
					# where is it in SSE ?
					leafs<-cbind(leafs,which(nodes==leaf/2)*2)
				}
			
				if(leaf%%2==1)
				{
					# if it is a right node
					# where is it in SSE ?
					leafs<-cbind(leafs,which(nodes==(leaf-1)/2)*2+1)
				}
			}
			
			if(length(leafs)!=2)
			{
				
				bef<-apply(SSE[leafs[-1],],FUN=sum,MARGIN=2) # sum of leaf
			}
			if(length(leafs)==2)
			{
				
				bef<-SSE[leafs[-1],] # sum of leaf
			}
		}
		if(i==1) bef<-SSE[1,]
		
		
		## leafs now ##
		now_pos<-which(Tree[i,]==2)
		leafs<-1
		for(j in 1:length(now_pos))
		{
			leaf<-as.numeric(colnames(Tree)[now_pos[j]])
			
			if(leaf%%2==0)
			{
				# if it is a left node
				# where is it in SSE ?
				leafs<-cbind(leafs,which(nodes==leaf/2)*2)
			}
			
			if(leaf%%2==1)
			{
				# if it is a left node
				# where is it in SSE ?
				leafs<-cbind(leafs,which(nodes==(leaf-1)/2)*2+1)
			}
		}
		
		if(length(leafs)!=2)
		{
			
			now<-apply(SSE[leafs[-1],],FUN=sum,MARGIN=2) # sum of leaf
		}
		if(length(leafs)==2)
		{
			
			now<-SSE[leafs[-1],] # sum of leaf
		}
		
		## R2 for this node ##
		R2[i,]<-(bef-now)/SST
	}
	
	# Now get the discriminant species #
	
	
	sumR2<-apply(R2, 1,sum)

	mat_sum<-rep(sumR2, ncol(R2))
	mat_sum<-matrix(mat_sum,nrow(R2),ncol(R2))
	pourct<-(R2/mat_sum)*100
	colnames(pourct)<-colnames(obj$y)
	
	
	
	# Give nodes in real order of R2
	
	#R2_perct<-R2*100
	#R2_pernode<-apply(R2_perct,MARGIN=1,FUN=sum)
	#ii<-order(R2_pernode)
	#nodes<-as.numeric(row.names(obj$frame)[splits][ii])
	

	# Return results needed for complete summary
	
	
	MOYs<-MOYs[-1,]
	
	# Build table as table 1 in Dea'th (2002) article
	
	# Species names
	
	ifelse(length(species)==0,col1names<-colnames(obj$y),col1names<-species)

	coli<-t(R2)
	
	col_total_tree<-rowSums(coli)
	col_total_species<-SSE[1,]/sum(SSE[1,])
	
	#TABLE1
	TABLE1<-cbind(coli,col_total_tree,col_total_species)
	
	rownames(TABLE1)<-col1names
	
	# Add sum of columns
	total_col<-colSums(TABLE1)
	names(total_col)<-'SUMS'
	
	TABLE1<-rbind(TABLE1,total_col)
	
	ii<-ORDER_NODES(obj,R2,splits)
	if(LABELS==FALSE) mat_labels_NOTINORDER<-mat_labels_fct(obj)
	if(!LABELS==FALSE) mat_labels_NOTINORDER<-LABELS
	mat_labels<-mat_labels_NOTINORDER[ii,]
	if(length(splits)>1) colnames(TABLE1)[1:(ncol(TABLE1)-2)]=sedit(mat_labels[,2],c('<','>','=',' '),c('','','',''))
	
	if(length(splits)==1) colnames(TABLE1)[1:(ncol(TABLE1)-2)]=sedit(mat_labels[2],c('<','>','=',' '),c('','','',''))
	
	TABLE1<-TABLE1*100
	
	# nŽcessaire pour plot :
	# typeplot,Cex,widthtree,heighttree,widthtriordi,heighttriordi,test.F=test.F,silent=silent
	
	res<-list(nodes,pourct,R2,obj,percent,MOYs,RWHERE,LWHERE,TABLE1,LABELS,mat_labels)
	
	
	names(res)<-c('nodes','pourct','R2','obj','percent','MOYs','RWHERE','LWHERE','TABLE1','LABELS','mat_labels')
	class(res)<-'MRT'
	return(res)
	

	
}

### ------------------------------------ ###
###              plotmvpart              ###
### Gives usual tree, TRIplot and BIplot ###
### ------------------------------------ ###

plot.MRT<-function(x,NodeMarking=TRUE,typeplot=c('tree'),Cex=0.5,widthtree=7, heighttree=9,R2A=FALSE,X=NULL,T=NULL,...)
#ADD test.F=FALSE,silent=TRUE,'triordi',widthtriordi=7, heighttriordi=9,
{
	if(R2A & is.null(X)) stop ('If R2A is TRUE, X must be provided (non NULL)')
	if(R2A & is.null(T)) stop ('If R2A is TRUE, T must be provided (non NULL)')
	# If we want the tree plot
	if(any(typeplot=='tree'))
	{
		plot_tree(x$obj,NodeMarking,x$R2,Cex,widthtree,heighttree,x$LABELS,R2A=R2A,X=X,T=T,...)
	}

	# If we want the ordination plot
#	if(any(typeplot=='triordi'))
#	{
#		
#		plot_triplot(x$obj,x$LWHERE,x$RWHERE,x$R2,Cex,widthtriordi,heighttriordi,test.F=test.F,silent=silent,x$LABELS)
#	}
	
}

plot_tree<-function(obj,NodeMarking,R2,Cex,widthtree,heighttree,LABELS,R2A=FALSE,X=NULL,T=NULL,...)
{

	#require(Hmisc)

	# points in node
	a=3
	# points between nodes
	b=2
	
	# nodes designated by numbers following the frame imposed by mvpart
	nodes<-as.numeric(rownames(obj$frame))
	
	maxs=max(nodes)
	
	# RElates to how long the x axis should be (don't remember exaclty what this is)
	# Beleive this is how many maximum splits there is from the root node (leaf counting as a node)
	s<-floor(log(maxs)/log(2))+1
	# Splits : which one are not leafs, in others words their position in the table obj$frame #
	splits<-which(obj$frame[,1]!="<leaf>")

	# length of x axis is (a+b)s+b   3 |2|  3  |2|  3
	#par(mar=c(5,5,4,2),cex.axis=1,col.axis=184)
	
	# Getting the real order of the nodes (correction for children explaining more variance than parents)
	ii<-ORDER_NODES(obj,R2,splits)
	if(LABELS==FALSE) mat_labels_NOTINORDER<-mat_labels_fct(obj)
	if(!LABELS==FALSE) mat_labels_NOTINORDER<-LABELS
	R2_node<-rowSums(R2)
	
	
	mat_labels<-mat_labels_NOTINORDER[ii,]
	
	
	# Cumulative complexity (cumulative R square of the tree from top to bottom)
	#cum_R<-cumsum(R2_node[ii[order(length(ii):1)]])[order(length(ii):1)][order(ii)]
	
	#cum_R<-cumsum(R2_node[order(-ii)])[order(-c(1:length(R2_node)))][ii]
	
	#cum_R<-cumsum(R2_node[ii])
	
	cum_R<-c(0,cumsum(R2_node[ii]))
	# If we want to show to up to the maximum value of R2
	#quartz(width=widthtree, height=heighttree)
	par(cex=Cex)#,mar=c(5,5,0,2))
	plot(c(0,(a+b)*2^(s-2)+b),c(max(cum_R*100),0),xlim=c(0,(a+b)*2^(s-2)+b),ylim=c(max(cum_R*100),0),xlab="",ylab=expression(R^2),type='n',xaxt='n', yaxt='n',bty='n',plt=c(0.9,0.9,0.9,0.9))

		
	
	# If we want to show the whole axis no matter what ??? not shure about this
	#else
	#plot(c(0,(a+b)*2^(s-2)+b),c(100,0),xlab="",ylab=expression(R^2),type='n',xaxt='n', yaxt='n')
	
	#axis(side=2, at = c(0,c(0,cum_R*100)), labels = as.character(round(c(0,cumsum(R2_node[ii[order(length(ii):1)]][c(length(R2_node):1)]),1)*100,digits=2))[(length(R2_node)+2):1],pos=-0.75,las=1)
	
	axis(side=2, at = cum_R*100, labels = as.character(round(cum_R*100,digits=2)),pos=-0.5,las=1)
	
	# Buils mat_labels and mat_done for pretty labels attibution #
	# See if thisis done before #
	# call matlabels #
	
	#print(length(R2_node))
	
	
	# Get rid of data.matrix once and for all
	if(length(R2_node)==1) mat_labels<-as.data.frame(t(mat_labels))
	mat_labels[,2]<-substring(mat_labels[,2],substring.location(mat_labels[1,2], ').')$last+1)
	mat_labels[,3]<-substring(mat_labels[,3],substring.location(mat_labels[1,3], ').')$last+1)
	
	# Build vertical lines at the R2 level , and build horiontal lines #
	
	compt<-1
	complexity<-obj$frame$complexity
	
	# What if we change the complexity vector right away ?
	
	for(i in 1:length(nodes))
	{
		if(obj$frame$var[i]!='<leaf>') 
		{
			complexity[i]<-R2_node[compt]
			compt<-compt+1
		}	
	}
	
	
	compt<-1
	for(i in 1:length(nodes))
	{	
		# If the node is not a leaf
		if(obj$frame$var[i]!='<leaf>')
		{
			
			# This node is at which level of the hierarchy ?
			si<-floor(log(nodes[i])/log(2))+1
			# At this level, which position from the left does it occupy ?
			posi<-nodes[i]-2^(si-1)+1
			# Change complexity for the right one
			#complexity[i]<-R2_node[compt]
			
			
			# Horizontale lines, R2 printed, size of node if leaf
			
			if(si==(s-1))
			{
				min<-(posi-1)*(a+b)+b
				max<-posi*(a+b)
				lines(c(min,max),c(cum_R[which(ii==compt)]*100,cum_R[ which(ii==compt)]*100))
				#points(-1,obj$frame$complexity[i],pch=23)
				
				# R under the node
				middle<-mean(c(min,max))
				
				aaa<-as.numeric(R2_node[compt])*100
				aaa<-format(aaa,digits=2,nsmall=2)
				
				pretty_middle<-paste('(',nodes[i],')\n',aaa,'%',sep='')
				text(middle,cum_R[which(ii==compt)]*100,pos=1,labels=pretty_middle,cex=1)
					
				# Here we add the labels on the node ?
				caract_min<-mat_labels[which(as.numeric(mat_labels[,1])==nodes[i]),2] #which(as.numeric(splits[,1])==nodes[i])
				caract_max<-mat_labels[which(as.numeric(mat_labels[,1])==nodes[i]),3]
				text(min,cum_R[which(ii==compt)]*100-0.2,caract_min,pos=2,cex=1,offset=0.5)
				text(max,cum_R[which(ii==compt)]*100-0.2,caract_max,pos=4,cex=1,offset=0.5)
				compt=compt+1
			}

			if(si==(s-2))
			{
				min<-(2*posi-1)*(a+b)-(1/2)*a
				max<-2*posi*(a+b)-(1/2)*a
				lines(c(min,max),c(cum_R[which(ii==compt)]*100,cum_R[which(ii==compt)]*100))
				# points(-1,obj$frame$complexity[i],pch=23)
				
				middle<-mean(c(min,max))
				
				aaa<-as.numeric(R2_node[compt]*100)
				aaa<-format(aaa,digits=2,nsmall=2)
				
				pretty_middle<-paste('(',nodes[i],')\n',aaa,'%',sep='')
				text(middle,cum_R[which(ii==compt)]*100,pos=1,labels=pretty_middle,cex=1)
				#text(-1,obj$frame$complexity[i],pos=4,labels=format(obj$frame$complexity[i],digits=3),cex=0.5)
				caract_min<-mat_labels[which(as.numeric(mat_labels[,1])==nodes[i]),2]
				caract_max<-mat_labels[which(as.numeric(mat_labels[,1])==nodes[i]),3]
				text(min,cum_R[which(ii==compt)]*100-0.2,caract_min,pos=2,cex=1,offset=0.5)
				text(max,cum_R[which(ii==compt)]*100-0.2,caract_max,pos=4,cex=1,offset=0.5)
				compt=compt+1
			}

			if(si<(s-2))
			{
				# DiffÃ©rence de 2^(s-si-3)
				min<-(2^(s-si-2)+(posi-1)*2^(s-si-1)-2^(s-si-3))*(a+b)+(1/2)*b
				max<-(2^(s-si-2)+(posi-1)*2^(s-si-1)+2^(s-si-3))*(a+b)+(1/2)*b
				# Vertical line
				lines(c(min,max),c(cum_R[which(ii==compt)]*100,cum_R[which(ii==compt)]*100))
				# points(-1,obj$frame$complexity[i],pch=23)
				
				middle<-mean(c(min,max))
				
				aaa<-as.numeric(R2_node[compt])*100
				aaa<-format(aaa,digits=2,nsmall=2)
				
				pretty_middle<-paste('(',nodes[i],')\n',aaa,'%',sep='')
				text(middle,cum_R[which(ii==compt)]*100,pos=1,labels=pretty_middle,cex=1)
				#text(-1,obj$frame$complexity[i],pos=4,labels=format(obj$frame$complexity[i],digits=3),cex=0.5)
				caract_min<-mat_labels[which(as.numeric(mat_labels[,1])==nodes[i]),2]
				caract_max<-mat_labels[which(as.numeric(mat_labels[,1])==nodes[i]),3]
				text(min,cum_R[which(ii==compt)]*100-0.2,caract_min,pos=2,cex=1,offset=0.5)
				text(max,cum_R[which(ii==compt)]*100-0.2,caract_max,pos=4,cex=1,offset=0.5)
				compt=compt+1
			}
			
			
			# And the verticals #
			
			# Find the level of R2 of the children, if it's not a leaf
			
			# Number of the child to the right and the child to the left
			child_right<-which(as.numeric(row.names(obj$frame))==as.numeric(row.names(obj$frame)[i])*2)
			child_left<-which(as.numeric(row.names(obj$frame))==as.numeric(row.names(obj$frame)[i])*2+1)
			
			R_child_right<-complexity[child_right]
			R_child_left<-complexity[child_left]
			
			# These are the complexity of the children, if they are not leafs.
			posR_child_right<-match(R_child_right,R2_node[ii],nomatch=999)
			posR_child_left<-match(R_child_left,R2_node[ii],nomatch=999)
			
		
			
			if(posR_child_right==999) R_cum_child_right=max(cum_R)
			if(posR_child_right!=999) R_cum_child_right<-cum_R[posR_child_right]
			if(posR_child_left==999) R_cum_child_left=max(cum_R)
			if(posR_child_left!=999) R_cum_child_left<-cum_R[posR_child_left]
			
			
			if(NodeMarking)
			{
				if(obj$frame[child_right,1]=='<leaf>') text(min,max(cum_R)*100,paste('n=',obj$frame[child_right,2],'\n(#',child_right,')',sep=''),cex=0.75,pos=1)
			
				if(obj$frame[child_left,1]=='<leaf>') text(max,max(cum_R)*100,paste('n=',obj$frame[child_left,2],'\n(#',child_left,')',sep=''),cex=0.75,pos=1)
			}
			
			if(!NodeMarking)
			{
				if(obj$frame[child_right,1]=='<leaf>') text(min,max(cum_R)*100,paste('n=',obj$frame[child_right,2],sep=''),cex=0.75,pos=1)
			
				if(obj$frame[child_left,1]=='<leaf>') text(max,max(cum_R)*100,paste('n=',obj$frame[child_left,2],sep=''),cex=0.75,pos=1)
			}
			
			
			
			# Draw both vertical lines
			# For right child
		
			lines(c(min,min),c(cum_R[which(ii==(compt-1))]*100,R_cum_child_right*100))
			# For left child
			lines(c(max,max),c(cum_R[which(ii==(compt-1))]*100,R_cum_child_left*100))
			
		}
		
		# If the node is a leaf #
		# if(obj$frame$var[i]=='<leaf>')
		# {
			# Add number of objects per node
			# This leaf is at which level of the hierarchy ?
		#	si<-floor(log(nodes[i])/log(2))+1
			# At this level, which position from the left does it occupy ?
		#	posi<-nodes[i]-2^(si-1)+1
			
		#	text(,0,paste('n= '.arbre$frame[i,2]),cex=1,offset=0.5,pos=1)
		#}
	}
	
	
	# Add adjusted R squared and other interesting values #
	
	len <- dim(obj$cptable)[1]
	exp1<- expression(R^2)
	bbb<-signif(1-obj$cptable[len, 3],digits=3)
	if(R2A)
	{
		R2adj<-R2aGDF(MRT(obj),T=T,X=X,tau_const=0.6,...)
		foot0 <- paste('R2a: ',round(R2adj*100,2),'%')
		foot1 <- paste('R2 : ',bbb*100,'%')
		foot2 <- paste("  Error : ", signif(obj$cptable[len, 3], 
                digits=3))
		foot <- paste(foot0,foot1,foot2, "  CV Error : ", signif(obj$cptable[len, 
                  4], digits=3), "  SE : ", signif(obj$cptable[len, 
                  5], digits=3))
         mtext(foot, side = 1, line = 3.5, cex = 0.85)
    }
    
    if(!R2A)
	{
		foot1 <- paste('R2 : ',bbb*100,'%')
		foot2 <- paste("  Error : ", signif(obj$cptable[len, 3], 
                digits=3))
		foot <- paste(foot1,foot2, "  CV Error : ", signif(obj$cptable[len, 
                  4], digits=3), "  SE : ", signif(obj$cptable[len, 
                  5], digits=3))
         mtext(foot, side = 1, line = 3.5, cex = 0.85)
    }
            
   
	
	
}


### -------------------------------- ###
###         mat_lables_fct           ###
### -------------------------------- ###

# Le nŽcessaire pour construire 'mat_labels'
mat_labels_fct<-function(obj)
{
	split<-which(obj$frame$var!='<leaf>')
	
	Labels<-labels(obj,pretty=pretty)[-1]
	
	# Les splits en ordre ...
	splits<-as.numeric(row.names(obj$frame)[split])
	
	# Ajout des colonnes d'enfants de gauche et droite
	D_child<-2*splits+1
	G_child<-2*splits
	
	mat_splits<-cbind(splits,D_child,G_child)
	
	# Get the imaginary children out
	for(i in 1:nrow(mat_splits))
	{
		for (j in 2:3)
		{
			if(!any(splits==mat_splits[i,j]))
			{
				mat_splits[i,j]=0
			}
		}
	}


	# Count the number of relatives on the left
	G_relatives<-mat.or.vec(nrow(mat_splits),nrow(mat_splits))
	mat_splits_relatives<-cbind(mat_splits,G_relatives)
	
	
	
	for(i in 1:nrow(mat_splits)) # Pour tout les noeuds
	{
		temp<-mat_splits_relatives[i,3] # the first left child
		for(j in 1:nrow(mat_splits)) # Pour tout les noeuds qui pourraient tre parents
		{
			if(any(mat_splits_relatives[j,1]==temp)) # On retrouve ce noeuds dans les parents de gauche
			{
				mat_splits_relatives[j,i+3]<-1
				temp<-c(temp,mat_splits_relatives[j,c(2,3)])
				# retire les 0
				temp<-unique(temp)
				
			}
		}
	}
	

	# How many much left relatives ? :o)
	
	if(nrow(mat_splits_relatives)==1) 
	{
		nrelatives=0
	}
	else 
	{
		nrelatives<-apply(mat_splits_relatives[,-c(1:3)],FUN=sum,MARGIN=2)
	}
	# Caculating the positions of the labels
	
	# Initialisation
	pos_labels<-mat.or.vec(length(Labels),2)
	
	compt<-1
	for(i in 1:nrow(mat_splits)) # Pour tout les noeuds
	{
	
		# First empty space
		zero<-which(pos_labels[,1]==0)[1]
		
		# Left child : the first empty space
		pos_labels[zero,1]<-compt
		compt<-compt+1
		pos_labels[zero,2]<-1
		
		
		# Right child : pos Left child + 2xnb relatives on left + 1
		pos_labels[zero+2*nrelatives[i]+1,1]<-compt
		compt<-compt+1
		pos_labels[zero+2*nrelatives[i]+1,2]<-2
		
	}
	
	# Buiding labels.. the real ones
	Labels_ordre<-Labels[order(pos_labels[,1])]
	
	# matrix()
	Labels<-matrix(Labels_ordre,nrow(mat_splits),2,byrow=TRUE)
	Labels<-cbind(splits, Labels)
	return(Labels)
}

# ------------------------------------------------------------ #
# this function makes the summary output for the MRT analsyis. #
# ------------------------------------------------------------ #
	
summary.MRT<-function(object,IndvalPART=TRUE,IndvalNODE=TRUE,...)
{
	# IndvalPART : add indval species for final partition
	# IndvalNODE : add indval species for each node
	
	if(object$LABELS==FALSE) mat_labels<-mat_labels_fct(object$obj)
	if(!object$LABELS==FALSE) mat_labels<-object$LABELS
	

	cat('Portion (%) of deviance explained by species for every particular node','\n','\n')
	
	for(i in 1:length(object$nodes))
	{
		cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
		# Print node number, it's complexity, it's discriminant species, it's objects
		
		# Node number + complexity
		cat('                    --- Node',object$nodes[i],'---\n','Complexity(R2)',sum(object$R2[i,])*100,'\n',mat_labels[i,-1],'\n\n')
		
		# Discriminant species #
		cat('~ Discriminant species :','\n')
		
		# Create table of results #
		table<-rbind(object$pourct[i,object$pourct[i,]>object$percent],object$MOYs[i*2-1,object$pourct[i,]>object$percent],object$MOYs[i*2,object$pourct[i,]>object$percent])
		rownames(table)<-c('% of expl. deviance','Mean on the left','Mean on the right')
	
		print(table)
		
		# Build table of indval if 
		if(IndvalNODE)
		{
			# INDVAL species #
			
			LWHEREnode<-object$LWHERE[[i]]
			RWHEREnode<-object$RWHERE[[i]]
			Ynode<-object$obj$y[c(LWHEREnode,RWHEREnode),]
		
			
		clustnode<-c(mat.or.vec(length(LWHEREnode),1)+1,mat.or.vec(length(RWHEREnode),1)+2)

			
			INDVALnode<-indval(Ynode,clustering=clustnode,numitr=1000)

			
			cat('\n','~ INDVAL species for this node: : left is 1, right is 2','\n',sep='')
			summary(INDVALnode, p=0.05, type='short', ...)
		}
		
		cat('\n')
		
		cat('\n')
	}
	
	# Objects in each leaf
	cat('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
	cat('               --- Final partition ','---\n','\n')
	for(i in 1:length(object$obj$frame$var))
	{
		if(object$obj$frame$var[i]=='<leaf>')
		{
			cat('Sites in the leaf #',i,'\n')
			print(row.names(object$obj$y)[object$obj$where==i])
			cat('\n','\n')
		
		}
	}
	
	# Build table of indval if 
	if(IndvalPART)
	{
		# INDVAL species #
		cat('~ INDVAL species of the final partition:','\n')

		INDVALpart<-indval(object$obj$y,clustering=object$obj$where,numitr=1000)
		summary(INDVALpart, p=0.05, type='short',...)
		
		# Correspondance between $where and numbers given by indval
		k=1
		cat('\n')
		cat('~ Corresponding MRT leaf number for Indval cluster number:','\n')
		ordered_where<-sort(unique(object$obj$where))
		for(j in 1:length(object$obj$frame$var))
		{ 
			
			if(object$obj$frame$var[j]=='<leaf>')
			{
				cat('MRT leaf #',ordered_where[k],'is Indval cluster no.',k,'\n')
				k=k+1
		
			}
		}
	}
}

# This function takes Lwhere (list of left objects in every node) and Rwhere (list of rigth object in every node), and the labels from mat_labels, to create what is described in the function.


code.tree <- function(obj,LWHERE,RWHERE,LABELS)
# This function codes the tree.
# Each coding veriable in matrix 'code' represents one of the splits of the MRT.
# The objects on the left and right of the split have weights corresponding to 
# the inverse of the number of observations in each group.
# The objects not concerned with a split have a weight of zero.
# As a consequence, each column sums to 0.
#
#	   Marie-HŽlne Ouellette, August 2009
#      Modified from Pierre Legendre, August 2009
#      
{
	
	n <- nrow(obj$y)

	ifelse(length(LWHERE)==1,nvar<-1,nvar <- nrow(LABELS))
	
	code <- matrix(0,n,nvar)

	
	# For all splits, create code

	for(i in 1:length(LWHERE))
	{
		vec1 <- LWHERE[[i]]
		
		n1 <- length(vec1)
		vec2 <- RWHERE[[i]]
		
		n2 <- length(vec2)
		code[vec1,i] <- 1/n1
		
		code[vec2,i] <- -1/n2
		
	}

	# Cut the crap out (< > =) of LABELS
	if(length(LWHERE)==1)LABELS<-as.data.frame(t(LABELS))
	colnames(code) = sedit(LABELS[,2],c('<','>','=',' '),c('','','',''))
	return(code)
}

# ----------------------- #
# This function takes a coded tree from code.tree and builds the TRiplot as suggested by PL
# ----------------------- #

#plot_triplot<-function(obj,LWHERE,RWHERE,R2,Cex,widthtriordi,heighttriordi,test.F=test.F,silent=silent,LABELS,...)
#{
#	require(rdaTest)
#	group<-obj$where
#	splits<-which(obj$frame[,1]!="<leaf>")
#	ii<-ORDER_NODES(obj,R2,splits)
#	if(LABELS==FALSE) LABELS_NOTINORDER<-mat_labels_fct(obj)
#	if(!LABELS==FALSE) LABELS_NOTINORDER<-LABELS
#	LABELS<-LABELS_NOTINORDER[ii,]
	
#	#dummy = model.matrix(~ group)       # group4 was eliminated by model.matrix
#	#grouplast = dummy[,1] - apply(dummy[,-1],1,sum)
#	#dummy2 = cbind(grouplast,dummy[,-1])   # Put back group4 into the matrix
	
#	coded.tree = code.tree(obj,LWHERE,RWHERE,LABELS)
#	apply(coded.tree,2,sum)  # All columns sum to 0
	
#	# Change row names to site membership for plot
#	rownames(obj$y)<-obj$where
#	rownames(coded.tree)<-obj$where
#	rda.res3 = rdaTest(obj$y, coded.tree,test.F=test.F,silent=silent)
	
#	# Plot without objects
#	plot(rda.res3, graph.type="Z",cex=Cex,width=widthtriordi, height=heighttriordi,scaling=1)

	# Bimultivariate redundancy statistic (canonical R-square):
	# R-square = 0.7856865;   adjusted R-square =  0.7484146
#}

# ----------------------- #
# This function takes a coded tree from code.tree and builds the BIplot as suggested by PL
# ----------------------- #

# This function returns ii, the order of the nodes in terms of R2 of splits
ORDER_NODES<-function(obj,R2,splits)
{
	ii<-1
	if(length(splits)>1)
	
	{
		# The order shall be given by the children matrix (correction for when children explain more variance than parent)
	
		parent<-1 # the
		nodes2<-1
		pos_nodes<-1
		pos_parent<-1
		
		R2_node<-rowSums(R2)
		nodes_name<-as.numeric(row.names(obj$frame))[splits]
	
		for(i in 1:(length(splits)-1))
		{
			# For this node, which children are nodes also (update nodes with children) ?
			if(any(as.numeric(row.names(obj$frame))[splits]==parent*2))
			{
				nodes2<-c(nodes2,parent*2)
				pos_nodes<-c(pos_nodes,which(nodes_name==parent*2))
			}
		
			if(any(as.numeric(row.names(obj$frame))[splits]==(parent*2+1)))
			{
				nodes2<-c(nodes2,parent*2+1)
				pos_nodes<-c(pos_nodes,which(nodes_name==(parent*2+1)))
			}
		
			# Get this parent (node) out of nodes and pos_nodes
			pos_nodes<-pos_nodes[-which(nodes2==parent)]
			pos_nodes
			nodes2<-nodes2[-which(nodes2==parent)]
			nodes2

		
			# Get position in R2 and splits of the nodes
		
			# Of these, which has the smallest R2 ? This is our next parent
			parent<-nodes_name[which(R2_node==max(R2_node[pos_nodes]))]
			parent
			pos_parent<-which(nodes_name==parent)
			pos_parent
		
			# Stock it's position as it is the next node
		
			ii<-c(ii,which(R2_node==max(R2_node[pos_nodes])))
	
		}
	}
return(ii)
}
