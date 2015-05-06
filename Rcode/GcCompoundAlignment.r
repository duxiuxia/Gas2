GcCompoundAlignment<-function(params)
{
	# parameters for feature detection
	library(tcltk)
	JobName<-params$JobName
	workingDir <- params$WorkDir
	Max_Compound_ET_shift<-as.numeric(params$Max_Compound_ET_shift)/60
	nDaughters<-as.numeric(params$nFragments)
	DataFilelist <- params$ProcessedFiles#get file list
	maxSpecAngle<-params$maxSpecAngle
	featureSets<-vector("list",length(DataFilelist))
	
	#get feature dataset
	for(fileindex in 1:length(DataFilelist))
	{
		inFilePath <- DataFilelist[[fileindex]]#get the full path of each file
		Features<-read.csv(inFilePath)
		colnames(Features)[1]<-"FeatureID"
		Features<-subset(Features,CompoundID>0)#remove the low intensity feature as noise
		featureSets[[fileindex]]<-cbind(Features,FileID=fileindex)#add file label to feature
	}
	featureSets<-do.call(rbind,featureSets)
	
	featureSets$mz<-featureSets$FeatureID#get mz from featureID
	mzVec<-unlist(strsplit(as.character(featureSets$mz)," "))#remove feature ID from mz
	mzVec<-mzVec[seq(from=1,to=length(mzVec),by=2)]
	featureSets$mz<-mzVec
	featureSets$FeatureID<-NULL#remove featuer ID
	
	featureSets$clustID<-0#init cluster id to each feature

	#init feature index
	clust_ind<-0  
	curSet<-subset(featureSets,clustID==0)
	nWindowCount<-0
	total <- nrow(curSet)# create progress bar
	pb <- tkProgressBar(title = paste(JobName,"alignment"), min = 0,
			                    max = total, width = 300)

	#iteritatively group feature apexes according to mz and et range
	while(nrow(curSet)>0)
	{
		nrow(curSet)
		
		curCompoundID<-curSet[1,]$CompoundID# get current reference scan
		curFileID<-curSet[1,]$FileID
		curComFeatures<-subset(curSet,CompoundID==curCompoundID&FileID==curFileID)
		curET<-mean(curComFeatures$ET)
		curProfile<-subset(curComFeatures,select=c("mz","int"))
		curProfile<-aggregate(curProfile$int,by=list(curProfile$mz),max)
		colnames(curProfile)[1]<-"mz"
		
		#remember feature index and mzProfile ID which contain current scan
		Featureind<-rownames(subset(featureSets,FileID==curFileID&CompoundID%in%curCompoundID))
		isMatched<-0
		
		#if(nrow(curProfile)>=nDaughters)
#		{
			fidRange<-1:length(DataFilelist)
			fidRange<-fidRange[-curFileID]#exclude the current file
		
			for(fid in fidRange)
			{
				#get scans from each file which is close to current compound ET
				curClust<-subset(curSet,FileID==fid)
				curClust<-subset(curClust,abs(ET-curET)<Max_Compound_ET_shift)
				AnalyizedCompoundIDs<-unique(curClust$CompoundID)
				
				candidates<-NULL
				for(i in AnalyizedCompoundIDs)
				{

					refProfile<-subset(curClust,CompoundID==i,select=c("mz","int"))
					refProfile<-aggregate(refProfile$int,by=list(refProfile$mz),max)
					colnames(refProfile)[1]<-"mz"
					colnames(refProfile)[2]<-"y"
					mzMergeSet<-merge(curProfile,refProfile,by.x="mz",by.y="mz",all.x=TRUE,all.y=TRUE)
					
					x<-mzMergeSet$x
					y<-mzMergeSet$y
					if(length(x[!is.na(x)])>=nDaughters&&length(y[!is.na(y)])>=nDaughters)
					{
						if(length(x[is.na(x)])>0)
						{x[is.na(x)]<-0}
						if(length(y[is.na(y)])>0)
						{y[is.na(y)]<-0}
						#normalize
						if(sum(x)>0)
						{x<-x/sum(x)}
						if(sum(y)>0)
						{y<-y/sum(y)}
						
						#distCorr<-abs(cor(x,y))
						#distEucl<-dist(t(cbind(x,y)))
						dotProd<-(x%*%y)/((x%*%x)^0.5*(y%*%y)^0.5)
						angle<-acos(as.numeric(dotProd[1]))*180/pi
						
						candidates<-rbind(candidates,c(ET=i,angle=angle))
					}
					#browser()	
				}
					
				if(is.null(candidates)==FALSE)
				{		#browser()					
					minDist<-min(candidates[,2])
					minCompoundID<-candidates[which.min(candidates[,2]),1]
					#maxETs<-AnalyizedETs[abs(AnalyizedETs-maxET)<=daughterShift]
				
					if(minDist<maxSpecAngle)
					{
						isMatched<-1
						ind<-rownames(subset(featureSets,FileID==fid&CompoundID==minCompoundID))
						Featureind<-c(Featureind,ind)
					}
				}			
			}
		#}
		
		#if(isMatched==1)
#		{
			clust_ind<-clust_ind+1#increase clustid
			#clustMarker<-clust_ind
		#}else
#		{
#			clustMarker<--1
#			
#			}
		
		#if(length(which(Featureind==3213))>0){browser()}
		
		#mark featureSets and mzProfiles with clust ID
		featureSets$clustID[rownames(featureSets)%in%Featureind]<-clust_ind
		
		curSet<-subset(featureSets,clustID==0)#update the curset
		setTkProgressBar(pb, total-nrow(curSet), label=paste( round((total-nrow(curSet))/total*100, 0),  "% done"))
	}

	#featureSets<-subset(featureSets,clustID>0)
	total<-clust_ind
	for(i in 1:clust_ind)
	{
		setTkProgressBar(pb, total-i, label=paste( round((total-i)/total*100, 0),  "% done"))
		logicInd<-featureSets$clustID==i
		featureSets$alignedET[logicInd]<-mean(featureSets$ET[logicInd])#set cluster mark		
	}
	
	featureSets<-featureSets[order(featureSets$alignedET),]
	write.csv(featureSets,file=paste(workingDir,"/Data/",JobName,"_aligned_Features.csv",sep=""))
	close(pb)#close progress bar
		#write.csv(mzProfiles,file=paste(workingDir,"/Data/",JobName,"_aligned_mzProfiles.csv",sep=""))
}

GcUniqueIonTable<-function(params)
{
	library(tcltk)
	library(snow)
	JobName<-params$JobName
	workingDir <- params$WorkDir

	#DataFilelist <- params$ProcessedFiles#get file list
	#nFile<-length(DataFilelist)

	#featureSets<-read.csv(file=paste(workingDir,"/Data/",JobName,"_aligned_Features.csv",sep=""))
	featureSets<-read.csv(file=params$AlignedFile)
	FileIDs<-unique(featureSets$FileID)
	nFile<-length(FileIDs)
	fileVsClust<-subset(featureSets,select=c("FileID","clustID"))#get common mother Ions which occur in 80% of datasets
	fileVsClust<-unique(fileVsClust)
	groupCount<-aggregate(fileVsClust$FileID,by=list(fileVsClust$clustID),length)#group count fileid by clustid
	nMajority<-round(as.integer(params$occuRatio)*nFile)
	CommonClustID<-subset(groupCount,x>=nMajority)$Group.1# get common clustid which occurs in over 4 datasets
	featureSets<-subset(featureSets,clustID%in%CommonClustID)#filter featureset by common clustid

	nCommonClustID <- length(CommonClustID)	


	#pb <- tkProgressBar(title = paste(JobName,"Unique Mass searching"), min = 0,
#			                    max = total, width = 300)
	c1<-makeCluster(as.integer(params$nNode),type=params$clustType)
	
	commonClustGroups<-clusterSplit(c1,CommonClustID)
	
	
	cat("step1:construct spectrum accross multiple datasets...\n")

	transFeatureSet<-function(commonClustGroup,FileIDs,featureSets,params)
	{
		nFile<-length(FileIDs)
		nMajority<-round(as.integer(params$occuRatio)*nFile)
		nClust<-length(commonClustGroup)
		localTrainingSet<-vector("list",nClust)
		for(ind in 1:nClust)
		{
			curClustID<-commonClustGroup[ind]
			curClust<-subset(featureSets,clustID==curClustID)
			mzVec<-sort(unique(curClust$mz))
			nMzVec<-length(mzVec)
			localTrainingSet[[ind]]<-vector("list",nMzVec)
			#check each fragment 
			for(i in 1:nMzVec)
			{
				curMz<-mzVec[i]
				subDs<-subset(curClust,mz==curMz,select=c("FileID","int"))
				subDs<-merge(subDs,FileIDs,by.y=1,by.x=1,all.y=TRUE)
				nNonZeroValue<-length(subDs$int[is.na(subDs$int)==FALSE])
				#if current fragment exists in most datasets, then add it to spectrum
				if(nNonZeroValue>=0)#nMajority) 
				{
				subDs$int[is.na(subDs$int)]<-0
				subDs<-aggregate(subDs$int,by=list(subDs$FileID),sum)#merge the same fragments from the same file
				#transpose the intensity 
				rownames(subDs)<-subDs$FileID
				subDs<-subset(subDs,select=c("x"))
				subDs<-as.data.frame(t(subDs))
				subDs$mz<-curMz
				localTrainingSet[[ind]][[i]]<-cbind(subDs,clustID=curClustID)
				}
			}
			tempTrainingSet<-do.call(rbind,localTrainingSet[[ind]])
			if(is.null(tempTrainingSet)==FALSE) #check if the result is NULL
				{
					if(nrow(tempTrainingSet)>as.integer(params$nFragments))
					{localTrainingSet[[ind]]<-tempTrainingSet #if satisfy the threshold, then assign the result set
					}else{
						localTrainingSet[[ind]]<-numeric() #if less number of fragments,then assign empty value 
						}
				}else{

				localTrainingSet[[ind]]<-numeric() #if NULL then assign the empty value
			}
					
					
		}
		localTrainingSet<-do.call(rbind,localTrainingSet)
	}
	trainingSet<-clusterApplyLB(c1,commonClustGroups,transFeatureSet,FileIDs,featureSets,params)
	trainingSet<-do.call(rbind,trainingSet)
	
	cat("step2:get common daughter ions among datasets for each compouond...\n")
	getCommonFrag<-function(commonClustGroup,FileIDs,featureSets)
	{	

		nClust<-length(commonClustGroup)
		commonMzResults<-vector("list",nClust)
		for(ind in 1:nClust)
		{
			curClustID<-commonClustGroup[ind]
			isFirstFile=TRUE
			for(fileindex in FileIDs)
			{
				curClust<-subset(featureSets,clustID==curClustID&FileID==fileindex,select=c("alignedET","clustID","mz"))
				if(nrow(curClust)>0)
				{
					if(isFirstFile)
					{
						commonMz<-curClust
					}else
					{
						commonMz<-merge(commonMz,curClust,by.x="mz",by.y="mz")
						isFirstFile=FALSE
					}
				}
				
			}
			
			commonMzResults[[ind]]<-commonMz			
		}
		commonMzResults<-do.call(rbind,commonMzResults)
		commonMzResults
	}
	commonMzs<-clusterApplyLB(c1,commonClustGroups,getCommonFrag,FileIDs,featureSets)
	commonMzs<-do.call(rbind,commonMzs)
		
	cat("step3:#indentify unique daughter ions by comparing with neigbour mother profiles...\n")
	uMassWindow<-as.integer(params$uMassWindow)/60
	trainingSet$isUnique<-0
	trainingSet$id<-1:nrow(trainingSet)
	trainingSet$isMaxUmass<-0
	searchUmass<-function(commonClustGroup,commonMzs,uMassWindow,FileIDs,trainingSet)
	{

		nClust<-length(commonClustGroup)
		UniqueIdList<-vector("list",nClust)
		MaxUniqueIDList<-vector("list",nClust)
		for(ind in 1:nClust)
		{
			curClustID<-commonClustGroup[ind]
			Aclust<-subset(commonMzs,clustID==curClustID)#get current cluster
			curET<-unique(Aclust$alignedET)
			#get neighbor compounds within uMassWindow
			neighborCompounds<-subset(commonMzs,clustID!=curClustID&abs(curET-alignedET)<uMassWindow)
			localClustIDs<-unique(neighborCompounds$clustID)
			for(cmClustID in localClustIDs)
			{
				Bclust<-subset(neighborCompounds,clustID==cmClustID)
				#compare two mothers within the range
				
				excludeMz<-merge(Aclust,Bclust,by.x="mz",by.y="mz")$mz
				Aclust<-Aclust[!(Aclust$mz%in%excludeMz),]
				if(nrow(Aclust)==0){break}
			
			}
			if(nrow(Aclust)==0){next}
			#get id in bioMarkerList table
			curTrainingSet<-merge(trainingSet,Aclust)
			UniqueIdList[[ind]]<-curTrainingSet$id
			
			#get max int accross datasets for each ion
			if(nrow(curTrainingSet)>0)
			{
				curTrainingSet$rowMedian<-0
				for(i in 1:nrow(curTrainingSet))
				{
					nonZeroInd<-curTrainingSet[i,as.character(FileIDs)]>0
					nonZeroInt<-curTrainingSet[i,as.character(FileIDs)][nonZeroInd]
					curTrainingSet[i,]$rowMedian<-median(nonZeroInt)
				}
				MaxUniqueIDList[[ind]]<-curTrainingSet[which.max(curTrainingSet$rowMedian),]$id		#mark max umass in trainingSet
			}

		}
		UniqueIdList<-unlist(UniqueIdList)
		MaxUniqueIDList<-unlist(MaxUniqueIDList)
		UmassList<-list(MaxUniqueIDList=MaxUniqueIDList,UniqueIdList=UniqueIdList)
	}
	UmassLists<-clusterApplyLB(c1,commonClustGroups,searchUmass,commonMzs,uMassWindow,FileIDs,trainingSet)
	
	cat("collect and merge the results...\n")
	nUmassList<-length(UmassLists)
	MaxUniqueIDLists<-vector("list",nUmassList)
	UniqueIdLists<-vector("list",nUmassList)
	for(i in 1:nUmassList)
	{
		MaxUniqueIDLists[[i]]<-UmassLists[[i]]$MaxUniqueIDList
		UniqueIdLists[[i]]<-UmassLists[[i]]$UniqueIdList
		
	}
	MaxUniqueIDLists<-unlist(MaxUniqueIDLists)
	UniqueIdLists<-unlist(UniqueIdLists)
	
	tInd<-trainingSet$id%in%UniqueIdLists# mark umass
	trainingSet$isUnique[tInd]<-1
	tInd<-trainingSet$id%in%MaxUniqueIDLists# mark max umass
	trainingSet[tInd,]$isMaxUmass=1
	
	
	
	#add ET information to it
	bioMarkerList<-unique(subset(featureSets,select=c("clustID","alignedET","mz")))
	bioMarkerList<-merge(bioMarkerList,trainingSet)
	NonZeroMedian<-function(bioMarkerListRow)
	{
		bioMarkerListRow<-bioMarkerListRow[bioMarkerListRow>0]
		median(bioMarkerListRow)
	}
	bioMarkerList$int<-apply(subset(bioMarkerList,select=as.character(FileIDs)),1,NonZeroMedian)
	bioMarkerList<-subset(bioMarkerList,select=c("clustID","alignedET","mz","int","isUnique","isMaxUmass",as.character(FileIDs)))
	bioMarkerList<-bioMarkerList[order(bioMarkerList$alignedET),]
	write.csv(bioMarkerList,file=paste(workingDir,"/output/",JobName,"_bioMarkerList.csv",sep=""),row.names=FALSE)
	stopCluster(c1)

	
	cat("done!\n")
	
}
GcAlignmentSpectrumViz<-function(params)
{
	# parameters for feature detection
	JobName<-params$JobName
	workingDir <- params$WorkDir
	DataFilelist <- params$DataFiles#get file list
	
	library(playwith)
	featureSets<-read.csv(file=paste(workingDir,"/Data/",JobName,"_aligned_Features.csv",sep=""))
	
	ScanIDStart<-min(featureSets$ET)
	ScanIDEnd<-max(featureSets$ET)
	mzStart<-min(featureSets$mz)
	mzEnd<-max(featureSets$mz)
	IntStart<-min(featureSets$int)
	IntEnd<-max(featureSets$int)

	#feature map before alignment
	#playwith({
#	plot(0,type="n",xlab="mz",ylab="int",ylim=c(0,100), xlim=c(mzStart,mzEnd),main="spectrum before alignment")
#	for(fileindex in 1:15)#1:length(DataFilelist))
#	{
		#for(i in 1:max(featureSets$clustID))
#		for(i in c(866))
#		{
#			Features<-subset(featureSets,FileID==fileindex&clustID==i,select=c("int","mz"))
#			Features$int<-100*Features$int/max(Features$int)
#			points(Features$mz, Features$int,type="h",col=rainbow(15)[fileindex])
#		}
#	}
#	
#	},data.points=subset(featureSets,select=c("mz","int")),labels=featureSets$mz,cex.lab=0.02,new=TRUE)
#
		playwith({
			Features1<-subset(featureSets,FileID==4&clustID==958,select=c("int","mz"))
			Features1$int<-100*Features1$int/max(Features1$int)
			plot(Features1$mz, Features1$int,type="h",ylim=c(-100,100), xlim=c(mzStart,400),
			xlab=alignedET,ylab="int",main=paste("file:",fileindex," clust:",i,sep=""),col="red")
			Features2<-subset(featureSets,FileID==15&clustID==960,select=c("int","mz"))
			Features2$int<-100*Features2$int/max(Features2$int)
			points(Features2$mz, -Features2$int,type="h",ylim=c(-100,100), xlim=c(mzStart,400),col="blue")
	
		})


	clustIDs<-c(958,960)
	Clust<-subset(featureSets,clustID%in%clustIDs,select=c("alignedET","clustID","FileID","int","mz"))
	x11()
	fileIDs<-unique(Clust$FileID)
	nlength<-length(fileIDs)
	par(mfrow=c(ceiling(nlength/2),2))
	colorVec<-c("blue","red")
	for(i in 1:length(clustIDs))
	{
		curClustID<-clustIDs[i]
		curClust<-subset(Clust,clustID==curClustID,select=c("alignedET","FileID","int","mz"))
		alignedET<-unique(curClust$alignedET)
		curfileIDs<-unique(curClust$FileID)
		for(fileindex in curfileIDs)#1:length(DataFilelist))
		{
			Features<-subset(curClust,FileID==fileindex,select=c("int","mz"))
			Features$int<-100*Features$int/max(Features$int)
			plot(Features$mz, Features$int,type="h",ylim=c(0,100), xlim=c(mzStart,mzEnd),
			xlab=alignedET,ylab="int",main=paste("file:",fileindex," clust:",curClustID,sep=""),col=colorVec[i])
		}
	}
	


	# mz vs ET map for mother Ions which include the common daughter Ions in 40% of datasets
	#bioMarkerList<-read.csv(file=paste(workingDir,"/Data/",JobName,"_bioMarkerList.csv",sep=""))
#	bioMarkerList<-subset(bioMarkerList,is.na(clustID)==FALSE)
#	ScanIDStart<-min(bioMarkerList$alignedET)
#	ScanIDEnd<-max(bioMarkerList$alignedET)
#	mzStart<-min(bioMarkerList$mz)
#	mzEnd<-max(bioMarkerList$mz)
#	playwith({
#	plot(0,type="n",xlab="ET",ylab="mz",main="mz vs ET map for mother Ions after alignment",xlim=c(ScanIDStart,ScanIDEnd), ylim=c(mzStart,mzEnd))
#	for(i in unique(bioMarkerList$clustID))
#	{
#		daughters<-subset(bioMarkerList,clustID==i,select=c("alignedET","mz"))
#		points(daughters,pch=i,cex=0.4)
#		uMasses<-subset(bioMarkerList,clustID==i&isUnique==1,select=c("alignedET","mz"))
#		if(nrow(uMasses)>0)
#		{
#			points(uMasses,pch=i,col="green",cex=0.4)
#			maxUMass<-subset(bioMarkerList,clustID==i&isMaxUmass==1,select=c("alignedET","mz"))
#			points(maxUMass,pch=i,col="red",cex=0.5)
#			text(maxUMass$alignedET+0.05,maxUMass$mz,labels=maxUMass$avgMz,cex=1)
#		}
#	}
#	},data.points=subset(bioMarkerList,select=c("alignedET","mz")),labels=bioMarkerList$mz,cex.lab=0.02,new=TRUE)
#
#	
}



GcAlignmentViz<-function(params)
{
	# parameters for feature detection
	JobName<-params$JobName
	workingDir <- params$WorkDir
	DataFilelist <- params$DataFiles#get file list
	
	library(playwith)
	featureSets<-read.csv(file=paste(workingDir,"/Data/",JobName,"_aligned_Features.csv",sep=""))
	
	ScanIDStart<-min(featureSets$ET)
	ScanIDEnd<-max(featureSets$ET)
	mzStart<-min(featureSets$mz)
	mzEnd<-max(featureSets$mz)
	IntStart<-min(featureSets$int)
	IntEnd<-max(featureSets$int)

	#feature map before alignment
	#playwith({
#	plot(0,type="n",xlab="ScanID",ylab="mz",xlim=c(ScanIDStart,ScanIDEnd), ylim=c(mzStart,mzEnd),main="Mz vs ScanID before alignment")
#	for(fileindex in 5:5)#1:length(DataFilelist))
#	{
#		for(i in 1:max(featureSets$clustID))
#		{
#			Features<-subset(featureSets,FileID==fileindex&clustID==i,select=c("ET","avgMz"))
#			points(Features$ET, Features$avgMz,pch=i,col=rainbow(14)[fileindex])
#		}
#	}
#	
#	},data.points=subset(featureSets,select=c("ET","avgMz")),labels=featureSets$avgMz,cex.lab=0.02,new=TRUE)
#
	# mz vs ET map for mother Ions which include the common daughter Ions in 40% of datasets
	bioMarkerList<-read.csv(file=paste(workingDir,"/Data/",JobName,"_bioMarkerList.csv",sep=""))
	bioMarkerList<-subset(bioMarkerList,is.na(clustID)==FALSE)
	ScanIDStart<-min(bioMarkerList$alignedET)
	ScanIDEnd<-max(bioMarkerList$alignedET)
	mzStart<-min(bioMarkerList$mz)
	mzEnd<-max(bioMarkerList$mz)
	playwith({
	plot(0,type="n",xlab="ET",ylab="mz",main="mz vs ET map for mother Ions after alignment",xlim=c(ScanIDStart,ScanIDEnd), ylim=c(mzStart,mzEnd))
	for(i in unique(bioMarkerList$clustID))
	{
		daughters<-subset(bioMarkerList,clustID==i,select=c("alignedET","mz"))
		points(daughters,pch=i,cex=0.4)
		uMasses<-subset(bioMarkerList,clustID==i&isUnique==1,select=c("alignedET","mz"))
		if(nrow(uMasses)>0)
		{
			points(uMasses,pch=i,col="green",cex=0.4)
			maxUMass<-subset(bioMarkerList,clustID==i&isMaxUmass==1,select=c("alignedET","mz"))
			points(maxUMass,pch=i,col="red",cex=0.5)
			text(maxUMass$alignedET+0.05,maxUMass$mz,labels=maxUMass$avgMz,cex=1)
		}
	}
	},data.points=subset(bioMarkerList,select=c("alignedET","mz")),labels=bioMarkerList$mz,cex.lab=0.02,new=TRUE)

	
}

