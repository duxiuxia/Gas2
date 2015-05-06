
GcUniqueIonTableSQLite<-function(params)
{
	library(RSQLite)
	library(snow)
	JobName<-params$JobName
	workingDir <- params$WorkDir
	paraDB<- params$dbFile
	m <- dbDriver("SQLite")
	con <- dbConnect(m, dbname = paraDB)
	tname<-paste("t",JobName,"aligned",sep="_")
	if(!dbExistsTable(con,tname))
		stop(paste("Table '",tname,"' does not exist!",sep=""))
	featureSets<-dbReadTable(con,tname)
	######change intensity column name to int keep consisitant with the following reference to int
	names(featureSets)[names(featureSets)=="intensity"]<-"int"

	
	
	FileIDs<-unique(featureSets$FileID)
	nFile<-length(FileIDs)
	fileVsClust<-subset(featureSets,select=c("FileID","clustID"))#get common mother Ions which occur in 80% of datasets
	fileVsClust<-unique(fileVsClust)
	groupCount<-aggregate(fileVsClust$FileID,by=list(fileVsClust$clustID),length)#group count fileid by clustid
	nMajority<-round(as.numeric(params$sampleRatio)*nFile)
	CommonClustID<-subset(groupCount,x>=nMajority)$Group.1# get common clustid which occurs in over 4 datasets
	featureSets<-subset(featureSets,clustID%in%CommonClustID)#filter featureset by common clustid

	nCommonClustID <- length(CommonClustID)	


	#pb <- tkProgressBar(title = paste(JobName,"Unique Mass searching"), min = 0,
#			                    max = total, width = 300)
	c1<-makeCluster(as.numeric(params$nNode),type=params$clustType)
	
	commonClustGroups<-clusterSplit(c1,CommonClustID)
	
	
	cat("step1:construct spectrum accross multiple datasets...\n")

	transFeatureSet<-function(commonClustGroup,FileIDs,featureSets,params)
	{
		nFile<-length(FileIDs)
		nMajority<-round(as.numeric(params$sampleRatio)*nFile)
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
				if(nNonZeroValue>=0) 
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
					if(nrow(tempTrainingSet)>as.integer(params$Fragments))
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
	uMassWindow<-as.numeric(params$uMassWindow)/60
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
	
	dbWriteTable(con,"bioMarkerList",bioMarkerList,row.names=F,overwrite=T)

	# clean up
   	dbDisconnect(con)
	write.csv(bioMarkerList,file=paste(workingDir,"/output/",JobName,"_bioMarkerList.csv",sep=""),row.names=FALSE)
	stopCluster(c1)
	cat("done!\n")
	
}




