decomposition_finalAug<-function(inFilePath,params,cl,isDistParallel,clustingType)
{
	source(paste(params$codeDir,"pipeline.r",sep="/"))	
	library(gdata)
	library(gtools)
	library(cluster)
	library("ncdf")	
	options(expressions=1e5)
	DataFilelist<-params$DataFiles
	delaytime<-params$delaytime
	ScanInterval<-params$ScanInterval
	WorkDir<-params$WorkDir
	
	##load library
	lib<-readMSP2Spec(filename="../StdsLibrary/KQC.txt",withRT=T)
	
	fileName<-parseFileName(inFilePath)
	
	###############
	#read TIC data
	################		
	TICfile<-paste(WorkDir,"output/TIC/denoised_",fileName,"_TIC.cdf",sep="") 
	# The denoised TIC are saved in output/TIC/
	#if there's no denoised TIC file, then read raw cdf data
	TICfile<-ifelse(file.exists(TICfile),TICfile,inFilePath)
	ncid <- open.ncdf(TICfile)
	TIC <- get.var.ncdf(ncid, varid="total_intensity")
	close.ncdf(ncid)
	remove(ncid)
	
	
	###############
	#read EIC data
	################
	rawEICFile<-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
	denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
	EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
	
	if(file.exists(EICFile)==FALSE)next
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	##orignial mz vector from EIC data
	vectorMz<-get.var.ncdf(ncid, varid="mzVec")
	close.ncdf(ncid)
	remove(ncid)	
	
	###############################################
	#get decon window from TIC peak picking results
	###############################################
	cat(fileName,"reading TIC peak data...\n")
	TICPeakFile<-paste(params$WorkDir,"output/peakpicking/",fileName,"_TIC_PeakList.csv",sep="")
	TICApexList <- read.csv(TICPeakFile)
	TICApexList$RT<-(TICApexList$pkInd-1)*ScanInterval+delaytime
	TICpeaklist<-subset(TICApexList,isApex==1)
	
	winIDs<-1:nrow(TICpeaklist)
	
	#############################
	#read EIC peak picking result
	#############################
	cat(fileName,"reading EIC peak data...\n")
	EICPeakFile<-paste(WorkDir,"/output/peakpicking/",fileName,"PeakList.csv",sep="")
	if(file.exists(EICPeakFile)==FALSE)next
	EICpeaklist <- read.csv(EICPeakFile)
	colnames(EICpeaklist)[which(colnames(EICpeaklist)=="curMz")]<-"mz"
	
	#########################
	#start decon
	#########################
	cat(fileName,"start deconvolution...\n")
	
	###################
	#sequential version
	####################
	componentResults<-NULL
	specResults<-NULL
	mdlPkIDVec<-NULL
	EICpeaklist$flag=0
	totalscan<-length(vecInt)/length(vectorMz)
	
	TICPkmean <- mean(TIC[TICpeaklist$pkInd])
	minPkHeight<-TICPkmean*0.05
	
	for(windowID in winIDs)
	{
		cat("window:",windowID,"\n")
		curTICPk<-TICpeaklist[windowID,]	
		minET=(curTICPk$lboundInd-1)*params$ScanInterval+params$delaytime
		maxET=(curTICPk$rboundInd-1)*params$ScanInterval+params$delaytime
		
		##decide the exlcuded ion for current TIC window
		##based on ET, for model peak selection
		if(minET>=8)
		{
			CurNonUmassVec<-c(params$NonUmassVec,51:100)
		}else
		{
			CurNonUmassVec<-params$NonUmassVec
		}
		
		##################################################
		#get all the EIC peaks within current decon window
		##################################################
		localTICApexList<-subset(TICApexList,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		#get the potential compound number for TIC peak apex list
		#there should be at least as many or more compounds than the number of apex 
		#minCmp<-nrow(localTICApexList)
		
		#get all EIC peaks within current decon window
		localEICPeakList<-subset(EICpeaklist,flag==0&pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		
		if(nrow(localEICPeakList)>0)
		{
#			minPkHeight <- mean(TIC[localTICApexList$pkInd])*0.025
#			#minPkHeight <- ifelse(minPkHeight <=1000,1000,minPkHeight)
#			#minPkHeight <- ifelse(minPkHeight >=100000,0.1*minPkHeight,minPkHeight)
#			minPkHeight <- ifelse(minPkHeight <=800,800,minPkHeight)
#			minPkHeight <- ifelse(minPkHeight >5000&&minPkHeight <=15000,5000+round((minPkHeight-5000)*0.5),minPkHeight)
#			minPkHeight <- ifelse(minPkHeight >10000,10000+round((minPkHeight-10000)*0.1),minPkHeight)
#			
			##Exclude those mass not abundent enough for models
			CandidatePeaklist<-subset(localEICPeakList,Intensity>minPkHeight)
			
			##Exclude those mass not eligible for models
			CandidatePeaklist<-subset(CandidatePeaklist,!mz%in%CurNonUmassVec)
			
			##keep the unique mass profile of local EIC list
			##original local EIC list consists of merged peaks
			localEICPeakList <- subset(localEICPeakList,isApex==1)
			allprofiles<-vector("list",nrow(localEICPeakList))
			for(i in 1:nrow(localEICPeakList))
			{
				##get intensity vector of current mz
				##from long intensity vector
				curEICpk<-localEICPeakList[i,]			
				startInd<-curEICpk$offset+curEICpk$lboundInd
				endInd<-curEICpk$offset+curEICpk$rboundInd
				
				##get EIC peak profile
				ET<-curEICpk$lboundInd:curEICpk$rboundInd
				curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd])
				allprofiles[[i]]<-curProfile
				names(allprofiles)[i]<-rownames(curEICpk)	
			}
			
			##remove potentially shared ions by GSS and Sharpness
			if (minPkHeight <=5000) 
			{
				gaussianCutoff<-2.6
			}else 
			{
				gaussianCutoff<-2
			}
			
			#gaussianCutoff<-2.6
			ShrpCutoff<-5
			goodShapePeaklist<-subset(CandidatePeaklist,gss<=gaussianCutoff&shrp>=ShrpCutoff)
			
#			if (nrow(goodShapePeaklist)>10)
#			{
#				importMasslist <- goodShapePeaklist[order(goodShapePeaklist[,"Intensity"],decreasing=T),][c(1:5),]
#				ImportantLowMass <- unique(subset(importMasslist,gss<=2)$mz)
#				CurNonUmassVec<- CurNonUmassVec[!CurNonUmassVec%in%ImportantLowMass]
#				# remove other low and non-unique mass
#				goodShapePeaklist<-subset(goodShapePeaklist,!mz%in%CurNonUmassVec)
#			}
#			
			
			nPks<-nrow(goodShapePeaklist)		
			modelPkList<-NULL
			isMultGroups<-FALSE
			mdlPkIDs<-NULL
			
			if(nPks>0)
			{
				##get good peak profiles
				profiles<-vector("list",nrow(goodShapePeaklist))
				for(i in 1:nrow(goodShapePeaklist))
				{
					##get intensity vector of current mz
					##from long intensity vector
					curEICpk<-goodShapePeaklist[i,]					
					startInd<-curEICpk$offset+curEICpk$lboundInd
					endInd<-curEICpk$offset+curEICpk$rboundInd
					
					##get EIC profile
					ET<-curEICpk$lboundInd:curEICpk$rboundInd
					curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd])
					profiles[[i]]<-curProfile
					names(profiles)[i]<-rownames(curEICpk)	
				}
				
				##update each profile by creating mirror images
				##chosse the one with better GSS
				for(i in 1:length(profiles))
				{
					curPkId<-names(profiles[i])
					curPk<-goodShapePeaklist[curPkId,]
					pos<-curPk$gssPos
					profiles[[i]]<-copyShape(p=profiles[[i]],curPk$pkInd,from=pos)
					goodShapePeaklist[curPkId,]$lboundInd<-profiles[[i]]$ET[1]
					goodShapePeaklist[curPkId,]$rboundInd<-profiles[[i]]$ET[nrow(profiles[[i]])]
				}
				
				#############################
				#Hierachical clustering
				#decide the number of groups 
				#within current decon window
				#############################
				if(nPks>=2)
				{
					##save EIC profies in global variable for 
					##broadcasting to slave nodes for parallel computing
					assign("Global.curProfs", value=profiles, envir = .GlobalEnv)
					
					##calculate pairwise distance among the EICs				
					if(isDistParallel)
					{
						r<-parDistCal4(cl,isUnion=F)##parallel verison
					}else
					{
						r<-DistCal4(isUnion=F)#non-parallel version
					}	
					
					##convert distance matrix to triangle matrix
					distance <- as.dist(r)
					maxIntraDist<- 15 #cutoff to decide if split
					
					if(clustingType=="h")
					{
						##################################
						#hierarchical clustering and cut 
						#into groups by the height
						##################################
						clustResut<-hclust(distance)
						FinalClustResut<-cutree(clustResut,h=maxIntraDist)
						Clusters<-unique(FinalClustResut)
						if(length(Clusters)>1)isMultGroups<-T
						
					}else
					{	###############################################################
						#partitional clustering and split into groups by the silehoute
						###############################################################
						if(max(distance)<=maxIntraDist)
						{###if the max distance is already small enought, no need to deconvolute
							isMultGroups<-FALSE
						}else
						{
							
							if(nPks>2)
							{
								########################################
								#try different k to decide the cluster number by intra dist
								#######################################
								
								maxClustNumber<-min(20,nPks-1)
								asw<-numeric(maxClustNumber)
								for (k in min(maxClustNumber,max(2,minCmp)):maxClustNumber)
								{
									clustResut <-pam(distance,k)
									Clusters<-unique(clustResut$clustering)								
									isSplit<-F
									
									for(curCluster in Clusters)
									{
										#####get masses belong to current cluster
										curGroup<-names(clustResut$clustering[clustResut$clustering==curCluster])
										#######does not average the single-ion cluster which has zero intra dist
										if(length(curGroup)>1)
										{#######calculate the average intra distance to see if the current group is a valid single group of fragments
											
											if(max(as.dist(r[curGroup,curGroup]))>maxIntraDist)
											{##if current group is not valid, skip to the next K,
												#and the current k value 
												isSplit<-T
												break
											}
											
										}
									}
									if(!isSplit)
										asw[k]<-clustResut$silinfo $ avg.width
								}
								k.best <- which.max(asw)
								if(k.best>1)
								{
									FinalClustResut<-pam(distance,k.best,cluster.only=TRUE,stand=TRUE)
									
									#################################################
									#calculate the average intra dist in each cluster
									###############################################
									isTooBigIntraDist<-NULL
									avgIntraDistVec<-NULL
									Clusters<-unique(FinalClustResut)
									for(curCluster in Clusters)
									{
										#####get masses belong to current cluster
										curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
										
										IntraDist<-mean(as.dist(r[curGroup,curGroup]))
										isTooBigIntraDist<-c(isTooBigIntraDist,ifelse(IntraDist>6,T,F))
										if(length(curGroup)>1)#######does not average the single-ion cluster which has zero intra dist
											avgIntraDistVec<-c(avgIntraDistVec,IntraDist)
										
									}
									#################################################################
									#if after splitting,intra dist is smaller than before splitting
									#then accept the best clustering schema, 
									#and earch for unique mass within each cluster
									###############################################################
									if(mean(avgIntraDistVec)<mean(distance))
										isMultGroups<-TRUE
								}
							}else#if ==2, then directly set flag=T
							{
								isMultGroups<-TRUE
							}
						}
					}
				}
				
				
				
				
				######################
				#Model peak selection 
				######################
				
				modelPkList<-NULL
				
				# single group
				if (!isMultGroups && nPks >=1)
				{
					if (nPks==1) {
						modelPkList<-goodShapePeaklist
					} else 
					{
						##choose best one as model peak
						mdlPkCandidates<-goodShapePeaklist
						
						##scoring system
						#mass score
						mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
						
						##gassian similarity score
						mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
						
						##peak height score
						mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))					
						mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
						mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
						mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
						
						##linear combination of 3 scores
						scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
						#scores<-(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3+0.4*mdlPkCandidates$f4)
						
						##rescale score to (0,1000)
						mdlPkCandidates$score<-scale(scores,F,0.001)
						
						##select the max score (>500)
						mdlPkCandidates <- subset(mdlPkCandidates,score>=500)
						if (nrow(mdlPkCandidates)>0)
						{
							modelPkList<- mdlPkCandidates[which.max(mdlPkCandidates$score),]
						}
					}
				}
				
				
				if(isMultGroups)
				{
					if(nPks==2)
					{
						##only two candidates, simply assign them as model peaks
						#modelPkVec<-rownames(goodShapePeaklist)
						modelPkList<-goodShapePeaklist
					}else
					{
						nMinFragment<-1
						if(length(Clusters)>0)
						{	
							##select model peak from each group[!isTooBigIntraDist]
							for(curCluster in Clusters)
							{
								##get masses belong to current cluster
								curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
								
								##get model peak candidates 
								mdlPkCandidates<-goodShapePeaklist[curGroup,]
								
								if(nrow(mdlPkCandidates)>0)
								{
									##mass score
									mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
									##gassian similarity score
									mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
									##peak height score
									mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))								
									
									##assign zero value when each factor is NA
									mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
									mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
									mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
									
									##linear combination of 3 scores
									scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
									
									##rescale score to (0,1000)
									mdlPkCandidates$score<-scale(scores,F,0.001)
									
									##select the max score
									currentMdlPk<-mdlPkCandidates[which.max(mdlPkCandidates$score),]
									modelPkList<-rbind(modelPkList,currentMdlPk)
								}
							}
						}
					}			
				}						
			}
			
			
			###################################
			#decomposition based on model peaks
			###################################
			
			decom = 1
			while (decom==1) 
			{
				modelPkDistance <- NULL
				specList<-NULL
				componentList<-NULL
				mdlPkIDs <- NULL
				
				##No co-eluting compounds (only one or no model peak detected)
				##Then extract spectrum directly
				if(nrow(modelPkList)<=1||is.null(modelPkList))
				{
					##add model peak to the spectrum
					specList[[1]]<-data.frame(mz=as.integer(),int=as.integer())
					
					##No model peak detected,then simply extract 
					##spectrum from the TIC peak apex time
					if(nrow(modelPkList)==0||is.null(modelPkList))
					{	
						##get 1st mass				
						curMz<-as.character(localEICPeakList$mz[1])
						
						##get TIC peak apex time and boundaries
						maxInd<-curTICPk$pkInd
						curLbound<-curTICPk$lboundInd
						curRbound<-curTICPk$rboundInd
						
						##get the pk height,id, profile and area of the 1st peak
						curPhHeight<-localEICPeakList$Intensity[1]
						curPkID<-rownames(localEICPeakList)[1]
						curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.x="ET",by.y="ET",all.x=T)
						curProfile1$int[is.na(curProfile1$int)]<-0
						curArea<-sum(curProfile1$int)
						
					}else
					{
						##single one model peak detected,then extract 
						##spectrum from the model peak apex time
						mdlPkIDs<-rownames(modelPkList)
						curMz<-modelPkList$mz[1]
						maxInd<-modelPkList$pkInd[1]
						curLbound<-modelPkList$lboundInd[1]
						curRbound<-modelPkList$rboundInd[1]
						curPhHeight<-modelPkList$Intensity[1]	
						curPkID<-rownames(modelPkList)[1]
						
						# profiles save good peaks profiles
						curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
						curProfile1$int[is.na(curProfile1$int)]<-0
						curArea<-sum(curProfile1$int)
						
					}
					
					##Initialize component list 
					componentList[[1]]<-data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,pkInd=maxInd,Intensity=curPhHeight,area=curArea,windowID=windowID,compoundID=1)
					
					#create spec list as NIST format for future identifiation
					clustID<-format((maxInd-1)*ScanInterval+delaytime,digit=4)
					names(specList)[1]<-paste(clustID,curMz)
					
					##add rest peaks from the local window 
					##to component list and spectral list
					for(i in 1:nrow(localEICPeakList))
					{
						curPkID<-rownames(localEICPeakList)[i]					
						curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
						curProfile1$int[is.na(curProfile1$int)]<-0
						curArea<-sum(curProfile1$int)
						
						curMz<-as.character(localEICPeakList$mz[i])
						curInt<-subset(curProfile1,ET==maxInd)$int[1]
						
						if(!is.na(curInt)&&curInt>0)
						{
							specList[[1]]<-rbind(specList[[1]],data.frame(mz=as.integer(curMz),int=curInt))#add umass to the spectrum					
							componentList[[1]]<-rbind(componentList[[1]],data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,pkInd=maxInd,Intensity=curInt,area=curArea,windowID=windowID,compoundID=1))
						}
					}
					decom=0
				}else
				{
					##More than one model peaks
					##Decompostion
					
					#######################################
					#Makeup 1 for splitting issue
					#Checking pairwise model peaks distance
					#######################################
					
					modelPkDistance <- as.matrix(modelPkDist(modelPkList,profiles))
					rmmodelPk <- NULL
					ModelPkNames <- rownames(modelPkDistance)
					
					##model same matrix as modelPKDistance with index
					psmatrix <- matrix(c(1:length(modelPkDistance)),ncol=nrow(modelPkDistance),byrow=T)
					psmodel <- which(modelPkDistance<=4&modelPkDistance>0)
					
					##get index of model peaks with smaller distance
					for (i in psmodel) {
						pos <- matrix.index(psmatrix, i)
						rmmodelPk <- c(rmmodelPk,ModelPkNames[pos[1]],ModelPkNames[pos[2]])
					}
					rmmodelPks <- unique(rmmodelPk)
					
					rmmodelPkList <- modelPkList[rmmodelPks,]
					modelPkList1 <- rmmodelPkList[which.max(rmmodelPkList$score),]
					modelPkList <- modelPkList[!rownames(modelPkList)%in%rmmodelPks,]
					
					##Only keep the one with highest score
					modelPkList <- rbind(modelPkList,modelPkList1)				
					mdlPkIDs<-rownames(modelPkList)
					
					##only one model peak after testing
					##GO back to the upper step
					if (nrow(modelPkList)==1) decom=1
					
					else {
						##still multiple models
						lbound<-min(modelPkList$lboundInd)
						rbound<-max(modelPkList$rboundInd)
						
						#################################
						#get all EIC profiles within the window 
						#X matrix format:column-EIC,row-scan
						#############################################
						
						X<-NULL
						mzVec<-unique(localEICPeakList$mz)
						nMz<-length(mzVec)
						for(i in 1:nMz)
						{
							mz<-mzVec[i]
							mzInd<-which(vectorMz==mz)
							startInd<-(mzInd-1)*totalscan+1
							endInd<-startInd+totalscan-1
							curVecInt <- vecInt[startInd:endInd]
							X<-cbind(X,curVecInt[lbound:rbound])
						}
						colnames(X)<-mzVec	
						nScans<-nrow(X)
						
						###########################################
						#Get model peak profiles as S matrix
						#each column is the EIC for one model peak
						###########################################
						
						S<-NULL
						modelPkList$area<-0
						for(i in 1:nrow(modelPkList))
						{
							curProfile1<-merge(data.frame(ET=as.integer(lbound:rbound)),profiles[[rownames(modelPkList[i,])]],by.y="ET",all.x=T)
							
							##expand the profile by filling zero values to those uncovered area
							curProfile1$int[is.na(curProfile1$int)]<-0
							S<-cbind(S,curProfile1$int)
							modelPkList[i,]$area<-sum(curProfile1$int)		
							
							##init the component & spec list	
							componentList[[i]]<-data.frame(mz=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),pkInd=as.integer(),Intensity=as.integer(),area=as.integer(),windowID=as.integer(),compoundID=as.integer())
							
							##add model peak to the spectrum
							specList[[i]]<-data.frame(mz=as.integer(),int=as.integer())
							curRT<-format((modelPkList$pkInd[i]-1)*ScanInterval+delaytime,digit=4)
							names(specList)[i]<-paste(curRT,modelPkList$mz[i])
						}					
						colnames(S)<-mdlPkIDs
						
						##Residual minimization
						for(i in 1:ncol(X))
						{
							curMz<-colnames(X)[i]						
							M<-X[,i]#mixture signal
							
							srcIDs<-mdlPkIDs
							A<-optim(par=rep(0,ncol(S)),fn=fr, gr = NULL,M,S,lower = 0, method="L-BFGS-B")$par
							A<-as.matrix(A)
							
							##restore each ion signals and save them into 
							##component list ans speclist
							for(m in 1:length(srcIDs))
							{
								curmdlPkID<-srcIDs[m]
								j<-which(mdlPkIDs==curmdlPkID)
								
								if(A[m]>0)
								{
									specList[[j]]<-rbind(specList[[j]],data.frame(mz=as.integer(curMz),int=A[m,]))#add the curent restored mixing factor to spec
									componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=as.integer(modelPkList[curmdlPkID,]$Intensity*A[m,]),area=as.integer(modelPkList[curmdlPkID,]$area*A[m,]),windowID=windowID,compoundID=j))
								}
								
							}
						}
						
						###############################################
						#Makeup 2 for splitting issue
						#Checking pairwise resovlved spectra similarity
						###############################################
						
						RepeatModelPk <- NULL
						Num <- nrow(modelPkList)
						indexPairs<-combinations(Num,2)
						for(i in 1:nrow(indexPairs))
						{
							index <- indexPairs[i,]
							index1 <- index[1]
							index2 <- index[2]
							
							##calculate time difference between two model peaks
							modelPKET <- modelPkList$pkInd[c(index1,index2)]
							
							score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
							if (score >=970||(abs(modelPKET[1]-modelPKET[2]) <=6&&score >=900))
							{
								selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index2],mdlPkIDs[index1])
								RepeatModelPk <- c(RepeatModelPk,selectedModelPkID)
							}
						}
						
						##select unique filtered model peaks 
						RepeatModelPk <- unique(RepeatModelPk)
						modelPkList <- modelPkList[!rownames(modelPkList)%in%RepeatModelPk,]
						
						##The modelPk number is reduced, RE-decompose
						if (nrow(modelPkList) == Num) decom = 0				
					} # else model PK num is >2 after testing distance	
				} # model PK num is >2 at first
			} # while statement
			
			##flag analyzed EIC peaks 
			EICpeaklist[rownames(localEICPeakList),]$flag<-1
			componentResults<-c(componentResults,componentList)
			specResults<-c(specResults,specList)
			
			##model peak id
			mdlPkIDVec<-c(mdlPkIDVec,mdlPkIDs)		
		}	
	}	
	
	componentResults<-do.call(rbind,componentResults)	
	cat(fileName,"finished deconvolution...\n")
	
	##################################
	#output lib matching results
	###################################
	
	##similarity cutoff=600
	libmatch <- libMatching(lib,specResults,600)
	write.csv(libmatch,file=paste(params$WorkDir,"/output/decon/",fileName,"_libmatch.csv",sep=""))
	
	##############################################################
	#collect spectrum result in dataframe for reverse lib matching
	##############################################################
	
	##exp spectra as lib, get top 1 hit for each std
	evaluateIdentification(refSpec=lib,inputSpec=specResults,RT_Tolerance=70,minSpecSimilarity=650,fileName,withRT=T,params)
	cat(fileName,"finish matching...\n")
	write.csv(mdlPkIDVec,file=paste(params$WorkDir,"/output/decon/",fileName,"mdlpklist.csv",sep=""))
	
	################################################
	#collect component result in dataframe for alignment
	###################################################
	
	##filter out zero-int fragments
	componentResults<-subset(componentResults,Intensity>0,select=c("mz","lboundInd","rboundInd","pkInd","Intensity","area","windowID","compoundID"))
	write.csv(componentResults,row.names=FALSE,file=paste(params$WorkDir,"/output/decon/",fileName,"Decon.csv",sep=""))
	cat(fileName,"writing component results...\n")
}



decomposition_June<-function(inFilePath,params,cl,isDistParallel,clustingType)
{
	source(paste(params$codeDir,"pipeline.r",sep="/"))	
	library(gdata)
	library(gtools)
	library(cluster)
	library("ncdf")	
	options(expressions=1e5)
	DataFilelist<-params$DataFiles
	delaytime<-params$delaytime
	ScanInterval<-params$ScanInterval
	WorkDir<-params$WorkDir
	
	#####load library
	lib<-readMSP2Spec(filename="../StdsLibrary/KQC.txt",withRT=T)
	
	fileName<-parseFileName(inFilePath)
	
	###############
	#read TIC data
	################		
	TICfile<-paste(WorkDir,"output/TIC/denoised_",fileName,"_TIC.cdf",sep="") 
	TICfile<-ifelse(file.exists(TICfile),TICfile,inFilePath)
	ncid <- open.ncdf(TICfile)
	TIC <- get.var.ncdf(ncid, varid="total_intensity")
	close.ncdf(ncid)
	remove(ncid)
	
	
	###############
	#read EIC data
	################
	rawEICFile<-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
	denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
	EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
	
	if(file.exists(EICFile)==FALSE)next
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	#####orignial mz vector from EIC data
	vectorMz<-get.var.ncdf(ncid, varid="mzVec")
	close.ncdf(ncid)
	remove(ncid)	
	
	###########################################
	#get decon window from TIC peak picking results
	###########################################
	cat(fileName,"reading TIC peak data...\n")
	TICPeakFile<-paste(params$WorkDir,"output/peakpicking/",fileName,"_TIC_PeakList.csv",sep="")
	TICApexList <- read.csv(TICPeakFile)
	TICApexList$RT<-(TICApexList$pkInd-1)*ScanInterval+delaytime
	TICpeaklist<-subset(TICApexList,isApex==1)
	
	winIDs<-1:nrow(TICpeaklist)
	
	#########################
	#read EIC peak picking result
	#########################
	cat(fileName,"reading EIC peak data...\n")
	EICPeakFile<-paste(WorkDir,"/output/peakpicking/",fileName,"PeakList.csv",sep="")
	if(file.exists(EICPeakFile)==FALSE)next
	EICpeaklist <- read.csv(EICPeakFile)
	colnames(EICpeaklist)[which(colnames(EICpeaklist)=="curMz")]<-"mz"
	
	#########################
	#start decon
	#########################
	cat(fileName,"start deconvolution...\n")
	
	###################
	#sequential version
	####################
	componentResults<-NULL
	specResults<-NULL
	mdlPkIDVec<-NULL
	EICpeaklist$flag=0
	totalscan<-length(vecInt)/length(vectorMz)
	
	for(windowID in winIDs)
	{
		cat("window:",windowID,"\n")
		curTICPk<-TICpeaklist[windowID,]
		#decide the exlcuded ion list for current TIC window
		minET=(curTICPk$lboundInd-1)*params$ScanInterval+params$delaytime
		maxET=(curTICPk$rboundInd-1)*params$ScanInterval+params$delaytime
		if(minET>=8)
		{
			CurNonUmassVec<-c(params$NonUmassVec,51:100)
			#CurNonUmassVec<-c(params$NonUmassVec,51:72,74:100)
		}else
		{
			CurNonUmassVec<-params$NonUmassVec
		}
		
		####################################################
		#get all the EIC peaks within current decon window
		##################################################
		localTICApexList<-subset(TICApexList,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		localEICPeakList<-subset(EICpeaklist,flag==0&pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		
		if(nrow(localEICPeakList)>0)
		{	
			allprofiles<-vector("list",nrow(localEICPeakList))
			for(i in 1:nrow(localEICPeakList))
			{
				##############################################
				#get intensity vector the current mz
				# from the long  intensity vector
				################################################
				curEICpk<-localEICPeakList[i,]
				
				startInd<-curEICpk$offset+curEICpk$lboundInd
				endInd<-curEICpk$offset+curEICpk$rboundInd
				
				#######get the profile defined by EIC peak boundary
				ET<-curEICpk$lboundInd:curEICpk$rboundInd
				curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd])
				allprofiles[[i]]<-curProfile
				names(allprofiles)[i]<-rownames(curEICpk)
				
			}
			
			########filter out noisy small peaks
			minPkHeight <- mean(TIC[localTICApexList$pkInd])*0.025
			# removing possible noise
			minPkHeight <- ifelse(minPkHeight <=1000,1000,minPkHeight)
			# for a window with large and small peaks together, trying to keep the small one
			# avoiding the dominance of large peaks
			minPkHeight <- ifelse(minPkHeight >=100000,0.1*minPkHeight,minPkHeight)
			CandidatePeaklist<-subset(localEICPeakList,Intensity>minPkHeight)
			
			###exlude the shared ion by exclude list
			#CandidatePeaklist<-subset(CandidatePeaklist,!mz%in%CurNonUmassVec)
			
			#################plot the profiles of all peaks
#				colVec<-rainbow(nrow(CandidatePeaklist))
#				names(colVec)<-CandidatePeaklist$mz
#				labelVec<-paste("mz:",CandidatePeaklist$mz," ,gss:",format(CandidatePeaklist$gss,digits=2)," ,shrp:",format(CandidatePeaklist$shrp,digits=2))
#				dataPoints<-subset(CandidatePeaklist,select=c("pkInd","Intensity"))
#				dataPoints$pkInd<-(dataPoints$pkInd-1)*ScanInterval+delaytime
#				playwith({
#							plot(x=0,type = "n",xlim=(c(min(CandidatePeaklist$lboundInd),max(CandidatePeaklist$rboundInd))-1)*ScanInterval+delaytime,ylim=c(0,maxInt1))
#							for(i in 1:nrow(CandidatePeaklist))
#							{
#								curMz<-as.character(CandidatePeaklist[i,]$mz)
#								curID<-rownames(CandidatePeaklist[i,])
#								points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, type = "l",col=colVec[curMz])
#								
#							}
#						},data.points=dataPoints,labels=labelVec)
			##################	
			
			gaussianCutoff<-2.6
			ShrpCutoff<-5
			goodShapePeaklist<-subset(CandidatePeaklist,gss<=gaussianCutoff&shrp>=ShrpCutoff)
			
			if (nrow(goodShapePeaklist)>10)
			{
				ImportantLowMass <- unique(goodShapePeaklist[order(goodShapePeaklist[,"Intensity"],decreasing=T),][c(1:5),]$mz)
				CurNonUmassVec<- CurNonUmassVec[!CurNonUmassVec%in%ImportantLowMass]
				# remove other low and non-unique mass
				goodShapePeaklist<-subset(goodShapePeaklist,!mz%in%CurNonUmassVec)
			}
			
			nPks<-nrow(goodShapePeaklist)		
			modelPkList<-NULL
			isMultGroups<-FALSE
			mdlPkIDs<-NULL
			
			if(nPks>0)
			{
				#maxInt2<-max(goodShapePeaklist$Intensity)
				
				####get the subset of the profiles for the goodshape peaks
				profiles<-allprofiles[as.character(rownames(goodShapePeaklist))]
				
				####update each profile by duplicating the half of peak that has good shape
				
				for(i in 1:length(profiles))
				{
					curPkId<-names(profiles[i])
					curPk<-goodShapePeaklist[curPkId,]
					pos<-curPk$gssPos
					profiles[[i]]<-copyShape(p=profiles[[i]],from=pos)
					goodShapePeaklist[curPkId,]$lboundInd<-profiles[[i]]$ET[1]
					goodShapePeaklist[curPkId,]$rboundInd<-profiles[[i]]$ET[nrow(profiles[[i]])]
				}
				
				#####################plot the profiles of selected peaks with good gaussian shape
#					colVec<-rainbow(nrow(goodShapePeaklist))
#					labelVec<-paste("pkID:",rownames(goodShapePeaklist),"mz:",goodShapePeaklist$mz," ,gss:",format(goodShapePeaklist$gss,digits=2)," ,shrp:",format(goodShapePeaklist$shrp,digits=2))
#					dataPoints<-subset(goodShapePeaklist,select=c("pkInd","Intensity"))
#					dataPoints$pkInd<-(dataPoints$pkInd-1)*ScanInterval+delaytime
#					playwith({
#					plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,max(maxInt2)))
#					for(i in 1:nrow(goodShapePeaklist))
#					{
#						curID<-rownames(goodShapePeaklist[i,])
#						points((profiles[[curID]]$ET-1)*ScanInterval+delaytime,profiles[[curID]]$int, type = "l",col=colVec[i])
#						
#					}
#					},data.points=dataPoints,labels=labelVec)
				#####################
				
				
				######################################
				#decide the number of groups in the current window
				########################################
				if(nPks>=2)
				{
					#save the EIC profies in global variable for broadcasting to slave nodes for parallel computing
					assign("Global.curProfs", value=profiles, envir = .GlobalEnv)
					#			assign("Global.vecInt", value=vecInt, envir = .GlobalEnv)
					#calculate the mutual distance among the EICs
					
					if(isDistParallel)
					{
						r<-parDistCal4(cl,isUnion=F)##parallel verison
					}else
					{
						r<-DistCal4(isUnion=F)#non-parallel version
					}	
					#convert the distance matrix to triangle matrix
					distance <- as.dist(r)
					maxIntraDist<-15 #cutoff to decide if split
					
					if(clustingType=="h")
					{
						#######################################################################################
						#hierarchical clustering and cut into groups by the height
						#####################################################################################
						clustResut<-hclust(distance)
#									plot(clustResut)
						FinalClustResut<-cutree(clustResut,h=maxIntraDist)
						Clusters<-unique(FinalClustResut)
						if(length(Clusters)>1)isMultGroups<-T
						
						
						##########################plot the profiles of cluster groups
#							colVec<-rainbow(length(Clusters))
#							labelVec<-paste("pkID:",rownames(goodShapePeaklist),"mz:",goodShapePeaklist$mz," ,gss:",format(goodShapePeaklist$gss,digits=2)," ,shrp:",format(goodShapePeaklist$shrp,digits=2))
#							playwith({
#							plot(x=0,type = "n",xlim=(c(minLbound,maxRbound)-1)*ScanInterval+delaytime,ylim=c(0,max(maxInt2)))
#							for(i in 1:length(Clusters))
#							{
#								######get current cluster ID
#								curCluster<-Clusters[i]
#								#####get masses belong to current cluster
#								curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
#								
#								for(j in 1:length(curGroup))
#								{
#									
#									points((profiles[[curGroup[j]]]$ET-1)*ScanInterval+delaytime,profiles[[curGroup[j]]]$int, type = "l",col=colVec[i])
#								}
#							}
#							},data.points=subset(goodShapePeaklist,select=c("pkInd","Intensity")),labels=labelVec)
#							
						###########################							
						
						
					}else
					{	#######################################################################################
						#partitional clustering and split into groups by the silehoute
						#####################################################################################
						if(max(distance)<=maxIntraDist)
						{###if the max distance is already small enought, no need to deconvolute
							isMultGroups<-FALSE
						}else
						{
							
							if(nPks>2)
							{
								########################################
								#try different k to decide the cluster number by intra dist
								#######################################
								
								maxClustNumber<-min(20,nPks-1)
								asw<-numeric(maxClustNumber)
								for (k in min(maxClustNumber,max(2,minCmp)):maxClustNumber)
								{
									clustResut <-pam(distance,k)
									Clusters<-unique(clustResut$clustering)
									
									###############################plot the profiles of cluster groups
									#										colVec<-rainbow(length(Clusters))
									#										labelVec<-paste("pkID:",rownames(goodShapePeaklist),"mz:",goodShapePeaklist$mz," ,gss:",format(goodShapePeaklist$gss,digits=2)," ,shrp:",format(goodShapePeaklist$shrp,digits=2))
									#										playwith({
									#										plot(x=0,type = "n",xlim=c(minLbound,maxRbound),ylim=c(0,max(maxInt2)))
									#										for(i in 1:length(Clusters))
									#										{
									#											######get current cluster ID
									#											curCluster<-Clusters[i]
									#											#####get masses belong to current cluster
									#											curGroup<-names(clustResut$clustering[clustResut$clustering==curCluster])
									#											
									#											for(j in 1:length(curGroup))
									#											{
									#												
									#												points(profiles[[curGroup[j]]], type = "l",col=colVec[i])
									#											}
									#										}
									#										},data.points=subset(goodShapePeaklist,select=c("pkInd","Intensity")),labels=labelVec)
									###############################									
									isSplit<-F
									
									for(curCluster in Clusters)
									{
										#####get masses belong to current cluster
										curGroup<-names(clustResut$clustering[clustResut$clustering==curCluster])
										#######does not average the single-ion cluster which has zero intra dist
										if(length(curGroup)>1)
										{#######calculate the average intra distance to see if the current group is a valid single group of fragments
											
											if(max(as.dist(r[curGroup,curGroup]))>maxIntraDist)
											{##if current group is not valid, skip to the next K,
												#and the current k value 
												isSplit<-T
												break
											}
											
										}
									}
									if(!isSplit)
										asw[k]<-clustResut$silinfo $ avg.width
								}
								k.best <- which.max(asw)
								if(k.best>1)
								{
									FinalClustResut<-pam(distance,k.best,cluster.only=TRUE,stand=TRUE)
									
									#################################################
									#calculate the average intra dist in each cluster
									###############################################
									isTooBigIntraDist<-NULL
									avgIntraDistVec<-NULL
									Clusters<-unique(FinalClustResut)
									for(curCluster in Clusters)
									{
										#####get masses belong to current cluster
										curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
										
										IntraDist<-mean(as.dist(r[curGroup,curGroup]))
										isTooBigIntraDist<-c(isTooBigIntraDist,ifelse(IntraDist>6,T,F))
										if(length(curGroup)>1)#######does not average the single-ion cluster which has zero intra dist
											avgIntraDistVec<-c(avgIntraDistVec,IntraDist)
										
									}
									#################################################################
									#if after splitting,intra dist is smaller than before splitting
									#then accept the best clustering schema, 
									#and earch for unique mass within each cluster
									###############################################################
									if(mean(avgIntraDistVec)<mean(distance))
										isMultGroups<-TRUE
								}
							}else#if ==2, then directly set flag=T
							{
								isMultGroups<-TRUE
							}
						}
					}
				}
				
				
				
				
				#############################
				#select model peak for each group
				#######################
				modelPkList<-NULL
#					logRTfactor<-log(mean(c(minLbound,maxRbound)),base=2)/log(totalscan,base=2)
				
				if (!isMultGroups && nPks >=1)
				{
					if (nPks==1) {
						modelPkList<-goodShapePeaklist
					} else 
					{
						mdlPkCandidates<-goodShapePeaklist
						#if the group size is larger than threshold then calcalue the score for each peak candidate
						
						#mass score
						mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
						#gassian similarity score
						mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
						#peak height score
						mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))					
						mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
						mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
						mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
						
						#cacluate the over score which takes mass,intensity and gss into account
						scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
						#										scores<-(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3+0.4*mdlPkCandidates$f4)
						#rescale score to (0,1000)
						mdlPkCandidates$score<-scale(scores,F,0.001)
						##select the maximum score
						mdlPkCandidates <- subset(mdlPkCandidates,score>=500)
						if (nrow(mdlPkCandidates)>0)
						{
							modelPkList<- mdlPkCandidates[which.max(mdlPkCandidates$score),]
						}
					}
				}
				
				
				if(isMultGroups)
				{
					if(nPks==2)
					{#if there are only two candidate then simply use them as model peaks
						#modelPkVec<-rownames(goodShapePeaklist)
						modelPkList<-goodShapePeaklist
					}else
					{
						#						scoreCutoff<-500#minimum score for model peak
						nMinFragment<-1
						if(length(Clusters)>0)
						{	#select model peak from each group[!isTooBigIntraDist]
							for(curCluster in Clusters)
							{
								#####get masses belong to current cluster
								curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
								
								###get model peak candidates 
								mdlPkCandidates<-goodShapePeaklist[curGroup,]
								
								################################plot the profiles of cluster groups
#									colVec<-rainbow(length(curGroup))
#									labelVec<-paste("pkID:",rownames(mdlPkCandidates),"mz:",mdlPkCandidates$mz," ,gss:",format(mdlPkCandidates$gss,digits=2)," ,shrp:",format(mdlPkCandidates$shrp,digits=2))
#									playwith({
#												plot(x=0,type = "n",xlim=c(minLbound,maxRbound),ylim=c(0,maxInt2))
#												
#												for(j in 1:length(curGroup))
#												{
#													
#													points(profiles[[curGroup[j]]], type = "l",col=colVec[j])
#												}
#												
#											},data.points=subset(mdlPkCandidates,select=c("pkInd","Intensity")),labels=labelVec)
								###############################
								
								#if the group size is larger than threshold then calcalue the score for each peak candidate
								if(nrow(mdlPkCandidates)>0)
								{
									#mass score
									mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
									#gassian similarity score
									mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
									#peak height score
									mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))					
									#sharpness score
#										mdlPkCandidates$f4<-scale(mdlPkCandidates$shrp,min(mdlPkCandidates$shrp),diff(range(mdlPkCandidates$shrp)))					
									
									###assign zero value when each factor is NA
									mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
									mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
									mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
#										mdlPkCandidates$f4[is.na(mdlPkCandidates$f4)]<-0
									
									#cacluate the over score which takes mass,intensity and gss into account
									scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
#										scores<-(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3+0.4*mdlPkCandidates$f4)
									#rescale score to (0,1000)
									mdlPkCandidates$score<-scale(scores,F,0.001)
									##select the maximum score
									currentMdlPk<-mdlPkCandidates[which.max(mdlPkCandidates$score),]
									
									##########################
									#if the selected model peak has the score larger than threshold 
									#then add it to the final model peak list
									#############################################################
									modelPkList<-rbind(modelPkList,currentMdlPk)
								}
							}
						}
					}
					
				}
				
				
				
				#####################plot the profiles of model peaks
#					colVec<-rainbow(nrow(modelPkList))
#					
#					labelVec<-paste("pkid:",rownames(modelPkList),"mz:",modelPkList$mz," ,gss:",format(modelPkList$gss,digits=2)," ,shrp:",format(modelPkList$shrp,digits=2))
				##					labelVec<-paste("gss=",format(modelPkList$gss,digits=2),"\n shrp=",format(modelPkList$shrp,digits=2))
#					dataPoints<-subset(modelPkList,select=c("pkInd","Intensity"))
#					dataPoints$pkInd<-(dataPoints$pkInd-1)*ScanInterval+delaytime
#					
#					playwith({
#					plot(x=0,type = "n",xlim=(c(minLbound,maxRbound)-1)*ScanInterval+delaytime,ylim=c(0,max(maxInt2)))
#					for(i in 1:nrow(modelPkList))
#					{
#						curPkid<-rownames(modelPkList)[i]
#						
#						points((profiles[[curPkid]]$ET-1)*ScanInterval+delaytime,profiles[[curPkid]]$int, type = "l",col=colVec[i])
#					}
#					},data.points=dataPoints,labels=labelVec)
				####################			
				
			}
			
			
			##################################
			#decompose based on detected model peaks
			##################################
			##################################
			decom = 1
			while (decom==1) 
			{
				modelPkDistance <- NULL
				specList<-NULL
				componentList<-NULL
				mdlPkIDs <- NULL
				#########if no co-eluting compounds (only one or no model peak detected),extract spectrum directly
				if(nrow(modelPkList)<=1||is.null(modelPkList))
				{
					######add model peak to the spectrum
					specList[[1]]<-data.frame(mz=as.integer(),int=as.integer())
					#if no model peak detected, simply extract spectrum from the TIC peak apex time
					if(nrow(modelPkList)==0||is.null(modelPkList))
					{	#get mass from the first peak of the local window
						
						curMz<-as.character(localEICPeakList$mz[1])
						#get the TIC peak apex time
						maxInd<-curTICPk$pkInd
						#get the boundary of TIC peak
						curLbound<-curTICPk$lboundInd
						curRbound<-curTICPk$rboundInd
						#get the pk height of the first peak
						curPhHeight<-localEICPeakList$Intensity[1]
						#get the peak id of the first peak
						curPkID<-rownames(localEICPeakList)[1]
						#get the profile of the first peak 
						curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.x="ET",by.y="ET",all.x=T)
						curProfile1$int[is.na(curProfile1$int)]<-0
						#get the area of the first peak 
						curArea<-sum(curProfile1$int)
						
					}else
					{
						#if there is only one model peak detected,then extract spectrum from the model peak apex time
						mdlPkIDs<-rownames(modelPkList)
						curMz<-modelPkList$mz[1]
						maxInd<-modelPkList$pkInd[1]
						curLbound<-modelPkList$lboundInd[1]
						curRbound<-modelPkList$rboundInd[1]
						curPhHeight<-modelPkList$Intensity[1]
						
						curPkID<-rownames(modelPkList)[1]
						curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
						curProfile1$int[is.na(curProfile1$int)]<-0
						curArea<-sum(curProfile1$int)
						
					}
					#construct the component list for csv output (later used for quantition)
					componentList[[1]]<-data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,pkInd=maxInd,Intensity=curPhHeight,area=curArea,windowID=windowID,compoundID=1)
					#construct the compound name by the RT and model peak mass information for NIST file output(used for Identifiation)
					clustID<-format((maxInd-1)*ScanInterval+delaytime,digit=4)
					names(specList)[1]<-paste(clustID,curMz)
					#add the rest peaks from the local window to the component list and spectral list
					for(i in 1:nrow(localEICPeakList))
					{
						curPkID<-rownames(localEICPeakList)[i]
						
						curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
						curProfile1$int[is.na(curProfile1$int)]<-0
						curArea<-sum(curProfile1$int)
						
						curMz<-as.character(localEICPeakList$mz[i])
						curInt<-subset(curProfile1,ET==maxInd)$int[1]
						
						if(!is.na(curInt)&&curInt>0)
						{
							specList[[1]]<-rbind(specList[[1]],data.frame(mz=as.integer(curMz),int=curInt))######add umass to the spectrum
							
							componentList[[1]]<-rbind(componentList[[1]],data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,pkInd=maxInd,Intensity=curInt,area=curArea,windowID=windowID,compoundID=1))
						}
					}
					decom=0
				}else
				{#if there are more than one model peaks detected then perform decomposition
					
					## to remove those model peaks with small distance
					modelPkDistance <- as.matrix(modelPkDist(modelPkList,profiles))
					## model same matrix as modelPKDistance with index
					psmatrix <- matrix(c(1:length(modelPkDistance)),ncol=nrow(modelPkDistance),byrow=T)
					psmodel <- which(modelPkDistance<=4&modelPkDistance>0)
					rmmodelPk <- NULL
					ModelPkNames <- rownames(modelPkDistance)
					## get index of model peaks with smaller distance
					for (i in psmodel) {
						pos <- matrix.index(psmatrix, i)
						rmmodelPk <- c(rmmodelPk,ModelPkNames[pos[1]],ModelPkNames[pos[2]])
					}
					rmmodelPks <- unique(rmmodelPk)
					
					rmmodelPkList <- modelPkList[rmmodelPks,]
					modelPkList1 <- rmmodelPkList[which.max(rmmodelPkList$score),]
					modelPkList <- modelPkList[!rownames(modelPkList)%in%rmmodelPks,]
					modelPkList <- rbind(modelPkList,modelPkList1)
					
					mdlPkIDs<-rownames(modelPkList)
					
					# if there's only one model peak after testing, go back to the upper step
					if (nrow(modelPkList)==1) decom=1
					
					else {
						lbound<-min(modelPkList$lboundInd)
						rbound<-max(modelPkList$rboundInd)
						#################################
						#get entire profiles within the window for each local peaks
						#and saved it into X matrix;each column is one EIC,each row is one scan
						#############################################
						X<-NULL
						mzVec<-unique(localEICPeakList$mz)
						nMz<-length(mzVec)
						for(i in 1:nMz)
						{
							mz<-mzVec[i]
							mzInd<-which(vectorMz==mz)
							startInd<-(mzInd-1)*totalscan+1
							endInd<-startInd+totalscan-1
							curVecInt <- vecInt[startInd:endInd]
							
							X<-cbind(X,curVecInt[lbound:rbound])
						}
						colnames(X)<-mzVec	
						nScans<-nrow(X)
						
						
						#####################plot entire profile of all masses within the window range defined by umass
#					colVec<-rainbow(nMz)
#					playwith({
#					plot(x=0,type = "n",xlim=(c(lbound,rbound)-1)*ScanInterval+delaytime,ylim=c(0,max(maxInt1)))
#					for(i in 1:nMz)
#					{
#						curMz<-as.character(mzVec[i])
#						points(((lbound:rbound)-1)*ScanInterval+delaytime,X[,curMz], type = "l",col=colVec[i])
#					}
#					})
						############plot entire profiles of model peak mass within the window range defined by umass
#						colVec<-rainbow(nrow(modelPkList))
#						playwith({
#						plot(x=0,type = "n",xlim=c(lbound,rbound),ylim=c(0,maxInt2))
#						for(i in 1:nrow(modelPkList))
#						{
#							curMz<-as.character(modelPkList[i,]$mz)
#							points(lbound:rbound,X[,curMz], type = "l",col=colVec[i])
#						}
#						})
						#####################					
						
						
						######decovolute other peaks by the model peak
						#get model peak profiles and save them into S matrix
						#each column is the EIC for one model peak
						###################################################
						
						S<-NULL
						modelPkList$area<-0
						for(i in 1:nrow(modelPkList))
						{
							#########get the source signal for model peak
							#expand the profile by filling zero values to those uncovered area
							###################################################################
							curProfile1<-merge(data.frame(ET=as.integer(lbound:rbound)),profiles[[rownames(modelPkList[i,])]],by.y="ET",all.x=T)
							curProfile1$int[is.na(curProfile1$int)]<-0
							#get information of model peak and save it into compoent list
							S<-cbind(S,curProfile1$int)
							modelPkList[i,]$area<-sum(curProfile1$int)		
							#####init the spec list	
							componentList[[i]]<-data.frame(mz=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),pkInd=as.integer(),Intensity=as.integer(),area=as.integer(),windowID=as.integer(),compoundID=as.integer())
							######add model peak to the spectrum
							specList[[i]]<-data.frame(mz=as.integer(),int=as.integer())
							curRT<-format((modelPkList$pkInd[i]-1)*ScanInterval+delaytime,digit=4)
							names(specList)[i]<-paste(curRT,modelPkList$mz[i])
						}
						colnames(S)<-mdlPkIDs
						#plot the decomposed signals
						#		playwith({
						#plot(x=0,type = "n",xlim=c(0,nScans),ylim=c(0,max(X)))
						for(i in 1:ncol(X))
						{
							curMz<-colnames(X)[i]
							
							M<-X[,i]#####mixture signal
							srcIDs<-mdlPkIDs
							A<-optim(par=rep(0,ncol(S)),fn=fr, gr = NULL,M,S,lower = 0, method="L-BFGS-B")$par
							A<-as.matrix(A)
							
#						repeat
#						{
#							#########get mixing matrix A by least square approximation
#							A<-solve(t(S[,srcIDs])%*%S[,srcIDs])%*%t(S[,srcIDs])%*%M
							##							A<-pseudoinverse(S)%*%M
#							
#		
#							##remove the model peak when it produces the negative coefficients
#							##and do the decomposition again
#							nInvalid<-length(which(A<0))
#							if(nInvalid>0)
#							{
#								curMdlPkList<-modelPkList[srcIDs,]
							##								rmPkID<-rownames(curMdlPkList[which.min(curMdlPkList$score),])
#								rmPkID<-names(A[A<0,])
#								srcIDs<-srcIDs[srcIDs!=rmPkID]
#							}else
#							{
#								break
#							}
#						}
							#			deconS<-NULL
							##restore the shared ion signals and save them into component list ans speclist
							for(m in 1:length(srcIDs))
							{
								curmdlPkID<-srcIDs[m]
								j<-which(mdlPkIDs==curmdlPkID)
								
								#					deconS<-S[,j]*A[j]######restore orignal signal by multiply maxUmass with mixing factor						
								#					points(deconS,type="l",col=rainbow(length(modelPkVec))[j])
								#########add the restored signals to the respective component list
								#					deconS<-as.matrix(deconS)
								#					colnames(deconS)<-curMz
								if(A[m]>0)
								{
									specList[[j]]<-rbind(specList[[j]],data.frame(mz=as.integer(curMz),int=A[m,]))###add the curent restored mixing factor to spec
									componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=as.integer(modelPkList[curmdlPkID,]$Intensity*A[m,]),area=as.integer(modelPkList[curmdlPkID,]$area*A[m,]),windowID=windowID,compoundID=j))
								}
								
							}
						}
						
#							representModelPk <- NULL
#							# get the pairwise combination index 
#							Num <- nrow(modelPkList)
#							indexPairs<-combinations(Num,2)
#							for(i in 1:nrow(indexPairs))
#							{
#								index <- indexPairs[i,]
#								index1 <- index[1]
#								index2 <- index[2]
#								score <- 0
#								# calculate the time difference between two model peaks
#								modelPKET <- (modelPkList$pkInd[c(index1,index2)]-1)*ScanInterval+delaytime
#								if (abs(modelPKET[1]-modelPKET[2]) <= 0.0042) { # accepting 5 scans difference
#									score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
#									if (score >=970) {
#										# select the one with higher intensity
#										selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index1],mdlPkIDs[index2])
#										representModelPk <- c(representModelPk,selectedModelPkID)
#									}
#									else representModelPk <- c(representModelPk,mdlPkIDs[c(index1,index2)])		
#								}			
#							}
#							# select the unique filtered model peaks 
#							representModelPk <- unique(representModelPk)
#							modelPkList <- modelPkList[rownames(modelPkList)%in%representModelPk,]
#							
#							# decide if re-decompositon is necessary: if the number of modelPk is reduced, then re-decompose
#							if (nrow(modelPkList) == Num) decom = 0
						
						# testing the mass similartiy of spectra after decompositon
						RepeatModelPk <- NULL
# get the pairwise combination index 
						Num <- nrow(modelPkList)
						indexPairs<-combinations(Num,2)
						for(i in 1:nrow(indexPairs))
						{
							index <- indexPairs[i,]
							index1 <- index[1]
							index2 <- index[2]
							# score <- 0
							# calculate the time difference between two model peaks
							modelPKET <- modelPkList$pkInd[c(index1,index2)]
#							if (abs(modelPKET[1]-modelPKET[2]) <= 2) {
#								score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
#								if (score >=970) {
#									# select the one with higher intensity
#									selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index2],mdlPkIDs[index1])
#									RepeatModelPk <- c(RepeatModelPk,selectedModelPkID)
#								}
#								#else representModelPk <- c(representModelPk,mdlPkIDs[c(index1,index2)])		
#							}
							
							score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
							if (score >=970||(abs(modelPKET[1]-modelPKET[2]) <=2&&score >=900))
							{
								selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index2],mdlPkIDs[index1])
								RepeatModelPk <- c(RepeatModelPk,selectedModelPkID)
							}
						}
# select the unique filtered model peaks 
#representModelPk <- unique(representModelPk)
						RepeatModelPk <- unique(RepeatModelPk)
						modelPkList <- modelPkList[!rownames(modelPkList)%in%RepeatModelPk,]
						
# decide if re-decompositon is necessary: if the number of modelPk is reduced, then re-decompose
						if (nrow(modelPkList) == Num) decom = 0
						
					} # else model PK num is >2 after testing distance	
				} # model PK num is >2 at first
			} # while statement
			
			###flag the EIC peaks 
			EICpeaklist[rownames(localEICPeakList),]$flag<-1
			componentResults<-c(componentResults,componentList)
			specResults<-c(specResults,specList)
			
			### comment the unnessary mdlPKIDVec  -- yan
			mdlPkIDVec<-c(mdlPkIDVec,mdlPkIDs)
			
		}	
	}	
	
	componentResults<-do.call(rbind,componentResults)	
	cat(fileName,"finished deconvolution...\n")
	
	##################################
	#collect spectrum result in NIST fomrat
	###################################
	libmatch <- libMatching(lib,specResults,600)
	write.csv(libmatch,file=paste(params$WorkDir,"/output/decon/",fileName,"_libmatch.csv",sep=""))
	
	#################################################
	#collect spectrum result in dataframe for lib matching
	###################################################
	
	evaluateIdentification(refSpec=lib,inputSpec=specResults,RT_Tolerance=70,minSpecSimilarity=650,fileName,withRT=T,params)
	cat(fileName,"finish matching...\n")
	write.csv(mdlPkIDVec,file=paste(params$WorkDir,"/output/decon/",fileName,"mdlpklist.csv",sep=""))
	
#		#################################################
#		#collect component result in dataframe for alignment
#		###################################################
#		##filter out zero-int fragments
	
	componentResults<-subset(componentResults,Intensity>0,select=c("mz","lboundInd","rboundInd","pkInd","Intensity","area","windowID","compoundID"))
	write.csv(componentResults,row.names=FALSE,file=paste(params$WorkDir,"/output/decon/",fileName,"Decon.csv",sep=""))
	cat(fileName,"writing component results...\n")
}


decomposition_old<-function(inFilePath,params,c1,isDistParallel,clustingType)
{
	
	source(paste(params$codeDir,"pipeline.r",sep="/"))
	
	library(gdata)
	library(gtools)
	library(cluster)
	library("ncdf")
	
	options(expressions=1e5)
	
	DataFilelist<-params$DataFiles
	#####load library
	lib<-readMSP2Spec(filename="../StdsLibrary/QcCompound.txt",withRT=T)
	
	delaytime<-params$delaytime
	ScanInterval<-params$ScanInterval
	WorkDir<-params$WorkDir
	
	fileName<-parseFileName(inFilePath)
	
	###########################################
	#get decon window from TIC peak picking results
	###########################################
	cat(fileName,"reading TIC peak data...\n")
	TICPeakFile<-paste(params$WorkDir,"output/peakpicking/",fileName,"_TIC_PeakList.csv",sep="")
	TICApexList <- read.csv(TICPeakFile)
	TICApexList$RT<-(TICApexList$pkInd-1)*ScanInterval+delaytime
	TICpeaklist<-subset(TICApexList,isApex==1)
	
	winIDs<-1:nrow(TICpeaklist)
	
	###########################################################
	#to save time,load the target RT info and filter the TIC peak list 
	#so that the decon is only performed in the selected windows
	#################################################
	
	winIDs<-sort(unique(winIDs))
	
	###############
	#read EIC data
	################
	rawEICFile<-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
	denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
	EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
	
	if(file.exists(EICFile)==FALSE)next
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	#####orignial mz vector from EIC data
	vectorMz<-get.var.ncdf(ncid, varid="mzVec")
	close.ncdf(ncid)
	remove(ncid)	
	
	#########################
	#read peak picking result
	#########################
	EICPeakFile<-paste(WorkDir,"/output/peakpicking/",fileName,"PeakList.csv",sep="")
	cat(fileName,"reading EIC peak data...\n")
	if(file.exists(EICPeakFile)==FALSE)next
	EICpeaklist <- read.csv(EICPeakFile)
	colnames(EICpeaklist)[which(colnames(EICpeaklist)=="curMz")]<-"mz"
	
	#########################
	#start decon
	#########################
	cat(fileName,"start deconvolution...\n")
	
	###################
	#sequential version
	####################
	totalscan<-length(vecInt)/length(vectorMz)
	
	componentResults<-NULL
	specResults<-NULL
	mdlPkIDVec<-NULL
	EICpeaklist$flag=0
	
	for(windowID in winIDs)
	{
		cat("window:",windowID,"\n")
		
		curTICPk<-TICpeaklist[windowID,]
		#decide the exlcuded ion list for current TIC window
		minET=(curTICPk$lboundInd-1)*params$ScanInterval+params$delaytime
		maxET=(curTICPk$rboundInd-1)*params$ScanInterval+params$delaytime
		if(minET>=8)
		{
			CurNonUmassVec<-c(params$NonUmassVec,51:100)
		}else
		{
			CurNonUmassVec<-params$NonUmassVec
		}
		
		####################################################
		#get all the EIC peaks within current decon window
		##################################################
		localTICApexList<-subset(TICApexList,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		#get the potential compound number for TIC peak apex list
		#there should be at least as many or more compounds than the number of apex 
		#minCmp<-nrow(localTICApexList)
		localEICPeakList<-subset(EICpeaklist,flag==0&pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		
		if(nrow(localEICPeakList)>0)
		{
			#maxInt1<-max(localEICPeakList$Intensity)
			
			#and adjust decon window by EIC peak profiles
			#minLbound<-min(localEICPeakList$lboundInd)
			#maxRbound<-max(localEICPeakList$rboundInd)	
			
			allprofiles<-vector("list",nrow(localEICPeakList))
			for(i in 1:nrow(localEICPeakList))
			{
				##############################################
				#get intensity vector the current mz
				# from the long  intensity vector
				################################################
				curEICpk<-localEICPeakList[i,]
				
				startInd<-curEICpk$offset+curEICpk$lboundInd
				endInd<-curEICpk$offset+curEICpk$rboundInd
				
				#######get the profile defined by EIC peak boundary
				ET<-curEICpk$lboundInd:curEICpk$rboundInd
				curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd])
				allprofiles[[i]]<-curProfile
				names(allprofiles)[i]<-rownames(curEICpk)
				
				
			}
			
			########filter out noisy small peaks
			minPkHeight<-2000
			CandidatePeaklist<-subset(localEICPeakList,Intensity>minPkHeight)
			
			###exlude the shared ion by exclude list
			CandidatePeaklist<-subset(CandidatePeaklist,!mz%in%CurNonUmassVec)			
			################################################################
			#remove shared Ion by gaussian similarity and sharpness
			##############################################################		
			
			gaussianCutoff<-1.5
			ShrpCutoff<-5
			
			goodShapePeaklist<-subset(CandidatePeaklist,gss<=gaussianCutoff&shrp>=ShrpCutoff)
			nPks<-nrow(goodShapePeaklist)
			
			modelPkList<-NULL
			isMultGroups<-FALSE
			mdlPkIDs<-NULL
			if(nPks>0)
			{
				#maxInt2<-max(goodShapePeaklist$Intensity)
				
				####get the subset of the profiles for the goodshape peaks
				profiles<-allprofiles[as.character(rownames(goodShapePeaklist))]
				
				####update each profile by duplicating the half of peak that has good shape
				
				for(i in 1:length(profiles))
				{
					curPkId<-names(profiles[i])
					curPk<-goodShapePeaklist[curPkId,]
					pos<-curPk$gssPos
					profiles[[i]]<-copyShape(p=profiles[[i]],from=pos)
					goodShapePeaklist[curPkId,]$lboundInd<-profiles[[i]]$ET[1]
					goodShapePeaklist[curPkId,]$rboundInd<-profiles[[i]]$ET[nrow(profiles[[i]])]
				}
				
				######################################
				#decide the number of groups in the current window
				########################################
				if(nPks>=2)
				{
					#save the EIC profies in global variable for broadcasting to slave nodes for parallel computing
					assign("Global.curProfs", value=profiles, envir = .GlobalEnv)
					#			assign("Global.vecInt", value=vecInt, envir = .GlobalEnv)
					#calculate the mutual distance among the EICs
					
					if(isDistParallel)
					{
						r<-parDistCal4(cl,isUnion=F)##parallel verison
					}else
					{
						r<-DistCal4(isUnion=F)#non-parallel version
					}	
					#convert the distance matrix to triangle matrix
					distance <- as.dist(r)
					maxIntraDist<-9#cutoff to decide if split
					
					if(clustingType=="h")
					{
						#######################################################################################
						#hierarchical clustering and cut into groups by the height
						#####################################################################################
						clustResut<-hclust(distance)
						
						FinalClustResut<-cutree(clustResut,h=maxIntraDist)
						Clusters<-unique(FinalClustResut)
						if(length(Clusters)>1)isMultGroups<-T
						
					}else
					{	#######################################################################################
						#partitional clustering and split into groups by the silehoute
						#####################################################################################
						if(max(distance)<=maxIntraDist)
						{###if the max distance is already small enought, no need to deconvolute
							isMultGroups<-FALSE
						}else
						{
							
							if(nPks>2)
							{
								########################################
								#try different k to decide the cluster number by intra dist
								#######################################
								
								maxClustNumber<-min(20,nPks-1)
								asw<-numeric(maxClustNumber)
								for (k in min(maxClustNumber,max(2,minCmp)):maxClustNumber)
								{
									clustResut <-pam(distance,k)
									Clusters<-unique(clustResut$clustering)
									
									isSplit<-F
									
									for(curCluster in Clusters)
									{
										#####get masses belong to current cluster
										curGroup<-names(clustResut$clustering[clustResut$clustering==curCluster])
										#######does not average the single-ion cluster which has zero intra dist
										if(length(curGroup)>1)
										{#######calculate the average intra distance to see if the current group is a valid single group of fragments
											
											if(max(as.dist(r[curGroup,curGroup]))>maxIntraDist)
											{##if current group is not valid, skip to the next K,
												#and the current k value 
												isSplit<-T
												break
											}
											
										}
									}
									if(!isSplit)
										asw[k]<-clustResut$silinfo $ avg.width
								}
								k.best <- which.max(asw)
								if(k.best>1)
								{
									FinalClustResut<-pam(distance,k.best,cluster.only=TRUE,stand=TRUE)
									
									#################################################
									#calculate the average intra dist in each cluster
									###############################################
									isTooBigIntraDist<-NULL
									avgIntraDistVec<-NULL
									Clusters<-unique(FinalClustResut)
									for(curCluster in Clusters)
									{
										#####get masses belong to current cluster
										curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
										
										IntraDist<-mean(as.dist(r[curGroup,curGroup]))
										isTooBigIntraDist<-c(isTooBigIntraDist,ifelse(IntraDist>6,T,F))
										if(length(curGroup)>1)#######does not average the single-ion cluster which has zero intra dist
											avgIntraDistVec<-c(avgIntraDistVec,IntraDist)
										
									}
									#################################################################
									#if after splitting,intra dist is smaller than before splitting
									#then accept the best clustering schema, 
									#and earch for unique mass within each cluster
									###############################################################
									if(mean(avgIntraDistVec)<mean(distance))
										isMultGroups<-TRUE
								}
							}else#if ==2, then directly set flag=T
							{
								isMultGroups<-TRUE
							}
						}
					}
				}
				
				
				
				
				#############################
				#select model peak for each group
				#######################
				modelPkList<-NULL
#					logRTfactor<-log(mean(c(minLbound,maxRbound)),base=2)/log(totalscan,base=2)
				
				if(isMultGroups)
				{
					if(nPks==2)
					{#if there are only two candidate then simply use them as model peaks
						modelPkVec<-rownames(goodShapePeaklist)
					}else
					{
						#						scoreCutoff<-500#minimum score for model peak
						nMinFragment<-1
						if(length(Clusters)>0)
						{	#select model peak from each group[!isTooBigIntraDist]
							for(curCluster in Clusters)
							{
								#####get masses belong to current cluster
								curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
								
								###get model peak candidates 
								mdlPkCandidates<-goodShapePeaklist[curGroup,]
								
								#if the group size is larger than threshold then calcalue the score for each peak candidate
								if(nrow(mdlPkCandidates)>0)
								{
									#mass score
									mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
									#gassian similarity score
									mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
									#peak height score
									mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))					
									#sharpness score
#										mdlPkCandidates$f4<-scale(mdlPkCandidates$shrp,min(mdlPkCandidates$shrp),diff(range(mdlPkCandidates$shrp)))					
									
									###assign zero value when each factor is NA
									mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
									mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
									mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
#										mdlPkCandidates$f4[is.na(mdlPkCandidates$f4)]<-0
									
									#cacluate the over score which takes mass,intensity and gss into account
									scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
#										scores<-(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3+0.4*mdlPkCandidates$f4)
									#rescale score to (0,1000)
									mdlPkCandidates$score<-scale(scores,F,0.001)
									##select the maximum score
									currentMdlPk<-mdlPkCandidates[which.max(mdlPkCandidates$score),]
									
									##########################
									#if the selected model peak has the score larger than threshold 
									#then add it to the final model peak list
									#############################################################
									modelPkList<-rbind(modelPkList,currentMdlPk)
								}
							}
						}
					}
					
				}
				mdlPkIDs<-rownames(modelPkList)			
				
			}
			##################################
			#decompose based on detected model peaks
			##################################
			##################################
			specList<-NULL
			componentList<-NULL
			#########if no co-eluting compounds (only one or no model peak detected),extract spectrum directly
			if(nrow(modelPkList)<=1||is.null(modelPkList))
			{
				######add model peak to the spectrum
				specList[[1]]<-data.frame(mz=as.integer(),int=as.integer())
				#if no model peak detected, simply extract spectrum from the TIC peak apex time
				if(nrow(modelPkList)==0||is.null(modelPkList))
				{	#get mass from the first peak of the local window
					
					curMz<-as.character(localEICPeakList$mz[1])
					#get the TIC peak apex time
					maxInd<-curTICPk$pkInd
					#get the boundary of TIC peak
					curLbound<-curTICPk$lboundInd
					curRbound<-curTICPk$rboundInd
					#get the pk height of the first peak
					curPhHeight<-localEICPeakList$Intensity[1]
					#get the peak id of the first peak
					curPkID<-rownames(localEICPeakList)[1]
					#get the profile of the first peak 
					curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.x="ET",by.y="ET",all.x=T)
					curProfile1$int[is.na(curProfile1$int)]<-0
					#get the area of the first peak 
					curArea<-sum(curProfile1$int)
					
				}else
				{
					#if there is only one model peak detected,then extract spectrum from the model peak apex time
					curMz<-modelPkList$mz[1]
					maxInd<-modelPkList$pkInd[1]
					curLbound<-modelPkList$lboundInd[1]
					curRbound<-modelPkList$rboundInd[1]
					curPhHeight<-modelPkList$Intensity[1]
					
					curPkID<-rownames(modelPkList)[1]
					curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
					curProfile1$int[is.na(curProfile1$int)]<-0
					curArea<-sum(curProfile1$int)
					
				}
				#construct the component list for csv output (later used for quantition)
				componentList[[1]]<-data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,pkInd=maxInd,Intensity=curPhHeight,area=curArea,windowID=windowID,compoundID=1)
				#construct the compound name by the RT and model peak mass information for NIST file output(used for Identifiation)
				clustID<-format((maxInd-1)*ScanInterval+delaytime,digit=4)
				names(specList)[1]<-paste(clustID,curMz)
				#add the rest peaks from the local window to the component list and spectral list
				for(i in 1:nrow(localEICPeakList))
				{
					curPkID<-rownames(localEICPeakList)[i]
					
					curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
					curProfile1$int[is.na(curProfile1$int)]<-0
					curArea<-sum(curProfile1$int)
					
					curMz<-as.character(localEICPeakList$mz[i])
					curInt<-subset(curProfile1,ET==maxInd)$int[1]
					
					if(!is.na(curInt)&&curInt>0)
					{
						specList[[1]]<-rbind(specList[[1]],data.frame(mz=as.integer(curMz),int=curInt))######add umass to the spectrum
						
						componentList[[1]]<-rbind(componentList[[1]],data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,pkInd=maxInd,Intensity=curInt,area=curArea,windowID=windowID,compoundID=1))
					}
				}
				
				
			}else
			{#if there are more than one model peaks detected then perform decomposition
				lbound<-min(modelPkList$lboundInd)
				rbound<-max(modelPkList$rboundInd)
				#################################
				#get entire profiles within the window for each local peaks
				#and saved it into X matrix;each column is one EIC,each row is one scan
				#############################################
				X<-NULL
				mzVec<-unique(localEICPeakList$mz)
				nMz<-length(mzVec)
				for(i in 1:nMz)
				{
					mz<-mzVec[i]
					mzInd<-which(vectorMz==mz)
					startInd<-(mzInd-1)*totalscan+1
					endInd<-startInd+totalscan-1
					curVecInt <- vecInt[startInd:endInd]
					
					X<-cbind(X,curVecInt[lbound:rbound])
				}
				colnames(X)<-mzVec	
				nScans<-nrow(X)
				######decovolute other peaks by the model peak
				#get model peak profiles and save them into S matrix
				#each column is the EIC for one model peak
				###################################################
				
				S<-NULL
				modelPkList$area<-0
				for(i in 1:nrow(modelPkList))
				{
					#########get the source signal for model peak
					#expand the profile by filling zero values to those uncovered area
					###################################################################
					curProfile1<-merge(data.frame(ET=as.integer(lbound:rbound)),profiles[[rownames(modelPkList[i,])]],by.y="ET",all.x=T)
					curProfile1$int[is.na(curProfile1$int)]<-0
					#get information of model peak and save it into compoent list
					S<-cbind(S,curProfile1$int)
					modelPkList[i,]$area<-sum(curProfile1$int)		
					#####init the spec list	
					componentList[[i]]<-data.frame(mz=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),pkInd=as.integer(),Intensity=as.integer(),area=as.integer(),windowID=as.integer(),compoundID=as.integer())
					######add model peak to the spectrum
					specList[[i]]<-data.frame(mz=as.integer(),int=as.integer())
					curRT<-format((modelPkList$pkInd[i]-1)*ScanInterval+delaytime,digit=4)
					names(specList)[i]<-paste(curRT,modelPkList$mz[i])
				}
				colnames(S)<-mdlPkIDs
				#plot the decomposed signals
				#		playwith({
				plot(x=0,type = "n",xlim=c(0,nScans),ylim=c(0,max(X)))
				for(i in 1:ncol(X))
				{
					curMz<-colnames(X)[i]
					
					M<-X[,i]#####mixture signal
					srcIDs<-mdlPkIDs
					A<-optim(par=rep(0,ncol(S)),fn=fr, gr = NULL,M,S,lower = 0, method="L-BFGS-B")$par
					A<-as.matrix(A)
					
					for(m in 1:length(srcIDs))
					{
						curmdlPkID<-srcIDs[m]
						j<-which(mdlPkIDs==curmdlPkID)
						
						#					deconS<-S[,j]*A[j]######restore orignal signal by multiply maxUmass with mixing factor
						
						#					points(deconS,type="l",col=rainbow(length(modelPkVec))[j])
						#########add the restored signals to the respective component list
						#					deconS<-as.matrix(deconS)
						#					colnames(deconS)<-curMz
						if(A[m]>0)
						{
							specList[[j]]<-rbind(specList[[j]],data.frame(mz=as.integer(curMz),int=A[m,]))###add the curent restored mixing factor to spec
							componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=as.integer(modelPkList[curmdlPkID,]$Intensity*A[m,]),area=as.integer(modelPkList[curmdlPkID,]$area*A[m,]),windowID=windowID,compoundID=j))
						}
						
					}
				}
				
			}
			###flag the EIC peaks 
			EICpeaklist[rownames(localEICPeakList),]$flag<-1
			componentResults<-c(componentResults,componentList)
			specResults<-c(specResults,specList)
			
			mdlPkIDVec<-c(mdlPkIDVec,mdlPkIDs)
			
		}	
		
	}	
	
	componentResults<-do.call(rbind,componentResults)
	
	cat(fileName,"finished deconvolution...\n")
	
	##################################
	#collect spectrum result in NIST fomrat
	###################################
#		specResult<-unlist(MatchResult,recursive=T)
#		write(MatchResult,file=paste(params$WorkDir,"/output/decon/",fileName,"_decon_NIST.txt",sep=""))
#		cat(fileName,"writing spec results...\n")
	
	#################################################
	#collect spectrum result in dataframe for lib matching
	###################################################
	
	evaluateIdentification(refSpec=lib,inputSpec=specResults,RT_Tolerance=20,minSpecSimilarity=650,fileName,withRT=T,params)
	cat(fileName,"finish matching...\n")
#		
	write.csv(mdlPkIDVec,file=paste(params$WorkDir,"/output/decon/",fileName,"mdlpklist.csv",sep=""))
	componentResults<-subset(componentResults,Intensity>0,select=c("mz","lboundInd","rboundInd","pkInd","Intensity","area","windowID","compoundID"))
	write.csv(componentResults,row.names=FALSE,file=paste(params$WorkDir,"/output/decon/",fileName,"Decon.csv",sep=""))
	cat(fileName,"writing component results...\n")
	
}


decomposition_SN<-function(inFilePath,params,cl,isDistParallel,clustingType)
{
	source(paste(params$codeDir,"pipeline.r",sep="/"))	
	library(gdata)
	library(gtools)
	library(cluster)
	library("ncdf")	
	options(expressions=1e5)
	DataFilelist<-params$DataFiles
	delaytime<-params$delaytime
	ScanInterval<-params$ScanInterval
	WorkDir<-params$WorkDir
	
	lib<-readMSP2Spec(filename="../StdsLibrary/KQC.txt",withRT=T)
	fileName<-parseFileName(inFilePath)
	
	###############
	#read TIC data
	################		
	TICfile<-paste(WorkDir,"/output/TIC/denoised_",fileName,"_TIC.cdf",sep="") 
	# The denoised TIC are saved in output/TIC/
	#if there's no denoised TIC file, then read raw cdf data
	TICfile<-ifelse(file.exists(TICfile),TICfile,inFilePath)
	ncid <- open.ncdf(TICfile)
	TIC <- get.var.ncdf(ncid, varid="total_intensity")
	close.ncdf(ncid)
	remove(ncid)
	
	###############
	#read EIC data
	################
	EICFile<-paste(WorkDir,"/output/EIC/",fileName,"EIC.cdf",sep="")
	if(file.exists(EICFile)==FALSE)next
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	#####orignial mz vector from EIC data
	vectorMz<-get.var.ncdf(ncid, varid="mzVec")
	close.ncdf(ncid)
	remove(ncid)	
	
	###########################################
	#get decon window from TIC peak picking results
	###########################################
	cat(fileName,"reading TIC peak data...\n")
	TICPeakFile<-paste(params$WorkDir,"output/peakpicking/",fileName,"_TIC_PeakList.csv",sep="")
	TICApexList <- read.csv(TICPeakFile)
	TICApexList$RT<-(TICApexList$pkInd-1)*ScanInterval+delaytime
	TICpeaklist<-subset(TICApexList,isApex==1)
	
	winIDs<-1:nrow(TICpeaklist)
	
#		RTVector<-read.csv(paste(WorkDir,"output/RTlist.csv",sep=""))$RT
#		TICRTvec<-(TICpeaklist$pkInd-1)*ScanInterval+delaytime
#		winIDs<-NULL		
#		for(i in RTVector)
#		{
#			winIDs<-c(winIDs,which(abs(TICRTvec-i)<=(12/60)))
#		}
#		
#		winIDs<-sort(unique(winIDs))
	
	#########################
	#read EIC peak picking result
	#########################
	cat(fileName,"reading EIC peak data...\n")
	EICPeakFile<-paste(WorkDir,"/output/peakpicking/",fileName,"PeakList.csv",sep="")
	if(file.exists(EICPeakFile)==FALSE)next
	EICpeaklist <- read.csv(EICPeakFile)
	colnames(EICpeaklist)[which(colnames(EICpeaklist)=="curMz")]<-"mz"
	
	#########################
	#start decon
	#########################
	cat(fileName,"start deconvolution...\n")
	
	###################
	#sequential version
	####################	
	componentResults<-NULL
	specResults<-NULL
	mdlPkIDVec<-NULL
	EICpeaklist$flag=0
	totalscan<-length(vecInt)/length(vectorMz)
	
	# to create high-frequency (noise) and low-freq (signal)
	cat(fileName,"calculating signal to noise...\n")
	noiseVec <- NULL
	smVec <- NULL
	for(mzInd in 1:length(vectorMz))
	{
		#mz<-vectorMz[mzInd]	
		startInd<-(mzInd-1)*totalscan+1
		endInd<-startInd+totalscan-1
		curVecInt <- vecInt[startInd:endInd]
		
		rawfeature <- curVecInt
		rawfeature[2:(totalscan-1)] <- (rawfeature[1:(totalscan-2)]+rawfeature[2:(totalscan-1)]+rawfeature[3:totalscan])/3	
		# select abosolute value of noise
		noise <- abs(curVecInt-rawfeature)
		noiseVec <- c(noiseVec,noise)
		smVec <- c(smVec,rawfeature)	
	}
	
	for(windowID in winIDs)
	#for(windowID in c(12:14))
	{
		cat("window:",windowID,"\n")
		curTICPk<-TICpeaklist[windowID,]
		#decide the exlcuded ion list for current TIC window
		minET=(curTICPk$lboundInd-1)*params$ScanInterval+params$delaytime
		maxET=(curTICPk$rboundInd-1)*params$ScanInterval+params$delaytime
		if(minET>=8)
		{
			CurNonUmassVec<-c(params$NonUmassVec,51:100)
			#CurNonUmassVec<-c(params$NonUmassVec,51:72,74:100)
		}else
		{
			CurNonUmassVec<-params$NonUmassVec
		}
		
		####################################################
		#get all the EIC peaks within current decon window
		##################################################
		localTICApexList<-subset(TICApexList,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		#get the potential compound number for TIC peak apex list
		#there should be at least as many or more compounds than the number of apex 
		#minCmp<-nrow(localTICApexList)
		localEICPeakList<-subset(EICpeaklist,flag==0&pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		
		if(nrow(localEICPeakList)>0)
		{
			# save all raw EIC features as allprofiles
#				allprofiles<-vector("list",nrow(localEICPeakList))
#				for(i in 1:nrow(localEICPeakList))
#				{
#					##############################################
#					#get intensity vector the current mz
#					# from the long  intensity vector
#					################################################
#					curEICpk<-localEICPeakList[i,]			
#					startInd<-curEICpk$offset+curEICpk$lboundInd
#					endInd<-curEICpk$offset+curEICpk$rboundInd
#					
#					#######get the profile defined by EIC peak boundary
#					ET<-curEICpk$lboundInd:curEICpk$rboundInd
#					curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd])
#					allprofiles[[i]]<-curProfile
#					names(allprofiles)[i]<-rownames(curEICpk)
#										
#				}
#				
			allprofiles<-vector("list",nrow(localEICPeakList))	
			sigNoise <- NULL
			
			for(i in 1:nrow(localEICPeakList))
			{
				##############################################
				#get intensity vector the current mz
				# from the long  intensity vector
				################################################
				curEICpk<-localEICPeakList[i,]
				
				startInd<-curEICpk$offset+curEICpk$lboundInd
				endInd<-curEICpk$offset+curEICpk$rboundInd
				
				noise <- noiseVec[startInd:endInd]
				signal <- smVec[startInd:endInd]
				sumSignal <- sum(signal)
				sumNoise <- sum(noise)
				sumNoise <- ifelse(sumNoise==0,1,sumNoise)
				sn <- sumSignal/sumNoise
				sigNoise <- c(sigNoise,sn)
				#######get the profile defined by EIC peak boundary
				ET<-curEICpk$lboundInd:curEICpk$rboundInd
				curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd])
				allprofiles[[i]]<-curProfile
				names(allprofiles)[i]<-rownames(curEICpk)
			}
			
			# add sigal to noise value to the local peak list
			localEICPeakList$sn <-sigNoise 
			
			################################################################
			#remove shared Ion by gaussian similarity and sharpness
			#And filter out noisy small peaks under 1000, or select top 50
			##############################################################		
#				gaussianCutoff<-2
#				ShrpCutoff<-5
#				#minPkHeight<- 1000 
#				minPkHeight <- mean(TIC[localTICApexList$pkInd])*0.025
#				
#				CandidatePeaklist<-subset(localEICPeakList,gss<=gaussianCutoff&shrp>=ShrpCutoff)
#				goodShapePeaklist<-subset(CandidatePeaklist,Intensity>minPkHeight)
#				goodShapePeaklist <- goodShapePeaklist[order(goodShapePeaklist[,"Intensity"],decreasing=T),]
#				# If top ten good peak mass has occurred in non-unique mass, keep them
#				ImportantLowMass <- goodShapePeaklist[order(goodShapePeaklist[,"Intensity"],decreasing=T),][c(1:10),]$mz
#				CurNonUmassVec<- CurNonUmassVec[!CurNonUmassVec%in%ImportantLowMass]
#				# remove other low and non-unique mass
#				goodShapePeaklist<-subset(goodShapePeaklist,!mz%in%CurNonUmassVec)
#				
#				if (nrow(goodShapePeaklist) >= 50) 
#				{
#					goodShapePeaklist <- goodShapePeaklist[order(goodShapePeaklist[,"Intensity"],decreasing=T),][c(1:50),]
#				}
#				nPks<-nrow(goodShapePeaklist)
			
			
			# create a scoring system combining GSS, Shrp and S/N by multiplication
			Gsslist <- 90-localEICPeakList$gss
			shrplist <- localEICPeakList$shrp
			snlist <- localEICPeakList$sn	
			localEICPeakList$Gssscore<-scale(Gsslist,min(Gsslist),diff(range(Gsslist)))
			localEICPeakList$shrpscore<-scale(shrplist,min(shrplist),diff(range(shrplist)))			
			localEICPeakList$snscore<-scale(snlist,min(snlist),diff(range(snlist)))
			###assign zero value when each factor is NA
			localEICPeakList$Gssscore[is.na(localEICPeakList$Gssscore)]<-0
			localEICPeakList$shrpscore[is.na(localEICPeakList$shrpscore)]<-0
			localEICPeakList$snscore[is.na(localEICPeakList$snscore)]<-0
			# combine three parameters together
			scores<-localEICPeakList$shrpscore*localEICPeakList$Gssscore*localEICPeakList$snscore
			localEICPeakList$scores<-scores
			
			# rough cutoff
			gaussianCutoff<-3
			ShrpCutoff<-5
			CandidatePeaklist<-subset(localEICPeakList,gss<=gaussianCutoff&shrp>=ShrpCutoff&Intensity>=500)				
			#CandidatePeaklist<-subset(localEICPeakList,gss<=gaussianCutoff&shrp>=ShrpCutoff)
			#goodShapePeaklist<-subset(CandidatePeaklist,scores>=0.1)
			
			if (nrow(CandidatePeaklist)>0) 
			{
#					b <- min(CandidatePeaklist$scores)
#					#cut <- b+0.25*diff(range(CandidatePeaklist$scores))		
#					cut <- b+0.2*diff(range(CandidatePeaklist$scores))
#					goodShapePeaklist<-subset(CandidatePeaklist,scores>=cut)
				
				a <- max(CandidatePeaklist$scores)
				b <- min(CandidatePeaklist$scores)
				if (a <= 0.02) goodShapePeaklist<- data.frame() # blank data frame
				
				if (a > 0.02) {
					# b >=0.2, keep all candidate peaks
					if (b < 0.2) {
						cut_temp <- b+0.25*diff(range(CandidatePeaklist$scores))
						cut <- ifelse(cut_temp < 0.05,0.05,cut_temp)
						CandidatePeaklist<-subset(CandidatePeaklist,scores>=cut)
					}
					
					if (nrow(CandidatePeaklist)>=50)
					{
						goodShapePeaklist <- CandidatePeaklist[order(CandidatePeaklist[,"scores"],decreasing=T),][c(1:50),]
					}
					else 
						goodShapePeaklist<- CandidatePeaklist
				}
				
				nPks<-nrow(goodShapePeaklist)
				modelPkList<-NULL
				isMultGroups<-FALSE
				mdlPkIDs<-NULL
				
				if(nPks>0)
				{
					####get the subset of the profiles for the goodshape peaks
					profiles<-allprofiles[as.character(rownames(goodShapePeaklist))]
					####update each profile by duplicating the half of peak that has good shape
					for(i in 1:length(profiles))
					{
						curPkId<-names(profiles[i])
						curPk<-goodShapePeaklist[curPkId,]
						pos<-curPk$gssPos
						profiles[[i]]<-copyShape(p=profiles[[i]],from=pos)
						goodShapePeaklist[curPkId,]$lboundInd<-profiles[[i]]$ET[1]
						goodShapePeaklist[curPkId,]$rboundInd<-profiles[[i]]$ET[nrow(profiles[[i]])]
					}
					
					######################################
					#decide the number of groups in the current window
					########################################
					if(nPks>=2)
					{
						#save the EIC profies in global variable for broadcasting to slave nodes for parallel computing
						assign("Global.curProfs", value=profiles, envir = .GlobalEnv)
						#			assign("Global.vecInt", value=vecInt, envir = .GlobalEnv)
						#calculate the mutual distance among the EICs
						
						if(isDistParallel)
						{
							r<-parDistCal4(cl,isUnion=F)##parallel verison
						}else
						{
							r<-DistCal4(isUnion=F)#non-parallel version
						}	
						#convert the distance matrix to triangle matrix
						distance <- as.dist(r)
						maxIntraDist<-15 #cutoff to decide if split
						
						if(clustingType=="h")
						{
							#######################################################################################
							#hierarchical clustering and cut into groups by the height
							#####################################################################################
							clustResut<-hclust(distance)
							FinalClustResut<-cutree(clustResut,h=maxIntraDist)
							Clusters<-unique(FinalClustResut)
							if(length(Clusters)>1)isMultGroups<-T
							
						}else
						{	#######################################################################################
							#partitional clustering and split into groups by the silehoute
							#####################################################################################
							if(max(distance)<=maxIntraDist)
							{###if the max distance is already small enought, no need to deconvolute
								isMultGroups<-FALSE
							}else
							{
								
								if(nPks>2)
								{
									########################################
									#try different k to decide the cluster number by intra dist
									#######################################
									
									maxClustNumber<-min(20,nPks-1)
									asw<-numeric(maxClustNumber)
									for (k in min(maxClustNumber,max(2,minCmp)):maxClustNumber)
									{
										clustResut <-pam(distance,k)
										Clusters<-unique(clustResut$clustering)								
										isSplit<-F
										
										for(curCluster in Clusters)
										{
											#####get masses belong to current cluster
											curGroup<-names(clustResut$clustering[clustResut$clustering==curCluster])
											#######does not average the single-ion cluster which has zero intra dist
											if(length(curGroup)>1)
											{#######calculate the average intra distance to see if the current group is a valid single group of fragments
												
												if(max(as.dist(r[curGroup,curGroup]))>maxIntraDist)
												{##if current group is not valid, skip to the next K,
													#and the current k value 
													isSplit<-T
													break
												}
												
											}
										}
										if(!isSplit)
											asw[k]<-clustResut$silinfo $ avg.width
									}
									k.best <- which.max(asw)
									if(k.best>1)
									{
										FinalClustResut<-pam(distance,k.best,cluster.only=TRUE,stand=TRUE)
										
										#################################################
										#calculate the average intra dist in each cluster
										###############################################
										isTooBigIntraDist<-NULL
										avgIntraDistVec<-NULL
										Clusters<-unique(FinalClustResut)
										for(curCluster in Clusters)
										{
											#####get masses belong to current cluster
											curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
											
											IntraDist<-mean(as.dist(r[curGroup,curGroup]))
											isTooBigIntraDist<-c(isTooBigIntraDist,ifelse(IntraDist>6,T,F))
											if(length(curGroup)>1)#######does not average the single-ion cluster which has zero intra dist
												avgIntraDistVec<-c(avgIntraDistVec,IntraDist)
											
										}
										#################################################################
										#if after splitting,intra dist is smaller than before splitting
										#then accept the best clustering schema, 
										#and earch for unique mass within each cluster
										###############################################################
										if(mean(avgIntraDistVec)<mean(distance))
											isMultGroups<-TRUE
									}
								}else#if ==2, then directly set flag=T
								{
									isMultGroups<-TRUE
								}
							}
						}
					}
					
					
					
					
					#############################
					#select model peak for each group
					#######################
					modelPkList<-NULL
					
					# if there's one compound, trying to select one peak instead of noise to represent this compound
					if (!isMultGroups && nPks >=1)
					{
						mdlPkCandidates<-goodShapePeaklist
						#if the group size is larger than threshold then calcalue the score for each peak candidate
						
						#mass score
						mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
						#gassian similarity score
						mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
						#peak height score
						mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))					
						mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
						mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
						mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
						
						#cacluate the over score which takes mass,intensity and gss into account
						scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
						#										scores<-(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3+0.4*mdlPkCandidates$f4)
						#rescale score to (0,1000)
						mdlPkCandidates$score<-scale(scores,F,0.001)
						##select the maximum score
						mdlPkCandidates <- subset(mdlPkCandidates,score>=500)
						if (nrow(mdlPkCandidates)>0)
						{
							modelPkList<- mdlPkCandidates[which.max(mdlPkCandidates$score),]
						}
						
					}
					
					if(isMultGroups)
					{
						if(nPks==2)
						{#if there are only two candidate then simply use them as model peaks
							#modelPkVec<-rownames(goodShapePeaklist)
							modelPkList<-goodShapePeaklist
						}else
						{
							#						scoreCutoff<-500#minimum score for model peak
							nMinFragment<-1
							if(length(Clusters)>0)
							{	#select model peak from each group[!isTooBigIntraDist]
								for(curCluster in Clusters)
								{
									#####get masses belong to current cluster
									curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
									
									###get model peak candidates 
									mdlPkCandidates<-goodShapePeaklist[curGroup,]
									
									#if the group size is larger than threshold then calcalue the score for each peak candidate
									if(nrow(mdlPkCandidates)>0)
									{
										#mass score
										mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
										#gassian similarity score
										mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
										#peak height score
										mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))					
										#sharpness score
#										mdlPkCandidates$f4<-scale(mdlPkCandidates$shrp,min(mdlPkCandidates$shrp),diff(range(mdlPkCandidates$shrp)))					
										
										###assign zero value when each factor is NA
										mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
										mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
										mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
#										mdlPkCandidates$f4[is.na(mdlPkCandidates$f4)]<-0
										
										#cacluate the over score which takes mass,intensity and gss into account
										scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
#										scores<-(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3+0.4*mdlPkCandidates$f4)
										#rescale score to (0,1000)
										mdlPkCandidates$score<-scale(scores,F,0.001)
										##select the maximum score
										
										##########################
										#if the selected model peak has the score larger than threshold 
										#then add it to the final model peak list
										#############################################################
										mdlPkCandidates <- subset(mdlPkCandidates,score>=500)
										if (nrow(mdlPkCandidates)>0)
										{
											currentMdlPk<-mdlPkCandidates[which.max(mdlPkCandidates$score),]
											modelPkList<-rbind(modelPkList,currentMdlPk)	
										}
										
									}
								}
							}
						}
						
					}								
				}
				
				
				##################################
				#decompose based on detected model peaks
				##################################
				
				decom = 1
				while (decom==1) 
				{
					modelPkDistance <- NULL
					specList<-NULL
					componentList<-NULL
					
					#########if no co-eluting compounds (only one or no model peak detected),extract spectrum directly
					if(nrow(modelPkList)<=1||is.null(modelPkList))
					{
						######add model peak to the spectrum
						specList[[1]]<-data.frame(mz=as.integer(),int=as.integer())
						#if no model peak detected, simply extract spectrum from the TIC peak apex time
						if(nrow(modelPkList)==0||is.null(modelPkList))
						{	#get mass from the first peak of the local window
							
							curMz<-as.character(localEICPeakList$mz[1])
							#get the TIC peak apex time
							maxInd<-curTICPk$pkInd
							#get the boundary of TIC peak
							curLbound<-curTICPk$lboundInd
							curRbound<-curTICPk$rboundInd
							#get the pk height of the first peak
							curPhHeight<-localEICPeakList$Intensity[1]
							#get the peak id of the first peak
							curPkID<-rownames(localEICPeakList)[1]
							#get the profile of the first peak 
							curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.x="ET",by.y="ET",all.x=T)
							curProfile1$int[is.na(curProfile1$int)]<-0
							#get the area of the first peak 
							curArea<-sum(curProfile1$int)
							
						}else
						{
							#if there is only one model peak detected,then extract spectrum from the model peak apex time
							curMz<-modelPkList$mz[1]
							maxInd<-modelPkList$pkInd[1]
							curLbound<-modelPkList$lboundInd[1]
							curRbound<-modelPkList$rboundInd[1]
							curPhHeight<-modelPkList$Intensity[1]
							
							curPkID<-rownames(modelPkList)[1]
							curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
							curProfile1$int[is.na(curProfile1$int)]<-0
							curArea<-sum(curProfile1$int)
							
						}
						#construct the component list for csv output (later used for quantition)
						componentList[[1]]<-data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,pkInd=maxInd,Intensity=curPhHeight,area=curArea,windowID=windowID,compoundID=1)
						#construct the compound name by the RT and model peak mass information for NIST file output(used for Identifiation)
						clustID<-format((maxInd-1)*ScanInterval+delaytime,digit=4)
						names(specList)[1]<-paste(clustID,curMz)
						#add the rest peaks from the local window to the component list and spectral list
						for(i in 1:nrow(localEICPeakList))
						{
							curPkID<-rownames(localEICPeakList)[i]
							
							curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
							curProfile1$int[is.na(curProfile1$int)]<-0
							curArea<-sum(curProfile1$int)
							
							curMz<-as.character(localEICPeakList$mz[i])
							curInt<-subset(curProfile1,ET==maxInd)$int[1]
							
							if(!is.na(curInt)&&curInt>0)
							{
								specList[[1]]<-rbind(specList[[1]],data.frame(mz=as.integer(curMz),int=curInt))######add umass to the spectrum
								
								componentList[[1]]<-rbind(componentList[[1]],data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,pkInd=maxInd,Intensity=curInt,area=curArea,windowID=windowID,compoundID=1))
							}
						}
						componentList[[1]]<- componentList[[1]][!duplicated(componentList[[1]]),]
						decom=0
					}else
					{#if there are more than one model peaks detected then perform decomposition
						
						## to remove those model peaks with small distance
						modelPkDistance <- as.matrix(modelPkDist(modelPkList,profiles))
						## model same matrix as modelPKDistance with index
						psmatrix <- matrix(c(1:length(modelPkDistance)),ncol=nrow(modelPkDistance),byrow=T)
						psmodel <- which(modelPkDistance<=4&modelPkDistance>0)
						rmmodelPk <- NULL
						ModelPkNames <- rownames(modelPkDistance)
						## get index of model peaks with smaller distance
						for (i in psmodel) {
							pos <- matrix.index(psmatrix, i)
							rmmodelPk <- c(rmmodelPk,ModelPkNames[pos[1]],ModelPkNames[pos[2]])
						}
						rmmodelPks <- unique(rmmodelPk)
						
						rmmodelPkList <- modelPkList[rmmodelPks,]
						modelPkList1 <- rmmodelPkList[which.max(rmmodelPkList$score),]
						modelPkList <- modelPkList[!rownames(modelPkList)%in%rmmodelPks,]
						modelPkList <- rbind(modelPkList,modelPkList1)
						
						mdlPkIDs<-rownames(modelPkList)
						
						# if there's only one model peak after testing, go back to the upper step
						if (nrow(modelPkList)==1) decom=1
						
						else {
							lbound<-min(modelPkList$lboundInd)
							rbound<-max(modelPkList$rboundInd)
							#################################
							#get entire profiles within the window for each local peaks
							#and saved it into X matrix;each column is one EIC,each row is one scan
							#############################################
							X<-NULL
							mzVec<-unique(localEICPeakList$mz)
							nMz<-length(mzVec)
							for(i in 1:nMz)
							{
								mz<-mzVec[i]
								mzInd<-which(vectorMz==mz)
								startInd<-(mzInd-1)*totalscan+1
								endInd<-startInd+totalscan-1
								curVecInt <- vecInt[startInd:endInd]
								
								X<-cbind(X,curVecInt[lbound:rbound])
							}
							colnames(X)<-mzVec	
							nScans<-nrow(X)
							
							######decovolute other peaks by the model peak
							#get model peak profiles and save them into S matrix
							#each column is the EIC for one model peak
							###################################################
							
							S<-NULL
							modelPkList$area<-0
							for(i in 1:nrow(modelPkList))
							{
								#########get the source signal for model peak
								#expand the profile by filling zero values to those uncovered area
								###################################################################
								curProfile1<-merge(data.frame(ET=as.integer(lbound:rbound)),profiles[[rownames(modelPkList[i,])]],by.y="ET",all.x=T)
								curProfile1$int[is.na(curProfile1$int)]<-0
								#get information of model peak and save it into compoent list
								S<-cbind(S,curProfile1$int)
								modelPkList[i,]$area<-sum(curProfile1$int)		
								#####init the spec list	
								componentList[[i]]<-data.frame(mz=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),pkInd=as.integer(),Intensity=as.integer(),area=as.integer(),windowID=as.integer(),compoundID=as.integer())
								######add model peak to the spectrum
								specList[[i]]<-data.frame(mz=as.integer(),int=as.integer())
								curRT<-format((modelPkList$pkInd[i]-1)*ScanInterval+delaytime,digit=4)
								names(specList)[i]<-paste(curRT,modelPkList$mz[i])
							}
							colnames(S)<-mdlPkIDs
							
							for(i in 1:ncol(X))
							{
								curMz<-colnames(X)[i]
								
								M<-X[,i]#####mixture signal
								srcIDs<-mdlPkIDs
								A<-optim(par=rep(0,ncol(S)),fn=fr, gr = NULL,M,S,lower = 0, method="L-BFGS-B")$par
								A<-as.matrix(A)
								
								
								##restore the shared ion signals and save them into component list ans speclist
								for(m in 1:length(srcIDs))
								{
									curmdlPkID<-srcIDs[m]
									j<-which(mdlPkIDs==curmdlPkID)
									
									if(A[m]>0)
									{
										specList[[j]]<-rbind(specList[[j]],data.frame(mz=as.integer(curMz),int=A[m,]))###add the curent restored mixing factor to spec
										componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=as.integer(modelPkList[curmdlPkID,]$Intensity*A[m,]),area=as.integer(modelPkList[curmdlPkID,]$area*A[m,]),windowID=windowID,compoundID=j))
									}
									componentList[[j]]<- componentList[[j]][!duplicated(componentList[[j]]),]
								}
							}
							
							# testing the mass similartiy of spectra after decompositon
							RepeatModelPk <- NULL
							# get the pairwise combination index 
							Num <- nrow(modelPkList)
							indexPairs<-combinations(Num,2)
							for(i in 1:nrow(indexPairs))
							{
								index <- indexPairs[i,]
								index1 <- index[1]
								index2 <- index[2]
								score <- 0
								# calculate the time difference between two model peaks
								modelPKET <- modelPkList$pkInd[c(index1,index2)]
								if (abs(modelPKET[1]-modelPKET[2]) <= 2) {
									score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
									#print(paste(score,":",i))
									if (!is.na(score) && score >=970) {
										# select the one with higher intensity
										selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index2],mdlPkIDs[index1])
										RepeatModelPk <- c(RepeatModelPk,selectedModelPkID)
									}
									#else representModelPk <- c(representModelPk,mdlPkIDs[c(index1,index2)])		
								}
							}
							# select the unique filtered model peaks 
							#representModelPk <- unique(representModelPk)
							RepeatModelPk <- unique(RepeatModelPk)
							modelPkList <- modelPkList[!rownames(modelPkList)%in%RepeatModelPk,]
							
							# decide if re-decompositon is necessary: if the number of modelPk is reduced, then re-decompose
							if (nrow(modelPkList) == Num) decom = 0
							
						} # else model PK num is >2 after testing distance	
					} # model PK num is >2 at first
				} # while statement
				
				###flag the EIC peaks 
				EICpeaklist[rownames(localEICPeakList),]$flag<-1
				componentResults<-c(componentResults,componentList)
				specResults<-c(specResults,specList)
				
				### comment the unnessary mdlPKIDVec  -- yan
				#mdlPkIDVec<-c(mdlPkIDVec,mdlPkIDs)
			}# candidate
			
		}	# localpeak
	} # winID for loop 		
	
	componentResults<-do.call(rbind,componentResults)	
	cat(fileName,"finished deconvolution...\n")
	
	##################################
	#collect spectrum result in NIST fomrat
	###################################
#		specResult<-unlist(MatchResult,recursive=T)
#		write(MatchResult,file=paste(params$WorkDir,"/output/decon/",fileName,"_decon_NIST.txt",sep=""))
#		cat(fileName,"writing spec results...\n")
	
	#################################################
	#collect spectrum result in dataframe for lib matching
	###################################################
	
	evaluateIdentification(refSpec=lib,inputSpec=specResults,RT_Tolerance=70,minSpecSimilarity=600,fileName,withRT=T,params)
	cat(fileName,"finish matching...\n")
	#write.csv(mdlPkIDVec,file=paste(params$WorkDir,"/output/decon/",fileName,"mdlpklist.csv",sep=""))
	
#		#################################################
#		#collect component result in dataframe for alignment
#		###################################################
#		##filter out zero-int fragments
	
	componentResults<-subset(componentResults,Intensity>0,select=c("mz","lboundInd","rboundInd","pkInd","Intensity","area","windowID","compoundID"))
	write.csv(componentResults,row.names=FALSE,file=paste(params$WorkDir,"/output/decon/",fileName,"Decon.csv",sep=""))
	cat(fileName,"writing component results...\n")
	
}


decomposition_raw_plot <-function(inFilePath,params,cl,isDistParallel,clustingType)
{
	source(paste(params$codeDir,"pipeline.r",sep="/"))	
	library(gdata)
	library(gtools)
	library(cluster)
	library("ncdf")	
	options(expressions=1e5)
	DataFilelist<-params$DataFiles
	delaytime<-params$delaytime
	ScanInterval<-params$ScanInterval
	WorkDir<-params$WorkDir
	
	#####load library
	##  save all libraries into StdsLibrary -yan
	lib<-NULL
	if(params$JobName=="SerumQC")
	{
		lib<-readMSP2Spec(filename="../StdsLibrary/QcCompound.txt",withRT=T)
	}else
	{
		lib<-readMSP2Spec(filename="../StdsLibrary/KQC.txt",withRT=T)
	}
	
	fileName<-parseFileName(inFilePath)
	###########################################
	#get decon window from TIC peak picking results
	###########################################
	cat(fileName,"reading TIC peak data...\n")
	TICPeakFile<-paste(params$WorkDir,"output/peakpicking/",fileName,"_TIC_PeakList.csv",sep="")
	TICApexList <- read.csv(TICPeakFile)
	TICApexList$RT<-(TICApexList$pkInd-1)*ScanInterval+delaytime
	TICpeaklist<-subset(TICApexList,isApex==1)
	
	###########################	
#		##read TIC 	
	TICfile<-paste(WorkDir,"output/TIC/denoised_",fileName,"_TIC.cdf",sep="") 
	# The denoised TIC are saved under output directory -yan
	#if TIC not existed in seperate file, then read it from raw cdf
	TICfile<-ifelse(file.exists(TICfile),TICfile,inFilePath)
	ncid <- open.ncdf(TICfile)
	TIC <- get.var.ncdf(ncid, varid="total_intensity")
	close.ncdf(ncid)
	remove(ncid)
	
	###########################		
	#TICPkmean <- mean(TIC[TICpeaklist$pkInd])
	winIDs<-1:nrow(TICpeaklist)
	
	###########################################################
	#to save time,load the target RT info and filter the TIC peak list 
	#so that the decon is only performed in the selected windows
	#################################################
	
	winIDs<-sort(unique(winIDs))
	
	###############
	#read EIC data
	################
	rawEICFile<-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
	denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
	EICfile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
	
	if(file.exists(EICFile)==FALSE)next
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	#####orignial mz vector from EIC data
	vectorMz<-get.var.ncdf(ncid, varid="mzVec")
	close.ncdf(ncid)
	remove(ncid)	
	#########################
	#read peak picking result
	#########################
	EICPeakFile<-paste(WorkDir,"output/peakpicking/",fileName,"PeakList.csv",sep="")
	cat(fileName,"reading EIC peak data...\n")
	if(file.exists(EICPeakFile)==FALSE)next
	EICpeaklist <- read.csv(EICPeakFile)
	colnames(EICpeaklist)[which(colnames(EICpeaklist)=="curMz")]<-"mz"
	
	#########################
	#start decon
	#########################
	cat(fileName,"start deconvolution...\n")
	
	###################
	#sequential version
	####################
	totalscan<-length(vecInt)/length(vectorMz)
	
	componentResults<-NULL
	specResults<-NULL
	mdlPkIDVec<-NULL
	EICpeaklist$flag=0
	
	for(windowID in winIDs)
	{
		cat("window:",windowID,"\n")
		curTICPk<-TICpeaklist[windowID,]
		#decide the exlcuded ion list for current TIC window
		minET=(curTICPk$lboundInd-1)*params$ScanInterval+params$delaytime
		maxET=(curTICPk$rboundInd-1)*params$ScanInterval+params$delaytime
		if(minET>=8)
		{
			CurNonUmassVec<-c(params$NonUmassVec,51:100)
			#CurNonUmassVec<-c(params$NonUmassVec,51:72,74:100)
		}else
		{
			CurNonUmassVec<-params$NonUmassVec
		}
		
		####################################################
		#get all the EIC peaks within current decon window
		##################################################
		localTICApexList<-subset(TICApexList,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		#get the potential compound number for TIC peak apex list
		#there should be at least as many or more compounds than the number of apex 
		minCmp<-nrow(localTICApexList)
		localEICPeakList<-subset(EICpeaklist,flag==0&pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		
		if(nrow(localEICPeakList)>0)
		{
			maxInt1<-max(localEICPeakList$Intensity)
			
			##and adjust decon window by EIC peak profiles
			## It seems to be unuseful, temporily comment --- yan
			
			#minLbound<-min(localEICPeakList$lboundInd)
			#maxRbound<-max(localEICPeakList$rboundInd)	
			
			allprofiles<-vector("list",nrow(localEICPeakList))
			for(i in 1:nrow(localEICPeakList))
			{
				##############################################
				#get intensity vector the current mz
				# from the long  intensity vector
				################################################
				curEICpk<-localEICPeakList[i,]
				
				startInd<-curEICpk$offset+curEICpk$lboundInd
				endInd<-curEICpk$offset+curEICpk$rboundInd
				
				#######get the profile defined by EIC peak boundary
				ET<-curEICpk$lboundInd:curEICpk$rboundInd
				curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd])
				allprofiles[[i]]<-curProfile
				names(allprofiles)[i]<-rownames(curEICpk)
				
			}
			
			########filter out noisy small peaks
			
			#minPkHeight<-TICPkmean*0.05
			minPkHeight <- mean(TIC[localTICApexList$pkInd])*0.025
			CandidatePeaklist<-subset(localEICPeakList,Intensity>minPkHeight)
			
			###exlude the shared ion by exclude list
			#CandidatePeaklist<-subset(CandidatePeaklist,!mz%in%CurNonUmassVec)
			
			#################plot the profiles of all peaks
#				colVec<-rainbow(nrow(CandidatePeaklist))
#				names(colVec)<-CandidatePeaklist$mz
#				labelVec<-paste("mz:",CandidatePeaklist$mz," ,gss:",format(CandidatePeaklist$gss,digits=2)," ,shrp:",format(CandidatePeaklist$shrp,digits=2))
#				dataPoints<-subset(CandidatePeaklist,select=c("pkInd","Intensity"))
#				dataPoints$pkInd<-(dataPoints$pkInd-1)*ScanInterval+delaytime
#				playwith({
#							plot(x=0,type = "n",xlim=(c(min(CandidatePeaklist$lboundInd),max(CandidatePeaklist$rboundInd))-1)*ScanInterval+delaytime,ylim=c(0,maxInt1))
#							for(i in 1:nrow(CandidatePeaklist))
#							{
#								curMz<-as.character(CandidatePeaklist[i,]$mz)
#								curID<-rownames(CandidatePeaklist[i,])
#								points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, type = "l",col=colVec[curMz])
#								
#							}
#						},data.points=dataPoints,labels=labelVec)
			##################	
			
			
			################################################################
			#remove shared Ion by gaussian similarity and sharpness
			##############################################################
			if (minPkHeight <=1500) 
			{
				gaussianCutoff<-3
			}else 
			{
				gaussianCutoff<-1.5
			}
			ShrpCutoff<-5
			goodShapePeaklist<-subset(CandidatePeaklist,gss<=gaussianCutoff&shrp>=ShrpCutoff)
			
			nPks<-nrow(goodShapePeaklist)
			
			modelPkList<-NULL
			isMultGroups<-FALSE
			mdlPkIDs<-NULL
			if(nPks>0)
			{
				maxInt2<-max(goodShapePeaklist$Intensity)
				
				####get the subset of the profiles for the goodshape peaks
				profiles<-allprofiles[as.character(rownames(goodShapePeaklist))]
				
				####update each profile by duplicating the half of peak that has good shape
				
				for(i in 1:length(profiles))
				{
					curPkId<-names(profiles[i])
					curPk<-goodShapePeaklist[curPkId,]
					pos<-curPk$gssPos
					profiles[[i]]<-copyShape(p=profiles[[i]],from=pos)
					goodShapePeaklist[curPkId,]$lboundInd<-profiles[[i]]$ET[1]
					goodShapePeaklist[curPkId,]$rboundInd<-profiles[[i]]$ET[nrow(profiles[[i]])]
				}
				
				#####################plot the profiles of selected peaks with good gaussian shape
#					colVec<-rainbow(nrow(goodShapePeaklist))
#					labelVec<-paste("pkID:",rownames(goodShapePeaklist),"mz:",goodShapePeaklist$mz," ,gss:",format(goodShapePeaklist$gss,digits=2)," ,shrp:",format(goodShapePeaklist$shrp,digits=2))
#					dataPoints<-subset(goodShapePeaklist,select=c("pkInd","Intensity"))
#					dataPoints$pkInd<-(dataPoints$pkInd-1)*ScanInterval+delaytime
#					playwith({
#					plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,max(maxInt2)))
#					for(i in 1:nrow(goodShapePeaklist))
#					{
#						curID<-rownames(goodShapePeaklist[i,])
#						points((profiles[[curID]]$ET-1)*ScanInterval+delaytime,profiles[[curID]]$int, type = "l",col=colVec[i])
#						
#					}
#					},data.points=dataPoints,labels=labelVec)
				#####################
				
				
				######################################
				#decide the number of groups in the current window
				########################################
				if(nPks>=2)
				{
					#save the EIC profies in global variable for broadcasting to slave nodes for parallel computing
					assign("Global.curProfs", value=profiles, envir = .GlobalEnv)
					#			assign("Global.vecInt", value=vecInt, envir = .GlobalEnv)
					#calculate the mutual distance among the EICs
					
					if(isDistParallel)
					{
						r<-parDistCal4(cl,isUnion=F)##parallel verison
					}else
					{
						r<-DistCal4(isUnion=F)#non-parallel version
					}	
					#convert the distance matrix to triangle matrix
					distance <- as.dist(r)
					maxIntraDist<-15 #cutoff to decide if split
					
					if(clustingType=="h")
					{
						#######################################################################################
						#hierarchical clustering and cut into groups by the height
						#####################################################################################
						clustResut<-hclust(distance)
#									plot(clustResut)
						FinalClustResut<-cutree(clustResut,h=maxIntraDist)
						Clusters<-unique(FinalClustResut)
						if(length(Clusters)>1)isMultGroups<-T
						
						
						##########################plot the profiles of cluster groups
#							colVec<-rainbow(length(Clusters))
#							labelVec<-paste("pkID:",rownames(goodShapePeaklist),"mz:",goodShapePeaklist$mz," ,gss:",format(goodShapePeaklist$gss,digits=2)," ,shrp:",format(goodShapePeaklist$shrp,digits=2))
#							playwith({
#							plot(x=0,type = "n",xlim=(c(minLbound,maxRbound)-1)*ScanInterval+delaytime,ylim=c(0,max(maxInt2)))
#							for(i in 1:length(Clusters))
#							{
#								######get current cluster ID
#								curCluster<-Clusters[i]
#								#####get masses belong to current cluster
#								curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
#								
#								for(j in 1:length(curGroup))
#								{
#									
#									points((profiles[[curGroup[j]]]$ET-1)*ScanInterval+delaytime,profiles[[curGroup[j]]]$int, type = "l",col=colVec[i])
#								}
#							}
#							},data.points=subset(goodShapePeaklist,select=c("pkInd","Intensity")),labels=labelVec)
#							
						###########################							
						
						
					}else
					{	#######################################################################################
						#partitional clustering and split into groups by the silehoute
						#####################################################################################
						if(max(distance)<=maxIntraDist)
						{###if the max distance is already small enought, no need to deconvolute
							isMultGroups<-FALSE
						}else
						{
							
							if(nPks>2)
							{
								########################################
								#try different k to decide the cluster number by intra dist
								#######################################
								
								maxClustNumber<-min(20,nPks-1)
								asw<-numeric(maxClustNumber)
								for (k in min(maxClustNumber,max(2,minCmp)):maxClustNumber)
								{
									clustResut <-pam(distance,k)
									Clusters<-unique(clustResut$clustering)
									
									###############################plot the profiles of cluster groups
									#										colVec<-rainbow(length(Clusters))
									#										labelVec<-paste("pkID:",rownames(goodShapePeaklist),"mz:",goodShapePeaklist$mz," ,gss:",format(goodShapePeaklist$gss,digits=2)," ,shrp:",format(goodShapePeaklist$shrp,digits=2))
									#										playwith({
									#										plot(x=0,type = "n",xlim=c(minLbound,maxRbound),ylim=c(0,max(maxInt2)))
									#										for(i in 1:length(Clusters))
									#										{
									#											######get current cluster ID
									#											curCluster<-Clusters[i]
									#											#####get masses belong to current cluster
									#											curGroup<-names(clustResut$clustering[clustResut$clustering==curCluster])
									#											
									#											for(j in 1:length(curGroup))
									#											{
									#												
									#												points(profiles[[curGroup[j]]], type = "l",col=colVec[i])
									#											}
									#										}
									#										},data.points=subset(goodShapePeaklist,select=c("pkInd","Intensity")),labels=labelVec)
									###############################									
									isSplit<-F
									
									for(curCluster in Clusters)
									{
										#####get masses belong to current cluster
										curGroup<-names(clustResut$clustering[clustResut$clustering==curCluster])
										#######does not average the single-ion cluster which has zero intra dist
										if(length(curGroup)>1)
										{#######calculate the average intra distance to see if the current group is a valid single group of fragments
											
											if(max(as.dist(r[curGroup,curGroup]))>maxIntraDist)
											{##if current group is not valid, skip to the next K,
												#and the current k value 
												isSplit<-T
												break
											}
											
										}
									}
									if(!isSplit)
										asw[k]<-clustResut$silinfo $ avg.width
								}
								k.best <- which.max(asw)
								if(k.best>1)
								{
									FinalClustResut<-pam(distance,k.best,cluster.only=TRUE,stand=TRUE)
									
									#################################################
									#calculate the average intra dist in each cluster
									###############################################
									isTooBigIntraDist<-NULL
									avgIntraDistVec<-NULL
									Clusters<-unique(FinalClustResut)
									for(curCluster in Clusters)
									{
										#####get masses belong to current cluster
										curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
										
										IntraDist<-mean(as.dist(r[curGroup,curGroup]))
										isTooBigIntraDist<-c(isTooBigIntraDist,ifelse(IntraDist>6,T,F))
										if(length(curGroup)>1)#######does not average the single-ion cluster which has zero intra dist
											avgIntraDistVec<-c(avgIntraDistVec,IntraDist)
										
									}
									#################################################################
									#if after splitting,intra dist is smaller than before splitting
									#then accept the best clustering schema, 
									#and earch for unique mass within each cluster
									###############################################################
									if(mean(avgIntraDistVec)<mean(distance))
										isMultGroups<-TRUE
								}
							}else#if ==2, then directly set flag=T
							{
								isMultGroups<-TRUE
							}
						}
					}
				}
				
				
				
				
				#############################
				#select model peak for each group
				#######################
				modelPkList<-NULL
#					logRTfactor<-log(mean(c(minLbound,maxRbound)),base=2)/log(totalscan,base=2)
				
				if (!isMultGroups && nPks >=1)
				{
					mdlPkCandidates<-goodShapePeaklist
					#if the group size is larger than threshold then calcalue the score for each peak candidate
					
					#mass score
					mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
					#gassian similarity score
					mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
					#peak height score
					mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))					
					mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
					mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
					mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
					
					#cacluate the over score which takes mass,intensity and gss into account
					scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
					#										scores<-(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3+0.4*mdlPkCandidates$f4)
					#rescale score to (0,1000)
					mdlPkCandidates$score<-scale(scores,F,0.001)
					##select the maximum score
					mdlPkCandidates <- subset(mdlPkCandidates,score>=500)
					if (nrow(mdlPkCandidates)>0)
					{
						modelPkList<- mdlPkCandidates[which.max(mdlPkCandidates$score),]
					}
					
				}
				
				
				if(isMultGroups)
				{
					if(nPks==2)
					{#if there are only two candidate then simply use them as model peaks
						#modelPkVec<-rownames(goodShapePeaklist)
						modelPkList<-goodShapePeaklist
					}else
					{
						#						scoreCutoff<-500#minimum score for model peak
						nMinFragment<-1
						if(length(Clusters)>0)
						{	#select model peak from each group[!isTooBigIntraDist]
							for(curCluster in Clusters)
							{
								#####get masses belong to current cluster
								curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
								
								###get model peak candidates 
								mdlPkCandidates<-goodShapePeaklist[curGroup,]
								
								################################plot the profiles of cluster groups
#									colVec<-rainbow(length(curGroup))
#									labelVec<-paste("pkID:",rownames(mdlPkCandidates),"mz:",mdlPkCandidates$mz," ,gss:",format(mdlPkCandidates$gss,digits=2)," ,shrp:",format(mdlPkCandidates$shrp,digits=2))
#									playwith({
#												plot(x=0,type = "n",xlim=c(minLbound,maxRbound),ylim=c(0,maxInt2))
#												
#												for(j in 1:length(curGroup))
#												{
#													
#													points(profiles[[curGroup[j]]], type = "l",col=colVec[j])
#												}
#												
#											},data.points=subset(mdlPkCandidates,select=c("pkInd","Intensity")),labels=labelVec)
								###############################
								
								#if the group size is larger than threshold then calcalue the score for each peak candidate
								if(nrow(mdlPkCandidates)>0)
								{
									#mass score
									mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
									#gassian similarity score
									mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
									#peak height score
									mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))					
									#sharpness score
#										mdlPkCandidates$f4<-scale(mdlPkCandidates$shrp,min(mdlPkCandidates$shrp),diff(range(mdlPkCandidates$shrp)))					
									
									###assign zero value when each factor is NA
									mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
									mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
									mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
#										mdlPkCandidates$f4[is.na(mdlPkCandidates$f4)]<-0
									
									#cacluate the over score which takes mass,intensity and gss into account
									scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
#										scores<-(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3+0.4*mdlPkCandidates$f4)
									#rescale score to (0,1000)
									mdlPkCandidates$score<-scale(scores,F,0.001)
									##select the maximum score
									currentMdlPk<-mdlPkCandidates[which.max(mdlPkCandidates$score),]
									
									##########################
									#if the selected model peak has the score larger than threshold 
									#then add it to the final model peak list
									#############################################################
									modelPkList<-rbind(modelPkList,currentMdlPk)
								}
							}
						}
					}
					
				}
				
				
				
				#####################plot the profiles of model peaks
#					colVec<-rainbow(nrow(modelPkList))
#					
#					labelVec<-paste("pkid:",rownames(modelPkList),"mz:",modelPkList$mz," ,gss:",format(modelPkList$gss,digits=2)," ,shrp:",format(modelPkList$shrp,digits=2))
				##					labelVec<-paste("gss=",format(modelPkList$gss,digits=2),"\n shrp=",format(modelPkList$shrp,digits=2))
#					dataPoints<-subset(modelPkList,select=c("pkInd","Intensity"))
#					dataPoints$pkInd<-(dataPoints$pkInd-1)*ScanInterval+delaytime
#					
#					playwith({
#					plot(x=0,type = "n",xlim=(c(minLbound,maxRbound)-1)*ScanInterval+delaytime,ylim=c(0,max(maxInt2)))
#					for(i in 1:nrow(modelPkList))
#					{
#						curPkid<-rownames(modelPkList)[i]
#						
#						points((profiles[[curPkid]]$ET-1)*ScanInterval+delaytime,profiles[[curPkid]]$int, type = "l",col=colVec[i])
#					}
#					},data.points=dataPoints,labels=labelVec)
				####################			
				
			}
			
			
			##################################
			#decompose based on detected model peaks
			##################################
			##################################
			decom = 1
			while (decom==1) 
			{
				modelPkDistance <- NULL
				specList<-NULL
				componentList<-NULL
				#########if no co-eluting compounds (only one or no model peak detected),extract spectrum directly
				if(nrow(modelPkList)<=1||is.null(modelPkList))
				{
					######add model peak to the spectrum
					specList[[1]]<-data.frame(mz=as.integer(),int=as.integer())
					#if no model peak detected, simply extract spectrum from the TIC peak apex time
					if(nrow(modelPkList)==0||is.null(modelPkList))
					{	#get mass from the first peak of the local window
						
						curMz<-as.character(localEICPeakList$mz[1])
						#get the TIC peak apex time
						maxInd<-curTICPk$pkInd
						#get the boundary of TIC peak
						curLbound<-curTICPk$lboundInd
						curRbound<-curTICPk$rboundInd
						#get the pk height of the first peak
						curPhHeight<-localEICPeakList$Intensity[1]
						#get the peak id of the first peak
						curPkID<-rownames(localEICPeakList)[1]
						#get the profile of the first peak 
						curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.x="ET",by.y="ET",all.x=T)
						curProfile1$int[is.na(curProfile1$int)]<-0
						#get the area of the first peak 
						curArea<-sum(curProfile1$int)
						
					}else
					{
						#if there is only one model peak detected,then extract spectrum from the model peak apex time
						curMz<-modelPkList$mz[1]
						maxInd<-modelPkList$pkInd[1]
						curLbound<-modelPkList$lboundInd[1]
						curRbound<-modelPkList$rboundInd[1]
						curPhHeight<-modelPkList$Intensity[1]
						
						curPkID<-rownames(modelPkList)[1]
						curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
						curProfile1$int[is.na(curProfile1$int)]<-0
						curArea<-sum(curProfile1$int)
						
					}
					#construct the component list for csv output (later used for quantition)
					componentList[[1]]<-data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,pkInd=maxInd,Intensity=curPhHeight,area=curArea,windowID=windowID,compoundID=1)
					#construct the compound name by the RT and model peak mass information for NIST file output(used for Identifiation)
					clustID<-format((maxInd-1)*ScanInterval+delaytime,digit=4)
					names(specList)[1]<-paste(clustID,curMz)
					#add the rest peaks from the local window to the component list and spectral list
					for(i in 1:nrow(localEICPeakList))
					{
						curPkID<-rownames(localEICPeakList)[i]
						
						curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
						curProfile1$int[is.na(curProfile1$int)]<-0
						curArea<-sum(curProfile1$int)
						
						curMz<-as.character(localEICPeakList$mz[i])
						curInt<-subset(curProfile1,ET==maxInd)$int[1]
						
						if(!is.na(curInt)&&curInt>0)
						{
							specList[[1]]<-rbind(specList[[1]],data.frame(mz=as.integer(curMz),int=curInt))######add umass to the spectrum
							
							componentList[[1]]<-rbind(componentList[[1]],data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,pkInd=maxInd,Intensity=curInt,area=curArea,windowID=windowID,compoundID=1))
						}
					}
					decom=0
				}else
				{#if there are more than one model peaks detected then perform decomposition
					
					## to remove those model peaks with small distance
					modelPkDistance <- as.matrix(modelPkDist(modelPkList,profiles))
					## model same matrix as modelPKDistance with index
					psmatrix <- matrix(c(1:length(modelPkDistance)),ncol=nrow(modelPkDistance),byrow=T)
					psmodel <- which(modelPkDistance<=4&modelPkDistance>0)
					rmmodelPk <- NULL
					ModelPkNames <- rownames(modelPkDistance)
					## get index of model peaks with smaller distance
					for (i in psmodel) {
						pos <- matrix.index(psmatrix, i)
						rmmodelPk <- c(rmmodelPk,ModelPkNames[pos[1]],ModelPkNames[pos[2]])
					}
					rmmodelPks <- unique(rmmodelPk)
					
					rmmodelPkList <- modelPkList[rmmodelPks,]
					modelPkList1 <- rmmodelPkList[which.max(rmmodelPkList$score),]
					modelPkList <- modelPkList[!rownames(modelPkList)%in%rmmodelPks,]
					modelPkList <- rbind(modelPkList,modelPkList1)
					
					mdlPkIDs<-rownames(modelPkList)
					
					# if there's only one model peak after testing, go back to the upper step
					if (nrow(modelPkList)==1) decom=1
					
					else {
						lbound<-min(modelPkList$lboundInd)
						rbound<-max(modelPkList$rboundInd)
						#################################
						#get entire profiles within the window for each local peaks
						#and saved it into X matrix;each column is one EIC,each row is one scan
						#############################################
						X<-NULL
						mzVec<-unique(localEICPeakList$mz)
						nMz<-length(mzVec)
						for(i in 1:nMz)
						{
							mz<-mzVec[i]
							mzInd<-which(vectorMz==mz)
							startInd<-(mzInd-1)*totalscan+1
							endInd<-startInd+totalscan-1
							curVecInt <- vecInt[startInd:endInd]
							
							X<-cbind(X,curVecInt[lbound:rbound])
						}
						colnames(X)<-mzVec	
						nScans<-nrow(X)
						
						
						#####################plot entire profile of all masses within the window range defined by umass
#					colVec<-rainbow(nMz)
#					playwith({
#					plot(x=0,type = "n",xlim=(c(lbound,rbound)-1)*ScanInterval+delaytime,ylim=c(0,max(maxInt1)))
#					for(i in 1:nMz)
#					{
#						curMz<-as.character(mzVec[i])
#						points(((lbound:rbound)-1)*ScanInterval+delaytime,X[,curMz], type = "l",col=colVec[i])
#					}
#					})
						############plot entire profiles of model peak mass within the window range defined by umass
#						colVec<-rainbow(nrow(modelPkList))
#						playwith({
#						plot(x=0,type = "n",xlim=c(lbound,rbound),ylim=c(0,maxInt2))
#						for(i in 1:nrow(modelPkList))
#						{
#							curMz<-as.character(modelPkList[i,]$mz)
#							points(lbound:rbound,X[,curMz], type = "l",col=colVec[i])
#						}
#						})
						#####################					
						
						
						######decovolute other peaks by the model peak
						#get model peak profiles and save them into S matrix
						#each column is the EIC for one model peak
						###################################################
						
						S<-NULL
						modelPkList$area<-0
						for(i in 1:nrow(modelPkList))
						{
							#########get the source signal for model peak
							#expand the profile by filling zero values to those uncovered area
							###################################################################
							curProfile1<-merge(data.frame(ET=as.integer(lbound:rbound)),profiles[[rownames(modelPkList[i,])]],by.y="ET",all.x=T)
							curProfile1$int[is.na(curProfile1$int)]<-0
							#get information of model peak and save it into compoent list
							S<-cbind(S,curProfile1$int)
							modelPkList[i,]$area<-sum(curProfile1$int)		
							#####init the spec list	
							componentList[[i]]<-data.frame(mz=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),pkInd=as.integer(),Intensity=as.integer(),area=as.integer(),windowID=as.integer(),compoundID=as.integer())
							######add model peak to the spectrum
							specList[[i]]<-data.frame(mz=as.integer(),int=as.integer())
							curRT<-format((modelPkList$pkInd[i]-1)*ScanInterval+delaytime,digit=4)
							names(specList)[i]<-paste(curRT,modelPkList$mz[i])
						}
						colnames(S)<-mdlPkIDs
						#plot the decomposed signals
						#		playwith({
						#plot(x=0,type = "n",xlim=c(0,nScans),ylim=c(0,max(X)))
						for(i in 1:ncol(X))
						{
							curMz<-colnames(X)[i]
							
							M<-X[,i]#####mixture signal
							srcIDs<-mdlPkIDs
							A<-optim(par=rep(0,ncol(S)),fn=fr, gr = NULL,M,S,lower = 0, method="L-BFGS-B")$par
							A<-as.matrix(A)
							
#						repeat
#						{
#							#########get mixing matrix A by least square approximation
#							A<-solve(t(S[,srcIDs])%*%S[,srcIDs])%*%t(S[,srcIDs])%*%M
							##							A<-pseudoinverse(S)%*%M
#							
#		
#							##remove the model peak when it produces the negative coefficients
#							##and do the decomposition again
#							nInvalid<-length(which(A<0))
#							if(nInvalid>0)
#							{
#								curMdlPkList<-modelPkList[srcIDs,]
							##								rmPkID<-rownames(curMdlPkList[which.min(curMdlPkList$score),])
#								rmPkID<-names(A[A<0,])
#								srcIDs<-srcIDs[srcIDs!=rmPkID]
#							}else
#							{
#								break
#							}
#						}
							#			deconS<-NULL
							##restore the shared ion signals and save them into component list ans speclist
							for(m in 1:length(srcIDs))
							{
								curmdlPkID<-srcIDs[m]
								j<-which(mdlPkIDs==curmdlPkID)
								
								#					deconS<-S[,j]*A[j]######restore orignal signal by multiply maxUmass with mixing factor						
								#					points(deconS,type="l",col=rainbow(length(modelPkVec))[j])
								#########add the restored signals to the respective component list
								#					deconS<-as.matrix(deconS)
								#					colnames(deconS)<-curMz
								if(A[m]>0)
								{
									specList[[j]]<-rbind(specList[[j]],data.frame(mz=as.integer(curMz),int=A[m,]))###add the curent restored mixing factor to spec
									componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=as.integer(modelPkList[curmdlPkID,]$Intensity*A[m,]),area=as.integer(modelPkList[curmdlPkID,]$area*A[m,]),windowID=windowID,compoundID=j))
								}
								
							}
						}
						
#							representModelPk <- NULL
#							# get the pairwise combination index 
#							Num <- nrow(modelPkList)
#							indexPairs<-combinations(Num,2)
#							for(i in 1:nrow(indexPairs))
#							{
#								index <- indexPairs[i,]
#								index1 <- index[1]
#								index2 <- index[2]
#								score <- 0
#								# calculate the time difference between two model peaks
#								modelPKET <- (modelPkList$pkInd[c(index1,index2)]-1)*ScanInterval+delaytime
#								if (abs(modelPKET[1]-modelPKET[2]) <= 0.0042) { # accepting 5 scans difference
#									score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
#									if (score >=970) {
#										# select the one with higher intensity
#										selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index1],mdlPkIDs[index2])
#										representModelPk <- c(representModelPk,selectedModelPkID)
#									}
#									else representModelPk <- c(representModelPk,mdlPkIDs[c(index1,index2)])		
#								}			
#							}
#							# select the unique filtered model peaks 
#							representModelPk <- unique(representModelPk)
#							modelPkList <- modelPkList[rownames(modelPkList)%in%representModelPk,]
#							
#							# decide if re-decompositon is necessary: if the number of modelPk is reduced, then re-decompose
#							if (nrow(modelPkList) == Num) decom = 0
						
						# testing the mass similartiy of spectra after decompositon
						RepeatModelPk <- NULL
# get the pairwise combination index 
						Num <- nrow(modelPkList)
						indexPairs<-combinations(Num,2)
						for(i in 1:nrow(indexPairs))
						{
							index <- indexPairs[i,]
							index1 <- index[1]
							index2 <- index[2]
							score <- 0
							# calculate the time difference between two model peaks
							modelPKET <- modelPkList$pkInd[c(index1,index2)]
							if (abs(modelPKET[1]-modelPKET[2]) <= 2) {
								score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
								if (score >=970) {
									# select the one with higher intensity
									selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index2],mdlPkIDs[index1])
									RepeatModelPk <- c(RepeatModelPk,selectedModelPkID)
								}
								#else representModelPk <- c(representModelPk,mdlPkIDs[c(index1,index2)])		
							}
						}
# select the unique filtered model peaks 
#representModelPk <- unique(representModelPk)
						RepeatModelPk <- unique(RepeatModelPk)
						modelPkList <- modelPkList[!rownames(modelPkList)%in%RepeatModelPk,]
						
# decide if re-decompositon is necessary: if the number of modelPk is reduced, then re-decompose
						if (nrow(modelPkList) == Num) decom = 0
						
					} # else model PK num is >2 after testing distance	
				} # model PK num is >2 at first
			} # while statement
			
			###flag the EIC peaks 
			EICpeaklist[rownames(localEICPeakList),]$flag<-1
			componentResults<-c(componentResults,componentList)
			specResults<-c(specResults,specList)
			
			### comment the unnessary mdlPKIDVec  -- yan
			mdlPkIDVec<-c(mdlPkIDVec,mdlPkIDs)
			
			
			###########################################
			#output decon result in component format (adap1.0 csv format)
			############################################
			#			componentList
			#################################
			#output spec in dataframe format
			################################
			#	specList
			
#	#			#	##	######search  lib
#				res<-libMatching(lib,inputSpec=specList,minSpecSimilarity=650)
#				res
			
			#################
			#output NIST format
			##################
			#	strOutput<-NULL
			#	
			#	for(j in 1:length(specList))
			#	{
			#		
			#		componentName<-paste(names(specList[j]),sep="")
			#		curRT<-as.numeric(unlist(strsplit(componentName," "))[1])
			#		strOutput<-paste(strOutput,getMSPBlock(subset(specList[[j]],int>0),curRT*60,componentName,dbind=j),sep="\n\n")
			#	}
			#	strOutput
			
		}	
	}	
	
	componentResults<-do.call(rbind,componentResults)	
	cat(fileName,"finished deconvolution...\n")
	
	##################################
	#collect spectrum result in NIST fomrat
	###################################
#		specResult<-unlist(MatchResult,recursive=T)
#		write(MatchResult,file=paste(params$WorkDir,"/output/decon/",fileName,"_decon_NIST.txt",sep=""))
#		cat(fileName,"writing spec results...\n")
	libmatch <- libMatching(lib,specResults,650)
	write.csv(libmatch,file=paste(params$WorkDir,"/output/decon/",fileName,"_libmatch.csv",sep=""))
	
	#################################################
	#collect spectrum result in dataframe for lib matching
	###################################################
	
	evaluateIdentification(refSpec=lib,inputSpec=specResults,RT_Tolerance=70,minSpecSimilarity=650,fileName,withRT=T,params)
	cat(fileName,"finish matching...\n")
	write.csv(mdlPkIDVec,file=paste(params$WorkDir,"/output/decon/",fileName,"mdlpklist.csv",sep=""))
	
#		#################################################
#		#collect component result in dataframe for alignment
#		###################################################
#		##filter out zero-int fragments
	
	componentResults<-subset(componentResults,Intensity>0,select=c("mz","lboundInd","rboundInd","pkInd","Intensity","area","windowID","compoundID"))
	write.csv(componentResults,row.names=FALSE,file=paste(params$WorkDir,"/output/decon/",fileName,"Decon.csv",sep=""))
	cat(fileName,"writing component results...\n")
	
}

getPeaks_old<-function(vecInt,params,mz=NULL,mzVec=NULL)
{
	###############################
	#if EIC, then get the int vector
	#from the mz position
	########################
	if(!is.null(mz))
	{
		totalscan<-length(vecInt)/length(mzVec)
		mzInd<-which(mzVec==mz)
		startInd<-(mzInd-1)*totalscan+1
		endInd<-startInd+totalscan-1
		curVecInt <- vecInt[startInd:endInd]
		BHR<-0.3##boundary/height raio
		edgeHightDiffRatio<-0.2
		offset1<-startInd-1
		
	}else
	{######if TIC ,int vector is directly from parameter
		curVecInt<-vecInt
		totalscan<-length(curVecInt)
		BHR<-0.4#0.2
		edgeHightDiffRatio<-0.3#0.1
		
	}
	
	
	WorkDir<-params$WorkDir
	EicPkRatio<-as.integer(params$EicPkRatio)
	nPoints<-length(curVecInt)
	delaytime<-params$delaytime
	ScanInterval<-params$ScanInterval
	
	###############
	#apex detection
	################
	
	peakSpan<-min(9,ifelse(nPoints%%2==1,nPoints,nPoints-1))#peak detection
	isMax <- peaks(x=curVecInt, span=peakSpan)
	isPeak <- isMax #& isSoN & isZeroThrsh	
	peakInd<-which(isPeak==TRUE)#				
	
	
#	##########plot peaks
#	playwith({
#				plot(((1:totalscan)-1)*ScanInterval+delaytime,curVecInt,type="l")
#				points((peakInd-1)*ScanInterval+delaytime,curVecInt[peakInd],col="red")
#					
#				})
	
	
	intCutoff<-0
	######valley detection
	valleySpan<-5
	isMin <- peaks(-curVecInt, valleySpan)
	
	#####include those low-int points as valleys too
	#since sometime the consecutive low points may
	#not be detected by local minima if they have the same values
	#################################
	isZeroThrsh<-curVecInt<=intCutoff
	isValley <- isMin|isZeroThrsh
	valleyInd<-which(isValley==TRUE)	
	
	curPeakList<-data.frame(Intensity=curVecInt[peakInd],pkInd=peakInd)
	#######order by peak index
	curPeakList<-curPeakList[order(peakInd),]
	
	
	########################
	#boudnary detection
	######################
	
	curPeakList$lboundInd<-0
	curPeakList$rboundInd<-0
	curPeakList$isApex<-0
	curPeakList$isShared<-0
	maxWindowlength<-350
	maxlocalPkInd<-0######init the local peak apex position
	lbound<-0
	rbound<-0#########init the rbound 
	preRbound<-0
	
	for(i in 1:nrow(curPeakList))
	{
		curPeakInd<-curPeakList$pkInd[i]
		curPeakHight<-curPeakList$Intensity[i]
		
		####check if current peak apex is within the previous right boundary
		####if so, skip it since it has already been merged
		if(curPeakInd>preRbound)
		{
			
			#######lbound detection
			######search lbound between previous rbound and current peak apex
			leftValleyIndVec<-valleyInd[valleyInd<curPeakInd&valleyInd>=preRbound]
			if(length(leftValleyIndVec)==0)
			{
				#######no left bound
				curPeakList[i,]$isApex=-1
			}else
			{
				########lboudnary has to satisfy a cerntain ratio of peak width/height
				leftValleyIndVec1<-leftValleyIndVec[curVecInt[leftValleyIndVec]/curPeakHight<=BHR]
				
				lbound<-ifelse(length(leftValleyIndVec1)>0,min(leftValleyIndVec1),min(leftValleyIndVec))
				lbound<-max(lbound,0)
				lboundHeight<-curVecInt[lbound]
				
				curPeakList[i,]$lboundInd<-lbound
				##if reach the last peak apex,
				#or simply no value at the right hand side of peak
				#directly set rbound as the last scan
				if(i==nrow(curPeakList)||length(valleyInd[valleyInd>curPeakInd])==0)
				{
					rbound=nPoints
					
					
					localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
					maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
					curPeakList[i,]$rboundInd<-rbound
					if(length(localPkInd)==1)
					{
						curPeakList[i,]$isApex=1
					}else############merge the peak apexes within current boundarys.
					{
						curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-1
						pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
						curPeakList[pkID,]$isApex=1
						curPeakList[pkID,]$lboundInd=lbound
						curPeakList[pkID,]$rboundInd=rbound
						curPeakList[pkID,]$isShared<-1
					}
					preRbound<-rbound
					
				}else
				{
					#########search rbound within the range of from current peakind to maxWindowlength
					rightValleyIndVec<-valleyInd[valleyInd>curPeakInd&valleyInd<lbound+maxWindowlength]
					if(length(rightValleyIndVec)>0)
					{
						######filter by peak height
						rightValleyIndVec1<-rightValleyIndVec[curVecInt[rightValleyIndVec]/curPeakHight<=BHR]
						if(length(rightValleyIndVec1)>0)
						{
							#########filter by lbound height
							rightValleyIndVec2<-rightValleyIndVec1[abs(curVecInt[rightValleyIndVec1]-lboundHeight)/curPeakHight<=edgeHightDiffRatio]
							if(length(rightValleyIndVec2)>0)
							{
								rbound<-min(rightValleyIndVec2)
								rbound<-min(rbound,nPoints)
								curPeakList[i,]$rboundInd<-rbound
								preRbound<-rbound
								
								localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
								maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
								if(length(localPkInd)==1)
								{
									
									curPeakList[i,]$isApex=1
								}else############merge the peak apexes within current boundarys.
								{
									curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-1
									
									pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
									curPeakList[pkID,]$isApex=1
									curPeakList[pkID,]$lboundInd=lbound
									curPeakList[pkID,]$rboundInd=rbound
									curPeakList[pkID,]$isShared<-1
								}
								
							}else
							{	########################################################################
								#if edge difference is significant,try to merge
								#########################################################################
								rbound<-rightValleyIndVec1[1]
								########if lbound higher than rbound, then merge to previous peak
								if(maxlocalPkInd>0&&(lboundHeight>curVecInt[rbound]))
								{
									previousLbound<-curPeakList[curPeakList$pkInd==maxlocalPkInd,]$lboundInd		
									if((rbound-previousLbound)<=maxWindowlength)
									{
										localPkInd<-peakInd[peakInd>previousLbound&peakInd<rbound]
										curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-1
										
										maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
										pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
										curPeakList[pkID,]$isApex=1
										curPeakList[pkID,]$lboundInd=previousLbound
										curPeakList[pkID,]$rboundInd=rbound
										curPeakList[pkID,]$isShared<-1
									}else
									{#######if exceed the max window size, then simply keep the closed boundary as the final cut
										localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
										maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
										curPeakList[i,]$rboundInd<-rbound
										if(length(localPkInd)==1)
										{
											
											curPeakList[i,]$isApex=1
										}else############merge the peak apexes within current boundarys.
										{
											curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-1
											
											pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
											curPeakList[pkID,]$isApex=1
											curPeakList[pkID,]$lboundInd=lbound
											curPeakList[pkID,]$rboundInd=rbound
											curPeakList[pkID,]$isShared<-1
										}
									}
									preRbound<-rbound
									
								}else
								{
									########if lbound lower than rbound, then merge to next peak
									if((rbound-lbound)<=maxWindowlength)
									{
										curPeakList[i,]$isApex<--1
										
									}else
									{
										##if the range exceeds the max windowsize, then force 
										##the current rbound as the final boundary
										localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
										maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
										curPeakList[i,]$rboundInd<-rbound
										if(length(localPkInd)==1)
										{
											
											curPeakList[i,]$isApex=1
										}else############merge the peak apexes within current boundarys.
										{
											curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-1
											
											pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
											curPeakList[pkID,]$isApex=1
											curPeakList[pkID,]$lboundInd=lbound
											curPeakList[pkID,]$rboundInd=rbound
											curPeakList[pkID,]$isShared<-1
										}								
										preRbound<-rbound
									}
								}
							}
						}else
						{
							##########################################
							#if no rbound pass the pkheight ratio filter
							#then check the closed rbound
							#######################################
							rbound<-rightValleyIndVec[1]
							#if edge difference is not significant
							if(abs(curVecInt[rbound]-lboundHeight)/curPeakHight<=edgeHightDiffRatio)
							{
								curPeakList[i,]$rboundInd<-rbound
								preRbound<-rbound
								
								localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
								maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
								if(length(localPkInd)==1)
								{
									
									curPeakList[i,]$isApex=1
								}else############merge the peak apexes within current boundarys.
								{
									curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-1
									
									pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
									curPeakList[pkID,]$isApex=1
									curPeakList[pkID,]$lboundInd=lbound
									curPeakList[pkID,]$rboundInd=rbound
									curPeakList[pkID,]$isShared<-1
								}
							}else#if edge diff is significant,try to merge
							{
								########if lbound higher than rbound, then merge to previous peak
								if(maxlocalPkInd>0&&(lboundHeight>curVecInt[rbound]))
								{
									previousLbound<-curPeakList[curPeakList$pkInd==maxlocalPkInd,]$lboundInd		
									if((rbound-previousLbound)<=maxWindowlength)
									{
										localPkInd<-peakInd[peakInd>previousLbound&peakInd<rbound]
										curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-1
										
										maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
										pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
										curPeakList[pkID,]$isApex=1
										curPeakList[pkID,]$lboundInd=previousLbound
										curPeakList[pkID,]$rboundInd=rbound
										curPeakList[pkID,]$isShared<-1
									}else
									{#######if exceed the max window size, then simply keep the closed boundary as the final cut
										localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
										maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
										curPeakList[i,]$rboundInd<-rbound
										if(length(localPkInd)==1)
										{
											
											curPeakList[i,]$isApex=1
										}else############merge the peak apexes within current boundarys.
										{
											curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-1
											
											pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
											curPeakList[pkID,]$isApex=1
											curPeakList[pkID,]$lboundInd=lbound
											curPeakList[pkID,]$rboundInd=rbound
											curPeakList[pkID,]$isShared<-1
										}
									}
									preRbound<-rbound
								}else
								{
									########if lbound lower than rbound, then merge to next peak
									if((rbound-lbound)<=maxWindowlength)
									{
										curPeakList[i,]$isApex<--1
										
									}else
									{
										##if the range exceeds the max windowsize, then force 
										##the current rbound as the final boundary
										localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
										maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
										curPeakList[i,]$rboundInd<-rbound
										if(length(localPkInd)==1)
										{
											
											curPeakList[i,]$isApex=1
										}else############merge the peak apexes within current boundarys.
										{
											curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-1
											
											pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
											curPeakList[pkID,]$isApex=1
											curPeakList[pkID,]$lboundInd=lbound
											curPeakList[pkID,]$rboundInd=rbound
											curPeakList[pkID,]$isShared<-1
										}								
										preRbound<-rbound
									}
								}
							}
						}
						
					}else
					{
						##when the nearest rbound exceed the max window length, just pick the nearest one
						rbound<-min(valleyInd[valleyInd>curPeakInd])
						curPeakList[i,]$rboundInd<-rbound
						
						localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
						maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
						if(length(localPkInd)==1)
						{
							
							curPeakList[i,]$isApex=1
						}else############merge the peak apexes within current boundarys.
						{
							curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-1
							
							
							pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
							curPeakList[pkID,]$isApex=1
							curPeakList[pkID,]$lboundInd=lbound
							curPeakList[pkID,]$rboundInd=rbound
							curPeakList[pkID,]$isShared<-1
						}
						
						preRbound<-rbound
						
					}
				}
				
			}
			
			
			
		}
#		playwith({
#					plot(curVecInt,type="l")
#					points(curPeakInd,curPeakHight,col="red")
#					#points(rightValleyIndVec,curVecInt[rightValleyIndVec],col="green")
#					
#					
#					abline(v=rbound,col="blue")
#					abline(v=lbound,col="blue")
#				})
	}
	
#	curPeakList$ET<-(curPeakList$pkInd-1)*ScanInterval+delaytime
#	playwith({
#				
#				plot(((1:nPoints)-1)*ScanInterval+delaytime,curVecInt,type="l")
#				points(subset(curPeakList,isApex==1&isShared==1,select=c("ET","Intensity")),col="orange")
#				points(subset(curPeakList,isApex==1&isShared==0,select=c("ET","Intensity")),col="red")
#				points(subset(curPeakList,isApex==-1,select=c("ET","Intensity")),col="yellow")
#				points(subset(curPeakList,isApex==0,select=c("ET","Intensity")),col="green")
#				
#				points(valleyInd,curVecInt[valleyInd],col="pink")
#				
#				abline(v=(subset(curPeakList,isApex==1,select=c("lboundInd"))$lboundInd-1)*ScanInterval+delaytime,col="blue")
#				
#				abline(v=(subset(curPeakList,isApex==1,select=c("rboundInd"))$rboundInd-1)*ScanInterval+delaytime,col="blue")
#				
#			})
	##	##	
	
	
	
	
	cat(paste("finish peak picking\n",sep=""))
	
	if(!is.null(mz))
	{##EIC peaks
		
		curPeakList$curMz=mz
		##############################################
		#save the offset of the current EIC block position 
		#in the long intensity vector (vecInt) 
		#for each EIC peak,so that peak profile could be easily indexed 
		# from vecInt in decon
		################################################
		
		
		curPeakList$offset<-offset1
		curPeakList<-subset(curPeakList,isApex==1,select=c("curMz","pkInd","Intensity","lboundInd","rboundInd","isShared","offset"))
		
		##calculate the gss and sharpness for each EIC peak profile
		curPeakList$gss<-90
		curPeakList$gssPos<-0
		curPeakList$shrp<-0
		for(i in 1:nrow(curPeakList))
		{
			curPk<-curPeakList[i,]
			startInd<-curPk$lboundInd+offset1
			endInd<-curPk$rboundInd+offset1
			curProfile<-vecInt[startInd:endInd]
			
			gss<-GaussianDisimilarity1(curProfile)
			gssPos<-which.min(gss)
			curPeakList[i,]$gssPos<-gssPos
			curPeakList[i,]$gss<-gss[gssPos]
			
			shrp<-Sharpness1(yy=curProfile)	
			curPeakList[i,]$shrp<-shrp[gssPos]
			
		}
		curPeakList
		
	}else
	{##TIC peaks, only need to output isApex rows -- yan
		subset(curPeakList,isApex==1,select=c("pkInd","lboundInd","rboundInd","isShared","isApex"))
	}
	
}

decomposition_showPlotting<-function(inFilePath,params,cl,isDistParallel,clustingType)
{
	source(paste(params$codeDir,"pipeline.r",sep="/"))	
	library(gdata)
	library(gtools)
	library(cluster)
	library("ncdf")	
	options(expressions=1e5)
	DataFilelist<-params$DataFiles
	delaytime<-params$delaytime
	ScanInterval<-params$ScanInterval
	WorkDir<-params$WorkDir
	
	##load library
	lib<-readMSP2Spec(filename="../StdsLibrary/KQC.txt",withRT=T)
	
	fileName<-parseFileName(inFilePath)
	
	###############
	#read TIC data
	################		
	TICfile<-paste(WorkDir,"output/TIC/denoised_",fileName,"_TIC.cdf",sep="") 
	# The denoised TIC are saved in output/TIC/
	#if there's no denoised TIC file, then read raw cdf data
	TICfile<-ifelse(file.exists(TICfile),TICfile,inFilePath)
	ncid <- open.ncdf(TICfile)
	TIC <- get.var.ncdf(ncid, varid="total_intensity")
	close.ncdf(ncid)
	remove(ncid)
	
	
	###############
	#read EIC data
	################
	rawEICFile<-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
	denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
	EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
	
	if(file.exists(EICFile)==FALSE)next
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	##orignial mz vector from EIC data
	vectorMz<-get.var.ncdf(ncid, varid="mzVec")
	close.ncdf(ncid)
	remove(ncid)	
	
	###############################################
	#get decon window from TIC peak picking results
	###############################################
	cat(fileName,"reading TIC peak data...\n")
	TICPeakFile<-paste(params$WorkDir,"output/peakpicking/",fileName,"_TIC_PeakList.csv",sep="")
	TICApexList <- read.csv(TICPeakFile)
	TICApexList$RT<-(TICApexList$pkInd-1)*ScanInterval+delaytime
	TICpeaklist<-subset(TICApexList,isApex==1)
	
	winIDs<-1:nrow(TICpeaklist)
	
	#############################
	#read EIC peak picking result
	#############################
	cat(fileName,"reading EIC peak data...\n")
	EICPeakFile<-paste(WorkDir,"/output/peakpicking/",fileName,"PeakList.csv",sep="")
	if(file.exists(EICPeakFile)==FALSE)next
	EICpeaklist <- read.csv(EICPeakFile)
	colnames(EICpeaklist)[which(colnames(EICpeaklist)=="curMz")]<-"mz"
	
	#########################
	#start decon
	#########################
	cat(fileName,"start deconvolution...\n")
	
	###################
	#sequential version
	####################
	componentResults<-NULL
	specResults<-NULL
	mdlPkIDVec<-NULL
	EICpeaklist$flag=0
	totalscan<-length(vecInt)/length(vectorMz)
	
#	TICPkmean <- mean(TIC[TICpeaklist$pkInd])
#	minPkHeight<-TICPkmean*0.05
	playwith({
				
				plot((1:length(TIC)-1)*ScanInterval+delaytime,TIC,type="l")
				points(TICpeaklist$RT,TIC[TICpeaklist$pkInd],col="red")
				text(TICpeaklist$RT,TIC[TICpeaklist$pkInd],c(1:nrow(TICpeaklist)))
				abline(v=(TICpeaklist$lboundInd-1)*ScanInterval+delaytime,col="blue")
				abline(v=(TICpeaklist$rboundInd-1)*ScanInterval+delaytime,col="blue")				
			})	

	for(windowID in winIDs)
	{
		cat("window:",windowID,"\n")
		curTICPk<-TICpeaklist[windowID,]	
		minET=(curTICPk$lboundInd-1)*params$ScanInterval+params$delaytime
		maxET=(curTICPk$rboundInd-1)*params$ScanInterval+params$delaytime
		
		##decide the exlcuded ion for current TIC window
		##based on ET, for model peak selection
		if(minET>=8)
		{
			CurNonUmassVec<-c(params$NonUmassVec,51:100)
		}else
		{
			CurNonUmassVec<-params$NonUmassVec
		}
		
		##################################################
		#get all the EIC peaks within current decon window
		##################################################
		localTICApexList<-subset(TICApexList,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		#get the potential compound number for TIC peak apex list
		#there should be at least as many or more compounds than the number of apex 
		#minCmp<-nrow(localTICApexList)
		
		#get all EIC peaks within current decon window
		localEICPeakList<-subset(EICpeaklist,flag==0&pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		
		##keep the unique mass profile of local EIC list
		##original local EIC list consists of merged peaks
		localEICPeakList <- subset(localEICPeakList,isApex==1)
		
		if(nrow(localEICPeakList)>0)
		{
			minPkHeight <- mean(TIC[localTICApexList$pkInd])*0.025
			#minPkHeight <- ifelse(minPkHeight <=1000,1000,minPkHeight)
			#minPkHeight <- ifelse(minPkHeight >=100000,0.1*minPkHeight,minPkHeight)
			minPkHeight <- ifelse(minPkHeight <=800,800,minPkHeight)
			minPkHeight <- ifelse(minPkHeight >5000&&minPkHeight <=15000,5000+round((minPkHeight-5000)*0.5),minPkHeight)
			minPkHeight <- ifelse(minPkHeight >10000,10000+round((minPkHeight-10000)*0.1),minPkHeight)
#			
			allprofiles<-vector("list",nrow(localEICPeakList))
			for(i in 1:nrow(localEICPeakList))
			{
				##get intensity vector of current mz
				##from long intensity vector
				curEICpk<-localEICPeakList[i,]			
				startInd<-curEICpk$offset+curEICpk$lboundInd
				endInd<-curEICpk$offset+curEICpk$rboundInd
				
				##get EIC peak profile
				ET<-curEICpk$lboundInd:curEICpk$rboundInd
				curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd])
				allprofiles[[i]]<-curProfile
				names(allprofiles)[i]<-rownames(curEICpk)	
			}
			

			##Exclude those mass not abundent enough for models
			CandidatePeaklist<-subset(localEICPeakList,Intensity>minPkHeight)
			
			##Exclude those mass not eligible for models
			CandidatePeaklist<-subset(CandidatePeaklist,!mz%in%CurNonUmassVec)
			
			#################plot the profiles of all peaks
			colVec<-rainbow(nrow(CandidatePeaklist))
			names(colVec)<-CandidatePeaklist$mz
			labelVec<-paste("mz:",CandidatePeaklist$mz," ,gss:",format(CandidatePeaklist$gss,digits=2)," ,shrp:",format(CandidatePeaklist$shrp,digits=2))
			dataPoints<-subset(CandidatePeaklist,select=c("pkInd","Intensity"))
			dataPoints$pkInd<-(dataPoints$pkInd-1)*ScanInterval+delaytime
			playwith({
						plot(x=0,type = "n",xlim=(c(min(CandidatePeaklist$lboundInd),max(CandidatePeaklist$rboundInd))-1)*ScanInterval+delaytime,ylim=c(0,maxInt1))
						for(i in 1:nrow(CandidatePeaklist))
						{
							curMz<-as.character(CandidatePeaklist[i,]$mz)
							curID<-rownames(CandidatePeaklist[i,])
							points((allprofiles[[curID]]$ET-1)*ScanInterval+delaytime,allprofiles[[curID]]$int, type = "l",lwd=2,col=colVec[curMz])
							points((CandidatePeaklist[i,]$pkInd-1)*ScanInterval+delaytime,CandidatePeaklist[i,]$Intensity,col="green")
							
						}
					},data.points=dataPoints,labels=labelVec)
			##################
			
			##remove potentially shared ions by GSS and Sharpness
			if (minPkHeight <=5000) 
			{
				gaussianCutoff<-2.6
			}else 
			{
				gaussianCutoff<-2
			}
			
			#gaussianCutoff<-2.6
			ShrpCutoff<-5
			goodShapePeaklist<-subset(CandidatePeaklist,gss<=gaussianCutoff&shrp>=ShrpCutoff)
			
#			if (nrow(goodShapePeaklist)>10)
#			{
#				importMasslist <- goodShapePeaklist[order(goodShapePeaklist[,"Intensity"],decreasing=T),][c(1:5),]
#				ImportantLowMass <- unique(subset(importMasslist,gss<=2)$mz)
#				CurNonUmassVec<- CurNonUmassVec[!CurNonUmassVec%in%ImportantLowMass]
#				# remove other low and non-unique mass
#				goodShapePeaklist<-subset(goodShapePeaklist,!mz%in%CurNonUmassVec)
#			}
#			

			#####################plot the profiles of selected peaks with good gaussian shape
			#				minLbound<-min(goodShapePeaklist$lboundInd)
			#				maxRbound<-max(goodShapePeaklist$rboundInd)	
			#				
			#					colVec<-rainbow(nrow(goodShapePeaklist))
			#					labelVec<-paste("pkID:",rownames(goodShapePeaklist),"mz:",goodShapePeaklist$mz," ,gss:",format(goodShapePeaklist$gss,digits=2)," ,shrp:",format(goodShapePeaklist$shrp,digits=2))
			#					dataPoints<-subset(goodShapePeaklist,select=c("pkInd","Intensity"))
			#					dataPoints$pkInd<-(dataPoints$pkInd-1)*ScanInterval+delaytime
			#					playwith({
			#					plot(x=0,type = "n",xlim=((c(minLbound,maxRbound)-1)*ScanInterval+delaytime),ylim=c(0,max(maxInt2)))
			#					for(i in 1:nrow(goodShapePeaklist))
			#					{
			#						curID<-rownames(goodShapePeaklist[i,])
			#						points((profiles[[curID]]$ET-1)*ScanInterval+delaytime,profiles[[curID]]$int, type = "l",col=colVec[i])
			#						
			#					}
			#					},data.points=dataPoints,labels=labelVec)
			#####################
			
			nPks<-nrow(goodShapePeaklist)		
			modelPkList<-NULL
			isMultGroups<-FALSE
			mdlPkIDs<-NULL
			
			if(nPks>0)
			{
				##get good peak profiles
				profiles<-vector("list",nrow(goodShapePeaklist))
				for(i in 1:nrow(goodShapePeaklist))
				{
					##get intensity vector of current mz
					##from long intensity vector
					curEICpk<-goodShapePeaklist[i,]					
					startInd<-curEICpk$offset+curEICpk$lboundInd
					endInd<-curEICpk$offset+curEICpk$rboundInd
					
					##get EIC profile
					ET<-curEICpk$lboundInd:curEICpk$rboundInd
					curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd])
					profiles[[i]]<-curProfile
					names(profiles)[i]<-rownames(curEICpk)	
				}
				
				##update each profile by creating mirror images
				##chosse the one with better GSS
				for(i in 1:length(profiles))
				{
					curPkId<-names(profiles[i])
					curPk<-goodShapePeaklist[curPkId,]
					pos<-curPk$gssPos
					profiles[[i]]<-copyShape(p=profiles[[i]],curPk$pkInd,from=pos)
					goodShapePeaklist[curPkId,]$lboundInd<-profiles[[i]]$ET[1]
					goodShapePeaklist[curPkId,]$rboundInd<-profiles[[i]]$ET[nrow(profiles[[i]])]
				}
				
				#############################
				#Hierachical clustering
				#decide the number of groups 
				#within current decon window
				#############################
				if(nPks>=2)
				{
					##save EIC profies in global variable for 
					##broadcasting to slave nodes for parallel computing
					assign("Global.curProfs", value=profiles, envir = .GlobalEnv)
					
					##calculate pairwise distance among the EICs				
					if(isDistParallel)
					{
						r<-parDistCal4(cl,isUnion=F)##parallel verison
					}else
					{
						r<-DistCal4(isUnion=F)#non-parallel version
					}	
					
					##convert distance matrix to triangle matrix
					distance <- as.dist(r)
					maxIntraDist<- 15 #cutoff to decide if split
					
					if(clustingType=="h")
					{
						##################################
						#hierarchical clustering and cut 
						#into groups by the height
						##################################
						clustResut<-hclust(distance)
						FinalClustResut<-cutree(clustResut,h=maxIntraDist)
						Clusters<-unique(FinalClustResut)
						if(length(Clusters)>1)isMultGroups<-T
						
						##########################plot the profiles of cluster groups
#							colVec<-rainbow(length(Clusters))
#							labelVec<-paste("pkID:",rownames(goodShapePeaklist),"mz:",goodShapePeaklist$mz," ,gss:",format(goodShapePeaklist$gss,digits=2)," ,shrp:",format(goodShapePeaklist$shrp,digits=2))
#							playwith({
#							plot(x=0,type = "n",xlim=(c(minLbound,maxRbound)-1)*ScanInterval+delaytime,ylim=c(0,max(maxInt2)))
#							for(i in 1:length(Clusters))
#							{
#								######get current cluster ID
#								curCluster<-Clusters[i]
#								#####get masses belong to current cluster
#								curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
#								
#								for(j in 1:length(curGroup))
#								{
#									
#									points((profiles[[curGroup[j]]]$ET-1)*ScanInterval+delaytime,profiles[[curGroup[j]]]$int, type = "l",col=colVec[i])
#								}
#							}
#							},data.points=subset(goodShapePeaklist,select=c("pkInd","Intensity")),labels=labelVec)
						
					##########################	
						
					}else
					{	###############################################################
						#partitional clustering and split into groups by the silehoute
						###############################################################
						if(max(distance)<=maxIntraDist)
						{###if the max distance is already small enought, no need to deconvolute
							isMultGroups<-FALSE
						}else
						{
							
							if(nPks>2)
							{
								########################################
								#try different k to decide the cluster number by intra dist
								#######################################
								
								maxClustNumber<-min(20,nPks-1)
								asw<-numeric(maxClustNumber)
								for (k in min(maxClustNumber,max(2,minCmp)):maxClustNumber)
								{
									clustResut <-pam(distance,k)
									Clusters<-unique(clustResut$clustering)								
									isSplit<-F
									
									for(curCluster in Clusters)
									{
										#####get masses belong to current cluster
										curGroup<-names(clustResut$clustering[clustResut$clustering==curCluster])
										#######does not average the single-ion cluster which has zero intra dist
										if(length(curGroup)>1)
										{#######calculate the average intra distance to see if the current group is a valid single group of fragments
											
											if(max(as.dist(r[curGroup,curGroup]))>maxIntraDist)
											{##if current group is not valid, skip to the next K,
												#and the current k value 
												isSplit<-T
												break
											}
											
										}
									}
									if(!isSplit)
										asw[k]<-clustResut$silinfo $ avg.width
								}
								k.best <- which.max(asw)
								if(k.best>1)
								{
									FinalClustResut<-pam(distance,k.best,cluster.only=TRUE,stand=TRUE)
									
									#################################################
									#calculate the average intra dist in each cluster
									###############################################
									isTooBigIntraDist<-NULL
									avgIntraDistVec<-NULL
									Clusters<-unique(FinalClustResut)
									for(curCluster in Clusters)
									{
										#####get masses belong to current cluster
										curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
										
										IntraDist<-mean(as.dist(r[curGroup,curGroup]))
										isTooBigIntraDist<-c(isTooBigIntraDist,ifelse(IntraDist>6,T,F))
										if(length(curGroup)>1)#######does not average the single-ion cluster which has zero intra dist
											avgIntraDistVec<-c(avgIntraDistVec,IntraDist)
										
									}
									#################################################################
									#if after splitting,intra dist is smaller than before splitting
									#then accept the best clustering schema, 
									#and earch for unique mass within each cluster
									###############################################################
									if(mean(avgIntraDistVec)<mean(distance))
										isMultGroups<-TRUE
								}
							}else#if ==2, then directly set flag=T
							{
								isMultGroups<-TRUE
							}
						}
					}
				}
				
				
				
				
				######################
				#Model peak selection 
				######################
				
				modelPkList<-NULL
				
				# single group
				if (!isMultGroups && nPks >=1)
				{
					if (nPks==1) {
						modelPkList<-goodShapePeaklist
					} else 
					{
						##choose best one as model peak
						mdlPkCandidates<-goodShapePeaklist
						
						##scoring system
						#mass score
						mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
						
						##gassian similarity score
						mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
						
						##peak height score
						mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))					
						mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
						mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
						mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
						
						##linear combination of 3 scores
						scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
						#scores<-(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3+0.4*mdlPkCandidates$f4)
						
						##rescale score to (0,1000)
						mdlPkCandidates$score<-scale(scores,F,0.001)
						
						##select the max score (>500)
						mdlPkCandidates <- subset(mdlPkCandidates,score>=500)
						if (nrow(mdlPkCandidates)>0)
						{
							modelPkList<- mdlPkCandidates[which.max(mdlPkCandidates$score),]
						}
					}
				}
				
				
				if(isMultGroups)
				{
					if(nPks==2)
					{
						##only two candidates, simply assign them as model peaks
						#modelPkVec<-rownames(goodShapePeaklist)
						modelPkList<-goodShapePeaklist
					}else
					{
						nMinFragment<-1
						if(length(Clusters)>0)
						{	
							##select model peak from each group[!isTooBigIntraDist]
							for(curCluster in Clusters)
							{
								##get masses belong to current cluster
								curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
								
								##get model peak candidates 
								mdlPkCandidates<-goodShapePeaklist[curGroup,]
								
								################################plot the profiles of cluster groups
#									colVec<-rainbow(length(curGroup))
#									labelVec<-paste("pkID:",rownames(mdlPkCandidates),"mz:",mdlPkCandidates$mz," ,gss:",format(mdlPkCandidates$gss,digits=2)," ,shrp:",format(mdlPkCandidates$shrp,digits=2))
#									playwith({
#												plot(x=0,type = "n",xlim=c(minLbound,maxRbound),ylim=c(0,maxInt2))
#												
#												for(j in 1:length(curGroup))
#												{
#													
#													points(profiles[[curGroup[j]]], type = "l",col=colVec[j])
#												}
#												
#											},data.points=subset(mdlPkCandidates,select=c("pkInd","Intensity")),labels=labelVec)
								###############################
								
								if(nrow(mdlPkCandidates)>0)
								{
									##mass score
									mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
									##gassian similarity score
									mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
									##peak height score
									mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))								
									
									##assign zero value when each factor is NA
									mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
									mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
									mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
									
									##linear combination of 3 scores
									scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
									
									##rescale score to (0,1000)
									mdlPkCandidates$score<-scale(scores,F,0.001)
									
									##select the max score
									currentMdlPk<-mdlPkCandidates[which.max(mdlPkCandidates$score),]
									modelPkList<-rbind(modelPkList,currentMdlPk)
								}
							}
						}
					}			
				}
				
				#####################plot the profiles of model peaks
#					colVec<-rainbow(nrow(modelPkList))
#					
#					labelVec<-paste("pkid:",rownames(modelPkList),"mz:",modelPkList$mz," ,gss:",format(modelPkList$gss,digits=2)," ,shrp:",format(modelPkList$shrp,digits=2))
#				#					labelVec<-paste("gss=",format(modelPkList$gss,digits=2),"\n shrp=",format(modelPkList$shrp,digits=2))
#					dataPoints<-subset(modelPkList,select=c("pkInd","Intensity"))
#					dataPoints$pkInd<-(dataPoints$pkInd-1)*ScanInterval+delaytime
#					
#					playwith({
#					plot(x=0,type = "n",xlim=(c(minLbound,maxRbound)-1)*ScanInterval+delaytime,ylim=c(0,max(maxInt2)))
#					for(i in 1:nrow(modelPkList))
#					{
#						curPkid<-rownames(modelPkList)[i]
#						
#						points((profiles[[curPkid]]$ET-1)*ScanInterval+delaytime,profiles[[curPkid]]$int, type = "l",col=colVec[i])
#					}
#					},data.points=dataPoints,labels=labelVec)
				####################	
			}
			
			
			###################################
			#decomposition based on model peaks
			###################################
			
			decom = 1
			while (decom==1) 
			{
				modelPkDistance <- NULL
				specList<-NULL
				componentList<-NULL
				mdlPkIDs <- NULL
				
				##No co-eluting compounds (only one or no model peak detected)
				##Then extract spectrum directly
				if(nrow(modelPkList)<=1||is.null(modelPkList))
				{
					##add model peak to the spectrum
					specList[[1]]<-data.frame(mz=as.integer(),int=as.integer())
					
					##No model peak detected,then simply extract 
					##spectrum from the TIC peak apex time
					if(nrow(modelPkList)==0||is.null(modelPkList))
					{	
						##get 1st mass				
						curMz<-as.character(localEICPeakList$mz[1])
						
						##get TIC peak apex time and boundaries
						maxInd<-curTICPk$pkInd
						curLbound<-curTICPk$lboundInd
						curRbound<-curTICPk$rboundInd
						
						##Init component list
						componentList[[1]]<-data.frame(mz=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),pkInd=as.integer(),Intensity=as.integer(),
								area=as.integer(),windowID=as.integer(),compoundID=as.integer(),isModel=as.integer(),isUnique=as.integer(),factor=as.numeric())
						
						##name spec list with time and mz
						clustID<-format((maxInd-1)*ScanInterval+delaytime,digit=4)
						names(specList)[1]<-paste(clustID,curMz)
						
						##add rest peaks from the local window 
						##to component list and spectral list
						for(i in 1:nrow(localEICPeakList))
						{
							curPkID<-rownames(localEICPeakList)[i]					
							curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
							curProfile1$int[is.na(curProfile1$int)]<-0
							curArea<-sum(curProfile1$int)
							
							curMz<-as.character(localEICPeakList$mz[i])
							curInt<-subset(curProfile1,ET==maxInd)$int[1]
							
							if(!is.na(curInt)&&curInt>0)
							{
								specList[[1]]<-rbind(specList[[1]],data.frame(mz=as.integer(curMz),int=curInt))#add umass to the spectrum					
								componentList[[1]]<-rbind(componentList[[1]],data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,
												pkInd=maxInd,Intensity=curInt,area=curArea,windowID=windowID,compoundID=1,isModel=0,isUnique=0,factor=0))
							}
						}
						
						##Mark the ion with maximal intensity as model
						componentList[[1]][which.max(componentList[[1]]$Intensity),]$isModel=1
						
					}else
					{
						##single one model peak detected,then extract 
						##spectrum from the model peak apex time
						##Use model peak apex,boundary index
						curMz<-modelPkList$mz[1]
						maxInd<-modelPkList$pkInd[1]
						curLbound<-modelPkList$lboundInd[1]
						curRbound<-modelPkList$rboundInd[1]
						curPhHeight<-modelPkList$Intensity[1]	
						curPkID<-rownames(modelPkList)[1]
						
						##model peaks could not exist in all profiles which save all peaks with isApex=1
						##sometimes a boundary peak of neighboring merged peak is inclued because of its 
						##peak index is inside the decon range
						curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),profiles[[curPkID]],by.y="ET",all.x=T)
						curProfile1$int[is.na(curProfile1$int)]<-0
						curArea<-sum(curProfile1$int)
						componentList[[1]]<-data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,pkInd=maxInd,
								Intensity=curPhHeight,area=curArea,windowID=windowID,compoundID=1,isModel=1,isUnique=0,factor=0)
						
						#create spec list as NIST format for future identifiation
						clustID<-format((maxInd-1)*ScanInterval+delaytime,digit=4)
						names(specList)[1]<-paste(clustID,curMz)						
						
						##peak list id to extract spectrum
						peakListIndex <-c(1:nrow(localEICPeakList))				
						if (curPkID%in%rownames(localEICPeakList)) {
							##remove already detected model peaks
							peakListIndex <- peakListIndex[-which(rownames(localEICPeakList)==curPkID)]
						}
						
						##add rest peaks from the local window 
						##to component list and spectral list
						for(i in peakListIndex)
						{
							curPkID<-rownames(localEICPeakList)[i]					
							curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
							curProfile1$int[is.na(curProfile1$int)]<-0
							curArea<-sum(curProfile1$int)						
							curMz<-as.character(localEICPeakList$mz[i])
							curInt<-subset(curProfile1,ET==maxInd)$int[1]
							
							if(!is.na(curInt)&&curInt>0)
							{
								specList[[1]]<-rbind(specList[[1]],data.frame(mz=as.integer(curMz),int=curInt))#add umass to the spectrum					
								componentList[[1]]<-rbind(componentList[[1]],data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,
												pkInd=maxInd,Intensity=curInt,area=curArea,windowID=windowID,compoundID=1,isModel=0,isUnique=0,factor=0))
							}
						}						
					}
					
					componentList[[1]][which.max(componentList[[1]]$Intensity),]$isUnique=1
					##replace decom as 0, stop while loop
					decom=0
				}else
				{
					##More than one model peaks
					##Decompostion
					
					#######################################
					#Makeup 1 for splitting issue
					#Checking pairwise model peaks distance
					#######################################
					
					modelPkDistance <- as.matrix(modelPkDist(modelPkList,profiles))
					rmmodelPk <- NULL
					ModelPkNames <- rownames(modelPkDistance)
					
					##model same matrix as modelPKDistance with index
					psmatrix <- matrix(c(1:length(modelPkDistance)),ncol=nrow(modelPkDistance),byrow=T)
					psmodel <- which(modelPkDistance<=4&modelPkDistance>0)
					
					##get index of model peaks with smaller distance
					for (i in psmodel) {
						pos <- matrix.index(psmatrix, i)
						rmmodelPk <- c(rmmodelPk,ModelPkNames[pos[1]],ModelPkNames[pos[2]])
					}
					rmmodelPks <- unique(rmmodelPk)
					
					rmmodelPkList <- modelPkList[rmmodelPks,]
					modelPkList1 <- rmmodelPkList[which.max(rmmodelPkList$score),]
					modelPkList <- modelPkList[!rownames(modelPkList)%in%rmmodelPks,]
					
					##Only keep the one with highest score
					modelPkList <- rbind(modelPkList,modelPkList1)				
					mdlPkIDs<-rownames(modelPkList)
					
					##only one model peak after testing
					##GO back to the upper step
					if (nrow(modelPkList)==1) decom=1
					
					else {
						##still multiple models
						lbound<-min(modelPkList$lboundInd)
						rbound<-max(modelPkList$rboundInd)
						
						#################################
						#get all EIC profiles within the window 
						#X matrix format:column-EIC,row-scan
						#############################################
						
						X<-NULL
						mzVec<-unique(localEICPeakList$mz)
						nMz<-length(mzVec)
						for(i in 1:nMz)
						{
							mz<-mzVec[i]
							mzInd<-which(vectorMz==mz)
							startInd<-(mzInd-1)*totalscan+1
							endInd<-startInd+totalscan-1
							curVecInt <- vecInt[startInd:endInd]
							X<-cbind(X,curVecInt[lbound:rbound])
						}
						colnames(X)<-mzVec	
						nScans<-nrow(X)
						
						###########################################
						#Get model peak profiles as S matrix
						#each column is the EIC for one model peak
						###########################################
						
						S<-NULL
						modelPkList$area<-0
						for(i in 1:nrow(modelPkList))
						{
							curProfile1<-merge(data.frame(ET=as.integer(lbound:rbound)),profiles[[rownames(modelPkList[i,])]],by.y="ET",all.x=T)
							
							##expand the profile by filling zero values to those uncovered area
							curProfile1$int[is.na(curProfile1$int)]<-0
							S<-cbind(S,curProfile1$int)
							modelPkList[i,]$area<-sum(curProfile1$int)		
							
							##init the component & spec list	
							componentList[[i]]<-data.frame(mz=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),pkInd=as.integer(),Intensity=as.integer(),
									area=as.integer(),windowID=as.integer(),compoundID=as.integer(),isModel=as.integer(),isUnique=as.integer(),factor=as.numeric())
							##add model peak to the spectrum
							specList[[i]]<-data.frame(mz=as.integer(),int=as.integer())
							curRT<-format((modelPkList$pkInd[i]-1)*ScanInterval+delaytime,digit=4)
							names(specList)[i]<-paste(curRT,modelPkList$mz[i])
						}					
						colnames(S)<-mdlPkIDs
						
						##Residual minimization
						for(i in 1:ncol(X))
						{
							curMz<-colnames(X)[i]						
							M<-X[,i]#mixture signal
							
							srcIDs<-mdlPkIDs
							A<-optim(par=rep(0,ncol(S)),fn=fr, gr = NULL,M,S,lower = 0, method="L-BFGS-B")$par
							A<-as.matrix(A)
							
							##restore each ion signals and save them into 
							##component list ans speclist
							for(m in 1:length(srcIDs))
							{
								curmdlPkID<-srcIDs[m]
								j<-which(mdlPkIDs==curmdlPkID)
								
								if(A[m]>0)
								{
									specList[[j]]<-rbind(specList[[j]],data.frame(mz=as.integer(curMz),int=A[m,]))###add the curent restored mixing factor to spec
									if (modelPkList[curmdlPkID,]$mz==curMz) {
										componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,
														rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=as.integer(modelPkList[curmdlPkID,]$Intensity*A[m,]),
														area=as.integer(modelPkList[curmdlPkID,]$area*A[m,]),windowID=windowID,compoundID=j,isModel=1,isUnique=0,factor=A[m]))
									}
									
									else {
										componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,
														rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=as.integer(modelPkList[curmdlPkID,]$Intensity*A[m,]),
														area=as.integer(modelPkList[curmdlPkID,]$area*A[m,]),windowID=windowID,compoundID=j,isModel=0,isUnique=0,factor=A[m]))
									}
								}							
							}
						}
						
						###############################################
						#Makeup 2 for splitting issue
						#Checking pairwise resovlved spectra similarity
						###############################################
						
						RepeatModelPk <- NULL
						Num <- nrow(modelPkList)
						indexPairs<-combinations(Num,2)
						for(i in 1:nrow(indexPairs))
						{
							index <- indexPairs[i,]
							index1 <- index[1]
							index2 <- index[2]
							
							##calculate time difference between two model peaks
							modelPKET <- modelPkList$pkInd[c(index1,index2)]
							
							score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
							if (score >=970||(abs(modelPKET[1]-modelPKET[2]) <=6&&score >=900))
							{
								selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index2],mdlPkIDs[index1])
								RepeatModelPk <- c(RepeatModelPk,selectedModelPkID)
							}
						}
						
						##select unique filtered model peaks 
						RepeatModelPk <- unique(RepeatModelPk)
						modelPkList <- modelPkList[!rownames(modelPkList)%in%RepeatModelPk,]
						
						##The modelPk number keep the same, do not RE-decompose
						if (nrow(modelPkList) == Num) {
							decom = 0
							
							if (nrow(modelPkList) == 2) {
								mzlist1 <- subset(componentList[[1]],Intensity>=1000)$mz
								mzlist2 <- subset(componentList[[2]],Intensity>=1000)$mz
								SharedMS1 <- unique(c(mzlist1[mzlist1%in%mzlist2],CurNonUmassVec))
								SharedMS2 <- unique(c(mzlist2[mzlist2%in%mzlist1],CurNonUmassVec))
								
								unique_comp1 <- componentList[[1]][!componentList[[1]]$mz%in%SharedMS1,]
								unique_comp2 <- componentList[[2]][!componentList[[2]]$mz%in%SharedMS2,]
								
								componentList[[1]][componentList[[1]]$mz==unique_comp1[which.max(unique_comp1$Intensity),]$mz,]$isUnique=1
								componentList[[2]][componentList[[2]]$mz==unique_comp2[which.max(unique_comp2$Intensity),]$mz,]$isUnique=1
								
							}else
							{
								for (i in c(1:nrow(modelPkList)))
								{
									if (i==1||i==nrow(modelPkList))
									{
										if (i==1) {
											mzlist1 <- subset(componentList[[1]],Intensity>=1000)$mz
											mzlist2 <- subset(componentList[[2]],Intensity>=1000)$mz
											NonUniqueMS1 <- unique(c(mzlist1[mzlist1%in%mzlist2],CurNonUmassVec))						
											unique_comp1 <- componentList[[1]][!componentList[[1]]$mz%in%NonUniqueMS1,]				
											componentList[[1]][componentList[[1]]$mz==unique_comp1[which.max(unique_comp1$Intensity),]$mz,]$isUnique=1
										}else{
											mzlist1 <- subset(componentList[[i]],Intensity>=1000)$mz
											mzlist2 <- subset(componentList[[i-1]],Intensity>=1000)$mz
											NonUniqueMS1 <- unique(c(mzlist1[mzlist1%in%mzlist2],CurNonUmassVec))						
											unique_comp1 <- componentList[[i]][!componentList[[i]]$mz%in%NonUniqueMS1,]				
											componentList[[i]][componentList[[i]]$mz==unique_comp1[which.max(unique_comp1$Intensity),]$mz,]$isUnique=1
										}
									}else{
										mzlist1 <- subset(componentList[[i-1]],Intensity>=1000)$mz
										mzlist2 <- subset(componentList[[i]],Intensity>=1000)$mz
										mzlist3 <- subset(componentList[[i+1]],Intensity>=1000)$mz
										
										NonUniqueMS2 <- unique(c(mzlist2[mzlist2%in%mzlist1],mzlist2[mzlist2%in%mzlist3],CurNonUmassVec))					
										unique_comp2 <- componentList[[i]][!componentList[[i]]$mz%in%NonUniqueMS2,]				
										componentList[[i]][componentList[[i]]$mz==unique_comp2[which.max(unique_comp2$Intensity),]$mz,]$isUnique=1
									}						
								}
							}		
						}			
					} # else model PK num is >2 after testing distance	
				} # model PK num is >2 at first
			} # while statement
			
			##flag analyzed EIC peaks 
			EICpeaklist[rownames(localEICPeakList),]$flag<-1
			componentResults<-c(componentResults,componentList)
			specResults<-c(specResults,specList)
			
			##model peak id
			mdlPkIDVec<-c(mdlPkIDVec,mdlPkIDs)		
		}	
	}	
	
	componentResults<-do.call(rbind,componentResults)	
	cat(fileName,"finished deconvolution...\n")
	
	##################################
	#output lib matching results
	###################################
	
	##similarity cutoff=600
	libmatch <- libMatching(lib,specResults,600)
	write.csv(libmatch,file=paste(params$WorkDir,"/output/decon/",fileName,"_libmatch.csv",sep=""))
	
	##############################################################
	#collect spectrum result in dataframe for reverse lib matching
	##############################################################
	
	##exp spectra as lib, get top 1 hit for each std
	evaluateIdentification(refSpec=lib,inputSpec=specResults,RT_Tolerance=70,minSpecSimilarity=650,fileName,withRT=T,params)
	cat(fileName,"finish matching...\n")
	write.csv(mdlPkIDVec,file=paste(params$WorkDir,"/output/decon/",fileName,"mdlpklist.csv",sep=""))
	
	################################################
	#collect component result in dataframe for alignment
	###################################################
	
	##filter out zero-int fragments
	componentResults<-subset(componentResults,Intensity>0,select=c("mz","lboundInd","rboundInd","pkInd","Intensity","area","windowID","compoundID","isModel","isUnique","factor"))
	write.csv(componentResults,row.names=FALSE,file=paste(params$WorkDir,"/output/decon/",fileName,"Decon.csv",sep=""))
	cat(fileName,"writing component results...\n")
}

####################################################
# Before implementing StN and Sig level calculation
# Intensity cutoff for good peaks selection
# Oct 12, 2011
####################################################

decomposition_Oct12<-function(inFilePath,params,cl,isDistParallel,clustingType)
{
	source(paste(params$codeDir,"pipeline.r",sep="/"))	
	library(gdata)
	library(gtools)
	library(cluster)
	library("ncdf")	
	options(expressions=1e5)
	DataFilelist<-params$DataFiles
	delaytime<-params$delaytime
	ScanInterval<-params$ScanInterval
	WorkDir<-params$WorkDir
	
	##load library
	lib<-readMSP2Spec(filename="../StdsLibrary/KQC.txt",withRT=T)
	
	fileName<-parseFileName(inFilePath)
	
	###############
	#read TIC data
	################		
	TICfile<-paste(WorkDir,"output/TIC/denoised_",fileName,"_TIC.cdf",sep="") 
	# The denoised TIC are saved in output/TIC/
	#if there's no denoised TIC file, then read raw cdf data
	TICfile<-ifelse(file.exists(TICfile),TICfile,inFilePath)
	ncid <- open.ncdf(TICfile)
	TIC <- get.var.ncdf(ncid, varid="total_intensity")
	close.ncdf(ncid)
	remove(ncid)
	
	
	###############
	#read EIC data
	################
	rawEICFile<-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
	denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
	EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)
	
	if(file.exists(EICFile)==FALSE)next
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	##orignial mz vector from EIC data
	vectorMz<-get.var.ncdf(ncid, varid="mzVec")
	close.ncdf(ncid)
	remove(ncid)	
	
	###############################################
	#get decon window from TIC peak picking results
	###############################################
	cat(fileName,"reading TIC peak data...\n")
	TICPeakFile<-paste(params$WorkDir,"output/peakpicking/",fileName,"_TIC_PeakList.csv",sep="")
	TICApexList <- read.csv(TICPeakFile)
	TICApexList$RT<-(TICApexList$pkInd-1)*ScanInterval+delaytime
	TICpeaklist<-subset(TICApexList,isApex==1)
	
	winIDs<-1:nrow(TICpeaklist)
	
	#############################
	#read EIC peak picking result
	#############################
	cat(fileName,"reading EIC peak data...\n")
	EICPeakFile<-paste(WorkDir,"/output/peakpicking/",fileName,"PeakList.csv",sep="")
	if(file.exists(EICPeakFile)==FALSE)next
	EICpeaklist <- read.csv(EICPeakFile)
	colnames(EICpeaklist)[which(colnames(EICpeaklist)=="curMz")]<-"mz"
	
	#########################
	#start decon
	#########################
	cat(fileName,"start deconvolution...\n")
	
	###################
	#sequential version
	####################
	componentResults<-NULL
	specResults<-NULL
	mdlPkIDVec<-NULL
	EICpeaklist$flag=0
	totalscan<-length(vecInt)/length(vectorMz)
	modelPkList3s <-NULL
	modelPkList2s <-NULL
	modelPkList1s <-NULL
#	TICPkmean <- mean(TIC[TICpeaklist$pkInd])
#	minPkHeight<-TICPkmean*0.05
	
	for(windowID in winIDs)
	{
		cat("window:",windowID,"\n")
		curTICPk<-TICpeaklist[windowID,]	
		minET=(curTICPk$lboundInd-1)*params$ScanInterval+params$delaytime
		maxET=(curTICPk$rboundInd-1)*params$ScanInterval+params$delaytime
		
		##decide the exlcuded ion for current TIC window
		##based on ET, for model peak selection
		if(minET>=8)
		{
			CurNonUmassVec<-c(params$NonUmassVec,51:100)
		}else
		{
			CurNonUmassVec<-params$NonUmassVec
		}
		
		##################################################
		#get all the EIC peaks within current decon window
		##################################################
		localTICApexList<-subset(TICApexList,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		#get the potential compound number for TIC peak apex list
		#there should be at least as many or more compounds than the number of apex 
		#minCmp<-nrow(localTICApexList)
		
		#get all EIC peaks within current decon window
		localEICPeakList<-subset(EICpeaklist,flag==0&pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		
		##keep the unique mass profile of local EIC list
		##original local EIC list consists of merged peaks
		localEICPeakList <- subset(localEICPeakList,isApex==1)
		
		if(nrow(localEICPeakList)>0)
		{
			minPkHeight <- mean(TIC[localTICApexList$pkInd])*0.025
			#minPkHeight <- ifelse(minPkHeight <=1000,1000,minPkHeight)
			#minPkHeight <- ifelse(minPkHeight >=100000,0.1*minPkHeight,minPkHeight)
			minPkHeight <- ifelse(minPkHeight <=800,800,minPkHeight)
			minPkHeight <- ifelse(minPkHeight >5000&&minPkHeight <=15000,5000+round((minPkHeight-5000)*0.5),minPkHeight)
			minPkHeight <- ifelse(minPkHeight >10000,10000+round((minPkHeight-10000)*0.1),minPkHeight)
#			
			allprofiles<-vector("list",nrow(localEICPeakList))
			for(i in 1:nrow(localEICPeakList))
			{
				##get intensity vector of current mz
				##from long intensity vector
				curEICpk<-localEICPeakList[i,]			
				startInd<-curEICpk$offset+curEICpk$lboundInd
				endInd<-curEICpk$offset+curEICpk$rboundInd
				
				##get EIC peak profile
				ET<-curEICpk$lboundInd:curEICpk$rboundInd
				curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd])
				allprofiles[[i]]<-curProfile
				names(allprofiles)[i]<-rownames(curEICpk)	
			}
			
			##Exclude those mass not abundent enough for models
			CandidatePeaklist<-subset(localEICPeakList,Intensity>minPkHeight)
			
			##Exclude those mass not eligible for models
			CandidatePeaklist<-subset(CandidatePeaklist,!mz%in%CurNonUmassVec)
			
			##remove potentially shared ions by GSS and Sharpness
			if (minPkHeight <=5000) 
			{
				gaussianCutoff<-2.6
			}else 
			{
				gaussianCutoff<-2
			}
			
			#gaussianCutoff<-2.6
			ShrpCutoff<-5
			goodShapePeaklist<-subset(CandidatePeaklist,gss<=gaussianCutoff&shrp>=ShrpCutoff)
			
#			if (nrow(goodShapePeaklist)>10)
#			{
#				importMasslist <- goodShapePeaklist[order(goodShapePeaklist[,"Intensity"],decreasing=T),][c(1:5),]
#				ImportantLowMass <- unique(subset(importMasslist,gss<=2)$mz)
#				CurNonUmassVec<- CurNonUmassVec[!CurNonUmassVec%in%ImportantLowMass]
#				# remove other low and non-unique mass
#				goodShapePeaklist<-subset(goodShapePeaklist,!mz%in%CurNonUmassVec)
#			}
#			
			
			nPks<-nrow(goodShapePeaklist)		
			modelPkList<-NULL
			isMultGroups<-FALSE
			mdlPkIDs<-NULL
			
			if(nPks>0)
			{
				##get good peak profiles
				profiles<-vector("list",nrow(goodShapePeaklist))
				for(i in 1:nrow(goodShapePeaklist))
				{
					##get intensity vector of current mz
					##from long intensity vector
					curEICpk<-goodShapePeaklist[i,]					
					startInd<-curEICpk$offset+curEICpk$lboundInd
					endInd<-curEICpk$offset+curEICpk$rboundInd
					
					##get EIC profile
					ET<-curEICpk$lboundInd:curEICpk$rboundInd
					curProfile<-data.frame(ET=ET,int=vecInt[startInd:endInd])
					profiles[[i]]<-curProfile
					names(profiles)[i]<-rownames(curEICpk)	
				}
				
				##update each profile by creating mirror images
				##chosse the one with better GSS
				for(i in 1:length(profiles))
				{
					curPkId<-names(profiles[i])
					curPk<-goodShapePeaklist[curPkId,]
					pos<-curPk$gssPos
					profiles[[i]]<-copyShape(p=profiles[[i]],curPk$pkInd,from=pos)
					goodShapePeaklist[curPkId,]$lboundInd<-profiles[[i]]$ET[1]
					goodShapePeaklist[curPkId,]$rboundInd<-profiles[[i]]$ET[nrow(profiles[[i]])]
				}
				
				#############################
				#Hierachical clustering
				#decide the number of groups 
				#within current decon window
				#############################
				if(nPks>=2)
				{
					##save EIC profies in global variable for 
					##broadcasting to slave nodes for parallel computing
					assign("Global.curProfs", value=profiles, envir = .GlobalEnv)
					
					##calculate pairwise distance among the EICs				
					if(isDistParallel)
					{
						r<-parDistCal4(cl,isUnion=F)##parallel verison
					}else
					{
						r<-DistCal4(isUnion=F)#non-parallel version
					}	
					
					##convert distance matrix to triangle matrix
					distance <- as.dist(r)
					maxIntraDist<- 15 #cutoff to decide if split
					
					if(clustingType=="h")
					{
						##################################
						#hierarchical clustering and cut 
						#into groups by the height
						##################################
						clustResut<-hclust(distance)
						FinalClustResut<-cutree(clustResut,h=maxIntraDist)
						Clusters<-unique(FinalClustResut)
						if(length(Clusters)>1)isMultGroups<-T
						
					}else
					{	###############################################################
						#partitional clustering and split into groups by the silehoute
						###############################################################
						if(max(distance)<=maxIntraDist)
						{###if the max distance is already small enought, no need to deconvolute
							isMultGroups<-FALSE
						}else
						{
							
							if(nPks>2)
							{
								########################################
								#try different k to decide the cluster number by intra dist
								#######################################
								
								maxClustNumber<-min(20,nPks-1)
								asw<-numeric(maxClustNumber)
								for (k in min(maxClustNumber,max(2,minCmp)):maxClustNumber)
								{
									clustResut <-pam(distance,k)
									Clusters<-unique(clustResut$clustering)								
									isSplit<-F
									
									for(curCluster in Clusters)
									{
										#####get masses belong to current cluster
										curGroup<-names(clustResut$clustering[clustResut$clustering==curCluster])
										#######does not average the single-ion cluster which has zero intra dist
										if(length(curGroup)>1)
										{#######calculate the average intra distance to see if the current group is a valid single group of fragments
											
											if(max(as.dist(r[curGroup,curGroup]))>maxIntraDist)
											{##if current group is not valid, skip to the next K,
												#and the current k value 
												isSplit<-T
												break
											}
											
										}
									}
									if(!isSplit)
										asw[k]<-clustResut$silinfo $ avg.width
								}
								k.best <- which.max(asw)
								if(k.best>1)
								{
									FinalClustResut<-pam(distance,k.best,cluster.only=TRUE,stand=TRUE)
									
									#################################################
									#calculate the average intra dist in each cluster
									###############################################
									isTooBigIntraDist<-NULL
									avgIntraDistVec<-NULL
									Clusters<-unique(FinalClustResut)
									for(curCluster in Clusters)
									{
										#####get masses belong to current cluster
										curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
										
										IntraDist<-mean(as.dist(r[curGroup,curGroup]))
										isTooBigIntraDist<-c(isTooBigIntraDist,ifelse(IntraDist>6,T,F))
										if(length(curGroup)>1)#######does not average the single-ion cluster which has zero intra dist
											avgIntraDistVec<-c(avgIntraDistVec,IntraDist)
										
									}
									#################################################################
									#if after splitting,intra dist is smaller than before splitting
									#then accept the best clustering schema, 
									#and earch for unique mass within each cluster
									###############################################################
									if(mean(avgIntraDistVec)<mean(distance))
										isMultGroups<-TRUE
								}
							}else#if ==2, then directly set flag=T
							{
								isMultGroups<-TRUE
							}
						}
					}
				}
				
				
				
				
				######################
				#Model peak selection 
				######################
				
				modelPkList<-NULL
				
				# single group
				if (!isMultGroups && nPks >=1)
				{
					if (nPks==1) {
						modelPkList<-goodShapePeaklist
					} else 
					{
						##choose best one as model peak
						mdlPkCandidates<-goodShapePeaklist
						
						##scoring system
						#mass score
						mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
						
						##gassian similarity score
						mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
						
						##peak height score
						mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))					
						mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
						mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
						mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
						
						##linear combination of 3 scores
						scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
						#scores<-(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3+0.4*mdlPkCandidates$f4)
						
						##rescale score to (0,1000)
						mdlPkCandidates$score<-scale(scores,F,0.001)
						
						##select the max score (>500)
						mdlPkCandidates <- subset(mdlPkCandidates,score>=500)
						if (nrow(mdlPkCandidates)>0)
						{
							modelPkList<- mdlPkCandidates[which.max(mdlPkCandidates$score),]
						}
					}
				}
				
				
				if(isMultGroups)
				{
					if(nPks==2)
					{
						##only two candidates, simply assign them as model peaks
						#modelPkVec<-rownames(goodShapePeaklist)
						modelPkList<-goodShapePeaklist
					}else
					{
						nMinFragment<-1
						if(length(Clusters)>0)
						{	
							##select model peak from each group[!isTooBigIntraDist]
							for(curCluster in Clusters)
							{
								##get masses belong to current cluster
								curGroup<-names(FinalClustResut[FinalClustResut==curCluster])
								
								##get model peak candidates 
								mdlPkCandidates<-goodShapePeaklist[curGroup,]
								
								if(nrow(mdlPkCandidates)>0)
								{
									##mass score
									mdlPkCandidates$f1<-scale(mdlPkCandidates$mz,min(mdlPkCandidates$mz),diff(range(mdlPkCandidates$mz)))	
									##gassian similarity score
									mdlPkCandidates$f2<-scale(cos(mdlPkCandidates$gss*pi/180),min(cos(mdlPkCandidates$gss*pi/180)),diff(range(cos(mdlPkCandidates$gss*pi/180))))
									##peak height score
									mdlPkCandidates$f3<-scale(log(mdlPkCandidates$Intensity),min(log(mdlPkCandidates$Intensity)),diff(range(log(mdlPkCandidates$Intensity))))								
									
									##assign zero value when each factor is NA
									mdlPkCandidates$f1[is.na(mdlPkCandidates$f1)]<-0
									mdlPkCandidates$f2[is.na(mdlPkCandidates$f2)]<-0
									mdlPkCandidates$f3[is.na(mdlPkCandidates$f3)]<-0
									
									##linear combination of 3 scores
									scores<-(10/7)*(0.2*mdlPkCandidates$f1+0.4*mdlPkCandidates$f2+0.1*mdlPkCandidates$f3)
									
									##rescale score to (0,1000)
									mdlPkCandidates$score<-scale(scores,F,0.001)
									
									##select the max score
									currentMdlPk<-mdlPkCandidates[which.max(mdlPkCandidates$score),]
									modelPkList<-rbind(modelPkList,currentMdlPk)
								}
							}
						}
					}			
				}						
			}
			
			
			###################################
			#decomposition based on model peaks
			###################################
			
			decom = 1
			while (decom==1) 
			{
				modelPkDistance <- NULL
				specList<-NULL
				componentList<-NULL
				mdlPkIDs <- NULL
				modelPkList3 <- NULL
				modelPkList2 <- NULL
				modelPkList1 <- NULL
				##No co-eluting compounds (only one or no model peak detected)
				##Then extract spectrum directly
				if(nrow(modelPkList)<=1||is.null(modelPkList))
				{
					##add model peak to the spectrum
					specList[[1]]<-data.frame(mz=as.integer(),int=as.integer())
					
					##No model peak detected,then simply extract 
					##spectrum from the TIC peak apex time
					if(nrow(modelPkList)==0||is.null(modelPkList))
					{	
						##get 1st mass				
						curMz<-as.character(localEICPeakList$mz[1])
						
						##get TIC peak apex time and boundaries
						maxInd<-curTICPk$pkInd
						curLbound<-curTICPk$lboundInd
						curRbound<-curTICPk$rboundInd
						
						##Init component list
						componentList[[1]]<-data.frame(mz=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),pkInd=as.integer(),Intensity=as.integer(),
								area=as.integer(),windowID=as.integer(),compoundID=as.integer(),isModel=as.integer(),isUnique=as.integer(),factor=as.numeric())
						
						##name spec list with time and mz
						clustID<-format((maxInd-1)*ScanInterval+delaytime,digit=4)
						names(specList)[1]<-paste(clustID,curMz)
						
						##add rest peaks from the local window 
						##to component list and spectral list
						for(i in 1:nrow(localEICPeakList))
						{
							curPkID<-rownames(localEICPeakList)[i]					
							curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
							curProfile1$int[is.na(curProfile1$int)]<-0
							curArea<-sum(curProfile1$int)
							
							curMz<-as.character(localEICPeakList$mz[i])
							curInt<-subset(curProfile1,ET==maxInd)$int[1]
							
							if(!is.na(curInt)&&curInt>0)
							{
								specList[[1]]<-rbind(specList[[1]],data.frame(mz=as.integer(curMz),int=curInt))#add umass to the spectrum					
								componentList[[1]]<-rbind(componentList[[1]],data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,
												pkInd=maxInd,Intensity=curInt,area=curArea,windowID=windowID,compoundID=1,isModel=0,isUnique=0,factor=0))
							}
						}
						
						##Mark the ion with maximal intensity as model
						componentList[[1]][which.max(componentList[[1]]$Intensity),]$isModel=1
						modelPkList1 <- curTICPk$pkInd
						
					}else
					{
						##single one model peak detected,then extract 
						##spectrum from the model peak apex time
						##Use model peak apex,boundary index
						curMz<-modelPkList$mz[1]
						maxInd<-modelPkList$pkInd[1]
						curLbound<-modelPkList$lboundInd[1]
						curRbound<-modelPkList$rboundInd[1]
						curPhHeight<-modelPkList$Intensity[1]	
						curPkID<-rownames(modelPkList)[1]
						
						##model peaks could not exist in all profiles which save all peaks with isApex=1
						##sometimes a boundary peak of neighboring merged peak is inclued because of its 
						##peak index is inside the decon range
						curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),profiles[[curPkID]],by.y="ET",all.x=T)
						curProfile1$int[is.na(curProfile1$int)]<-0
						curArea<-sum(curProfile1$int)
						componentList[[1]]<-data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,pkInd=maxInd,
								Intensity=curPhHeight,area=curArea,windowID=windowID,compoundID=1,isModel=1,isUnique=0,factor=0)
						
						#create spec list as NIST format for future identifiation
						clustID<-format((maxInd-1)*ScanInterval+delaytime,digit=4)
						names(specList)[1]<-paste(clustID,curMz)						
						
						##peak list id to extract spectrum
						peakListIndex <-c(1:nrow(localEICPeakList))				
						if (curPkID%in%rownames(localEICPeakList)) {
							##remove already detected model peaks
							peakListIndex <- peakListIndex[-which(rownames(localEICPeakList)==curPkID)]
						}
						
						##add rest peaks from the local window 
						##to component list and spectral list
						for(i in peakListIndex)
						{
							curPkID<-rownames(localEICPeakList)[i]					
							curProfile1<-merge(data.frame(ET=as.integer(curLbound:curRbound)),allprofiles[[curPkID]],by.y="ET",all.x=T)
							curProfile1$int[is.na(curProfile1$int)]<-0
							curArea<-sum(curProfile1$int)						
							curMz<-as.character(localEICPeakList$mz[i])
							curInt<-subset(curProfile1,ET==maxInd)$int[1]
							
							if(!is.na(curInt)&&curInt>0)
							{
								specList[[1]]<-rbind(specList[[1]],data.frame(mz=as.integer(curMz),int=curInt))#add umass to the spectrum					
								componentList[[1]]<-rbind(componentList[[1]],data.frame(mz=as.integer(curMz),lboundInd=curLbound,rboundInd=curRbound,
												pkInd=maxInd,Intensity=curInt,area=curArea,windowID=windowID,compoundID=1,isModel=0,isUnique=0,factor=0))
							}
						}						
					}
					
					componentList[[1]][which.max(componentList[[1]]$Intensity),]$isUnique=1
					##replace decom as 0, stop while loop
					decom=0
					modelPkList2 <- modelPkList$pkInd
				}else
				{
					##More than one model peaks
					##Decompostion
					
					#######################################
					#Makeup 1 for splitting issue
					#Checking pairwise model peaks distance
					#######################################
					
					modelPkDistance <- as.matrix(modelPkDist(modelPkList,profiles))
					rmmodelPk <- NULL
					ModelPkNames <- rownames(modelPkDistance)
					
					##model same matrix as modelPKDistance with index
					psmatrix <- matrix(c(1:length(modelPkDistance)),ncol=nrow(modelPkDistance),byrow=T)
					psmodel <- which(modelPkDistance<=4&modelPkDistance>0)
					
					##get index of model peaks with smaller distance
					for (i in psmodel) {
						pos <- matrix.index(psmatrix, i)
						rmmodelPk <- c(rmmodelPk,ModelPkNames[pos[1]],ModelPkNames[pos[2]])
					}
					rmmodelPks <- unique(rmmodelPk)
					
					rmmodelPkList <- modelPkList[rmmodelPks,]
					modelPkList_temp <- rmmodelPkList[which.max(rmmodelPkList$score),]
					modelPkList <- modelPkList[!rownames(modelPkList)%in%rmmodelPks,]
					
					##Only keep the one with highest score
					modelPkList <- rbind(modelPkList,modelPkList_temp)				
					mdlPkIDs<-rownames(modelPkList)
					
					##only one model peak after testing
					##GO back to the upper step
					if (nrow(modelPkList)==1) decom=1
					
					else {
						##still multiple models
						lbound<-min(modelPkList$lboundInd)
						rbound<-max(modelPkList$rboundInd)
						
						#################################
						#get all EIC profiles within the window 
						#X matrix format:column-EIC,row-scan
						#############################################
						
						X<-NULL
						mzVec<-unique(localEICPeakList$mz)
						nMz<-length(mzVec)
						for(i in 1:nMz)
						{
							mz<-mzVec[i]
							mzInd<-which(vectorMz==mz)
							startInd<-(mzInd-1)*totalscan+1
							endInd<-startInd+totalscan-1
							curVecInt <- vecInt[startInd:endInd]
							X<-cbind(X,curVecInt[lbound:rbound])
						}
						colnames(X)<-mzVec	
						nScans<-nrow(X)
						
						###########################################
						#Get model peak profiles as S matrix
						#each column is the EIC for one model peak
						###########################################
						
						S<-NULL
						modelPkList$area<-0
						for(i in 1:nrow(modelPkList))
						{
							curProfile1<-merge(data.frame(ET=as.integer(lbound:rbound)),profiles[[rownames(modelPkList[i,])]],by.y="ET",all.x=T)
							
							##expand the profile by filling zero values to those uncovered area
							curProfile1$int[is.na(curProfile1$int)]<-0
							S<-cbind(S,curProfile1$int)
							modelPkList[i,]$area<-sum(curProfile1$int)		
							
							##init the component & spec list	
							componentList[[i]]<-data.frame(mz=as.integer(),lboundInd=as.integer(),rboundInd=as.integer(),pkInd=as.integer(),Intensity=as.integer(),
									area=as.integer(),windowID=as.integer(),compoundID=as.integer(),isModel=as.integer(),isUnique=as.integer(),factor=as.numeric())
							##add model peak to the spectrum
							specList[[i]]<-data.frame(mz=as.integer(),int=as.integer())
							curRT<-format((modelPkList$pkInd[i]-1)*ScanInterval+delaytime,digit=4)
							names(specList)[i]<-paste(curRT,modelPkList$mz[i])
						}					
						colnames(S)<-mdlPkIDs
						
						##Residual minimization
						for(i in 1:ncol(X))
						{
							curMz<-colnames(X)[i]						
							M<-X[,i]#mixture signal
							
							srcIDs<-mdlPkIDs
							A<-optim(par=rep(0,ncol(S)),fn=fr, gr = NULL,M,S,lower = 0, method="L-BFGS-B")$par
							A<-as.matrix(A)
							
							##restore each ion signals and save them into 
							##component list ans speclist
							for(m in 1:length(srcIDs))
							{
								curmdlPkID<-srcIDs[m]
								j<-which(mdlPkIDs==curmdlPkID)
								
								if(A[m]>0)
								{
									specList[[j]]<-rbind(specList[[j]],data.frame(mz=as.integer(curMz),int=A[m,]))###add the curent restored mixing factor to spec
									if (modelPkList[curmdlPkID,]$mz==curMz) {
										componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,
														rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=as.integer(modelPkList[curmdlPkID,]$Intensity*A[m,]),
														area=as.integer(modelPkList[curmdlPkID,]$area*A[m,]),windowID=windowID,compoundID=j,isModel=1,isUnique=0,factor=A[m]))
									}
									
									else {
										componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,
														rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=as.integer(modelPkList[curmdlPkID,]$Intensity*A[m,]),
														area=as.integer(modelPkList[curmdlPkID,]$area*A[m,]),windowID=windowID,compoundID=j,isModel=0,isUnique=0,factor=A[m]))
									}
								}							
							}
						}
						
						###############################################
						#Makeup 2 for splitting issue
						#Checking pairwise resovlved spectra similarity
						###############################################
						
						RepeatModelPk <- NULL
						Num <- nrow(modelPkList)
						indexPairs<-combinations(Num,2)
						for(i in 1:nrow(indexPairs))
						{
							index <- indexPairs[i,]
							index1 <- index[1]
							index2 <- index[2]
							
							##calculate time difference between two model peaks
							modelPKET <- modelPkList$pkInd[c(index1,index2)]
							
							score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
							if (score >=970||(abs(modelPKET[1]-modelPKET[2]) <=6&&score >=900))
							{
								selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index2],mdlPkIDs[index1])
								RepeatModelPk <- c(RepeatModelPk,selectedModelPkID)
							}
						}
						
						##select unique filtered model peaks 
						RepeatModelPk <- unique(RepeatModelPk)
						modelPkList <- modelPkList[!rownames(modelPkList)%in%RepeatModelPk,]
						modelPkList3 <- modelPkList$pkInd
						##The modelPk number keep the same, do not RE-decompose
						if (nrow(modelPkList) == Num) {
							decom = 0
							
							if (nrow(modelPkList) == 2) {
								mzlist1 <- subset(componentList[[1]],Intensity>=1000)$mz
								mzlist2 <- subset(componentList[[2]],Intensity>=1000)$mz
								SharedMS1 <- unique(c(mzlist1[mzlist1%in%mzlist2],CurNonUmassVec))
								SharedMS2 <- unique(c(mzlist2[mzlist2%in%mzlist1],CurNonUmassVec))
								
								unique_comp1 <- componentList[[1]][!componentList[[1]]$mz%in%SharedMS1,]
								unique_comp2 <- componentList[[2]][!componentList[[2]]$mz%in%SharedMS2,]
								
								componentList[[1]][componentList[[1]]$mz==unique_comp1[which.max(unique_comp1$Intensity),]$mz,]$isUnique=1
								componentList[[2]][componentList[[2]]$mz==unique_comp2[which.max(unique_comp2$Intensity),]$mz,]$isUnique=1
								
							}else
							{
								for (i in c(1:nrow(modelPkList)))
								{
									if (i==1||i==nrow(modelPkList))
									{
										if (i==1) {
											mzlist1 <- subset(componentList[[1]],Intensity>=1000)$mz
											mzlist2 <- subset(componentList[[2]],Intensity>=1000)$mz
											NonUniqueMS1 <- unique(c(mzlist1[mzlist1%in%mzlist2],CurNonUmassVec))						
											unique_comp1 <- componentList[[1]][!componentList[[1]]$mz%in%NonUniqueMS1,]				
											componentList[[1]][componentList[[1]]$mz==unique_comp1[which.max(unique_comp1$Intensity),]$mz,]$isUnique=1
										}else{
											mzlist1 <- subset(componentList[[i]],Intensity>=1000)$mz
											mzlist2 <- subset(componentList[[i-1]],Intensity>=1000)$mz
											NonUniqueMS1 <- unique(c(mzlist1[mzlist1%in%mzlist2],CurNonUmassVec))						
											unique_comp1 <- componentList[[i]][!componentList[[i]]$mz%in%NonUniqueMS1,]				
											componentList[[i]][componentList[[i]]$mz==unique_comp1[which.max(unique_comp1$Intensity),]$mz,]$isUnique=1
										}
									}else{
										mzlist1 <- subset(componentList[[i-1]],Intensity>=1000)$mz
										mzlist2 <- subset(componentList[[i]],Intensity>=1000)$mz
										mzlist3 <- subset(componentList[[i+1]],Intensity>=1000)$mz
										
										NonUniqueMS2 <- unique(c(mzlist2[mzlist2%in%mzlist1],mzlist2[mzlist2%in%mzlist3],CurNonUmassVec))					
										unique_comp2 <- componentList[[i]][!componentList[[i]]$mz%in%NonUniqueMS2,]				
										componentList[[i]][componentList[[i]]$mz==unique_comp2[which.max(unique_comp2$Intensity),]$mz,]$isUnique=1
									}						
								}
							}		
						}			
					} # else model PK num is >2 after testing distance	
				} # model PK num is >2 at first
			} # while statement
			
			##flag analyzed EIC peaks 
			EICpeaklist[rownames(localEICPeakList),]$flag<-1
			componentResults<-c(componentResults,componentList)
			specResults<-c(specResults,specList)
			modelPkList3s <- c(modelPkList3s,modelPkList3)
			modelPkList2s <- c(modelPkList2s,modelPkList2)
			modelPkList1s <- c(modelPkList1s,modelPkList1)
			
			##model peak id
			mdlPkIDVec<-c(mdlPkIDVec,mdlPkIDs)		
		}	
	}
	
	componentResults<-do.call(rbind,componentResults)	
	cat(fileName,"finished deconvolution...\n")
	
	##################################
	#output lib matching results
	###################################
	
	##similarity cutoff=600
	libmatch <- libMatching(lib,specResults,600)
	write.csv(libmatch,file=paste(params$WorkDir,"/output/decon/",fileName,"_libmatch.csv",sep=""))
	
	write.csv(modelPkList3s,file=paste(params$WorkDir,"/output/decon/",fileName,"_modelPkList3s.csv",sep=""))
	write.csv(modelPkList2s,file=paste(params$WorkDir,"/output/decon/",fileName,"_modelPkList2s.csv",sep=""))
	write.csv(modelPkList1s,file=paste(params$WorkDir,"/output/decon/",fileName,"_modelPkList1s.csv",sep=""))
	##############################################################
	#collect spectrum result in dataframe for reverse lib matching
	##############################################################
	
	##exp spectra as lib, get top 1 hit for each std
	evaluateIdentification(refSpec=lib,inputSpec=specResults,RT_Tolerance=70,minSpecSimilarity=650,fileName,withRT=T,params)
	cat(fileName,"finish matching...\n")
	write.csv(mdlPkIDVec,file=paste(params$WorkDir,"/output/decon/",fileName,"mdlpklist.csv",sep=""))
	
	################################################
	#collect component result in dataframe for alignment
	###################################################
	
	##filter out zero-int fragments
	componentResults<-subset(componentResults,Intensity>0,select=c("mz","lboundInd","rboundInd","pkInd","Intensity","area","windowID","compoundID","isModel","isUnique","factor"))
	write.csv(componentResults,row.names=FALSE,file=paste(params$WorkDir,"/output/decon/",fileName,"Decon.csv",sep=""))
	cat(fileName,"writing component results...\n")
}
