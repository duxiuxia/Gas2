# correct missing value in biomarkerlist
# 
# Author: yni
###############################################################################

#dbFile<-"../testing/yni/yni.db"
readCmpChkList<-function(dbFile)
{
	library(RSQLite)
	m <- dbDriver("SQLite")
	con <- dbConnect(m, dbname = dbFile)
	# query, MissingCmpList generates the zero/missing values
	# three columns: clustID, Qmass, FileID (#)
	rs <- dbGetQuery(con, "select distinct * from MissingCmpList")
	# clean up
	dbDisconnect(con)
	rs	
}


readMatchResult<-function(dbFile)
{
	library(RSQLite)
	m <- dbDriver("SQLite")
	con <- dbConnect(m, dbname = dbFile)
	# query, read data from libSearchResult
	rs <- dbGetQuery(con, "select * from LibSearchResult")
	# clean up
	dbDisconnect(con)	
	rs	
}

CompoundRedetection<-function(params,scoreTolerance=750)
{
	dbFile<-params$dbFile
	##########################
	#read compound check list
	##########################
	
	####  Load Missing data list
	CmpList<-readCmpChkList(dbFile) 

	if(nrow(CmpList)>0)
	{
		#### load match result
		matchResult<-readMatchResult(dbFile) # lib match result
		#####load library
		lib<-NULL
		lib<-loadSpecSQL("Spec_QcCompound","../JavaGUI/lib.db")
		
		#####load RT Data
		fileInfo<-read.csv(paste(params$WorkDir,"/output/",params$JobName,"_RT.csv",sep=""))
		params$delaytime<-fileInfo$firstRT[1]
		params$ScanInterval<-fileInfo$ScanInterval[1]
		
		#####get missing data fileID
		FileIDList<-unique(CmpList$FileID)	
		
		##################################
		# Re-check out the Qmass Intensity
		##################################
		
		for(i in 1:length(FileIDList))
		{
			###### File list Loop
			curFileID<-FileIDList[[i]]	
			# fileID should keep the same with curFileID
			fileName<-as.character(subset(fileInfo,FileID==curFileID)$FileName)					
			
			###read decon result using fileName	
			curCompoundList<-subset(CmpList,FileID==curFileID)
			curDeconData <- read.csv(paste(params$WorkDir,"output/decon/",fileName,"Decon.csv",sep=""))
			# RT unit: s
			curDeconData$RT <- ((curDeconData$pkInd-1)*params$ScanInterval+params$delaytime)*60
			
			# initialization
			isUpdate<-F
			
			### in case there exists multiple missing data of each file
			for(j in 1:nrow(curCompoundList)) 
			{
				###### Qmass Loop for each file
				curQmass <- curCompoundList$Qmass[j]
				curclustID <- curCompoundList$clustID[j]
				# get the comp name of Qmass
				##there multiple match result for each component depending on the top n parameters set on Java GUI
				##but here only the top match is considered as the correct identification
				matchedCmps<-subset(matchResult,clustID==curclustID)				
				curCompoundName<-matchedCmps[order(matchedCmps$score,decreasing=T),]$compoundName[1]
				## extracted the lib spectra for lib matching
				curLibSpec<-lib[[curCompoundName]]
				
				#plot(curLibSpec$`Dodecanoic acid`,type="h")
				
				######################################
				# get component spectra from decon result
				######################################
	
				curRT <- subset(matchResult,clustID==curclustID)$RT
				### time interval is 4s
				CandidateCompListFull <- subset(curDeconData,RT> (curRT-2)& RT <(curRT+2))
				###  compoundID <0 means fragments filtered out in pkpicking and decon step
				### -1: low intensity; 
				### -2: not enough # of fragments for a window
				### -3: not enough # for a cluster
				### -4: failed to satisfy s/n cutoff in pkpicking
				
				CandidateCompListFilter <- subset(CandidateCompListFull,compoundID>0)
				# Default maxScore is 0, so if FilterList is NULL, that keep the zero value;
				# Or maxScore will be calculated & compared.
				maxScore <- 0
				if (nrow(CandidateCompListFilter) > 0)
				{
					# create a new ID to represent comp
					CandidateCompListFilter$ID <- paste(CandidateCompListFilter$windowID,CandidateCompListFilter$compoundID)					
					cmpIDs<-unique(CandidateCompListFilter$ID)
					
					### Match curLibSpec with CurSpect List
					### Output match score list, find maximal score index and coresponding spectra list (int & mz)					
					scoreList <- NULL
					for (m in 1:length(cmpIDs))
					{
						##### spectra matching comparison Loop
						curSpec<-subset(CandidateCompListFilter,ID==cmpIDs[m],select=c("mz","Intensity"))
						names(curSpec)[2]<-"int"
						curScore <- specDistCal(curLibSpec,curSpec,isWeight=TRUE,isNist=TRUE)
						scoreList <- c(scoreList,curScore)
					}
					maxScore <- max(scoreList) 
				}
				
				##if the compound is detected in current file
				if(maxScore >=scoreTolerance)
				{
					maxCompID<-cmpIDs[which.max(scoreList)]					
					### Find coresponding intensity of Qmass
					matchedComponent<-subset(CandidateCompListFilter,ID==maxCompID)
					#there may be more than one peak of the same mz for one component, which is a bug in the current version ADAP1.0
					QmassInt<-subset(matchedComponent,mz==curQmass)$Intensity[1]
					### get the mean RT of component mass, compared with qmass RT list later
					MeanCurCompRT <- mean(matchedComponent$RT)
					
					if(is.na(QmassInt))
					{				
						#if qmass does not exist in the matched component
						#it is in decon result with -1,-2,-3 value, then directly recover it 
						#by looking for the nearest peak (compared to the mean RT of matchedComponent)
				
						curWindID <- unique(matchedComponent$windowID)
						curMzFilteredList <- subset(CandidateCompListFull,windowID==curWindID & compoundID < 0 & mz==curQmass)
						
						if(nrow(curMzFilteredList)>0)
						{
							QmassPos <- which.min(abs(curMzFilteredList$RT-MeanCurCompRT))
							QmassInt <- curMzFilteredList$Intensity[QmassPos]
							isUpdate=T
							
						}else
						{
							#if Qmass is still missing,then skip
							#then relax the window to search peak list for the Qmass peak
													
							curMzPks<-subset(CandidateCompListFull,mz==curQmass)
							if(nrow(curMzPks)>0)
							{	#if qmass is found in relaxed window, add it to spec
								QmassPos<-which.min(abs(curMzPks$RT-MeanCurCompRT))
								QmassInt<-curMzPks[QmassPos,]$Intensity
								curSpec<-rbind(curSpec,data.frame(mz=as.integer(curQmass),int=QmassInt))
								#rematch the lib
								curScore <- specDistCal(curLibSpec,curSpec,isWeight=TRUE,isNist=TRUE)
								#if score is improved after new found qmass, then set update as True
								if(curScore>=max(scoreList))
								{	
										isUpdate=T
								}else
								{
									cat("Qmass detected, but does not improve score!")
										
								}
							}else
							{
								cat("Qmass is not detected in this compound!")
							}
								
						}
					}else
					{
						isUpdate=T
					}
					
				}
				else  
				{
					#if not component is matched with a high score
					#then go back to peak list to gather all the fragment peaks
					newSpec<-NULL
					for(curMz in curLibSpec$mz)
					{
						curMzPks<-subset(CandidateCompListFull,mz==curMz)
						if(nrow(curMzPks)>0)
						{
							nPkpos<-which.min(abs(curMzPks$RT-curRT))
							newSpec<-rbind(newSpec,curMzPks[nPkpos,])
						}
						
					}
					newSpec<-subset(newSpec,select=c("mz","Intensity"))
					names(newSpec)[2]<-"int"
					curScore <- specDistCal(curLibSpec,newSpec,isWeight=TRUE,isNist=TRUE)
					#match the detected component to check if it is confident
					if(curScore>=scoreTolerance)
					{
						#there may be more than one peak of the same mz for one component, 
			            #which is a bug in the current version ADAP1.0
						QmassInt<-subset(newSpec,mz==curQmass)$int[1]
						if(!is.na(QmassInt))
						{
							isUpdate=T
						}else
						{
							cat("Mass",curQmass,"is not observed in file (yet component detected with matching score higher than tolerance)",fileName,"\n")
							
						}
						
					}
				}
				
				if(isUpdate)
				{
					###update biomarker table in SQL (Qmass Intensity)
					m <- dbDriver("SQLite")
					con <- dbConnect(m, dbname = dbFile)
					
					dbBeginTransaction(con)
					res <- dbSendQuery(con, paste("update bioMarkerList set X",curFileID,"=",QmassInt," where mz=",curQmass," and clustID=",curclustID,sep=""))
					dbCommit(con)
					# clean up
					dbDisconnect(con)	
				}
				
				
				
				##############################
				#if alignment is successful,only Qmass is missing, 
				#then get the Qmass peak back from decon or peakpicking result
				##########################################################
				
				
				
				
				##########################
				#if the whole component is missing
				#then perform target analysis
				##########################
	#			missingComponent<-targetAnalysis(cmpSpec,RTwindow,CurFile)
				
			
				
			}
		}
	}else
	{
		cat("no compound is selected!\n")
	}
}


