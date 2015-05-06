
## to select either Qmass or Summation of all fragments
## to represent a component/potential compound

data_select <- function(params) 
{
	choice<-params$QuantitationType
	paraDB<- params$dbFile
	JobName<- params$JobName
	workDir <- params$WorkDir
	m <- dbDriver("SQLite")
	con <- dbConnect(m, dbname = paraDB)
	tname<-"bioMarkerList"

	if(!dbExistsTable(con,tname))
		stop(paste("Table '",tname,"' does not exist!",sep=""))
	trainingSet<-dbReadTable(con,tname)
	Ncol <- ncol(trainingSet)
	
	if (choice=="S") {
		
		##to summarize all ions of each component
		##Keep two columns:clustID and alignedET 
		trainingSet<-trainingSet[,c(1,2,8:Ncol)]
		Ncol <- ncol(trainingSet)
		trainingSet <- aggregate(trainingSet[,3:Ncol],by=list(clustID=trainingSet$clustID,alignedET=trainingSet$alignedET),FUN=sum)
		
		## to select Qmass occuring more than nMarjority of samples
		nFile<-length(params$DataFiles)
		nMajority<-round(as.numeric(params$sampleRatio)*nFile)
		
		
		CountNum <- NULL		
		for (i in 1:nrow(trainingSet)) {
			count <- sum(trainingSet[i,][,c(3:Ncol)]>0)
			CountNum <- c(CountNum,count)
		}
		
		##add count number to filter and then remove it
		trainingSet$"count"<- CountNum
		trainingSet <- subset(trainingSet,count>=nMajority)
		trainingSet <- trainingSet[,-which(colnames(trainingSet)=="count")]
	}else {
		
		##Other three conditions: keep three columns-clustID,alignedET and MZ	
		if (choice=="Q") {
			trainingSet<-subset(trainingSet,isMaxUmass==1)
			trainingSet<-trainingSet[,c(1:3,8:Ncol)]
		}
		
		##To select unique ions
		if (choice=="U") {
			trainingSet<-subset(trainingSet,isUnique==1)
			trainingSet<-trainingSet[,c(1:3,8:Ncol)]
		}
		
		##to select model peak for each compoenent
		if (choice=="M") {
			trainingSet<-subset(trainingSet,isModel==1)
			trainingSet<-trainingSet[,c(1:3,8:Ncol)]
		}
		
		## to select Qmass occuring more than nMarjority of samples
		nFile<-length(params$DataFiles)
		nMajority<-round(as.numeric(params$sampleRatio)*nFile)
		
		Ncol <- ncol(trainingSet)
		CountNum <- NULL
		
		for (i in 1:nrow(trainingSet)) {
			count <- sum(trainingSet[i,][,c(4:Ncol)]>0)
			CountNum <- c(CountNum,count)
		}
		# add count number to filter and then remove it
		trainingSet$"count"<- CountNum
		trainingSet <- subset(trainingSet,count>=nMajority)
		trainingSet <- trainingSet[,-which(colnames(trainingSet)=="count")]
	}
	print("selection is done!\n")
	# clean up
	dbDisconnect(con)
	# save the csv file under output dir.
	write.csv(trainingSet,paste(workDir,"output/",JobName,"_",choice,"_","BiomarkerTable.csv",sep=""),row.names=F)
}	
