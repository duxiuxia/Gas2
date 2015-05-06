
##############################################
# write each model pk prfile into standard format
# ready for outputing text file
# Nov. 22, 2011
##############################################

getModelPkProfile_old<-function(curModel,curModelInt,windowID,lbound,rbound,model,compNum)
{
	
	strBlock<-sprintf("peakInd:%s",curModel$pkInd)
	strBlock<-paste(strBlock,sprintf("model mass:%d",curModel$mz),sep="\n")
	strBlock<-paste(strBlock,sprintf("model gss:%.2f",curModel$gss),sep="\n")
	strBlock<-paste(strBlock,sprintf("windowID:%d",windowID),sep="\n")
	strBlock<-paste(strBlock,sprintf("compoundID:%d",model),sep="\n")
	strBlock<-paste(strBlock,sprintf("window lboundInd:%d",lbound),sep="\n")
	strBlock<-paste(strBlock,sprintf("window rboundInd:%d",rbound),sep="\n")
	
	ModelIntBlock<-NULL
	for(j in 1:length(curModelInt))
	{
		strPair<-paste(curModelInt[j]," ",sep="")
		ModelIntBlock<-paste(ModelIntBlock,strPair)
	}
	strBlock<-paste(strBlock,ModelIntBlock,sep="\n")
	strBlock
}	

getModelPkProfile<-function(curModel,curModelInt,windowID,lbound,rbound,model,compNum)
{
	
	strBlock<-sprintf("peakInd:%s",curModel$pkInd)
	strBlock<-paste(strBlock,sprintf("model mass:%d",curModel$mz),sep="\n")
	strBlock<-paste(strBlock,sprintf("model gss:%.2f",curModel$gss),sep="\n")
	strBlock<-paste(strBlock,sprintf("windowID:%d",windowID),sep="\n")
	strBlock<-paste(strBlock,sprintf("compoundID:%d",model),sep="\n")
	strBlock<-paste(strBlock,sprintf("compID:%d",compNum+model-1),sep="\n")
	strBlock<-paste(strBlock,sprintf("window lboundInd:%d",lbound),sep="\n")
	strBlock<-paste(strBlock,sprintf("window rboundInd:%d",rbound),sep="\n")
	
	ModelIntBlock<-NULL
	for(j in 1:length(curModelInt))
	{
		strPair<-paste(curModelInt[j]," ",sep="")
		ModelIntBlock<-paste(ModelIntBlock,strPair)
	}
	strBlock<-paste(strBlock,ModelIntBlock,sep="\n")
	strBlock
}


localMin<-function(curVecInt,localwindowsize)
{
	
	minInd<-NULL
	for(i in 1:length(curVecInt))
	{	
		j<-i+localwindowsize-1
		localwindow<-curVecInt[i:j]
		localMinInd<-which.min(localwindow)
		minInd<-c(minInd,localMinInd+i-1)
		i<-j+1
		
	}
}

subProfile<-function(curFeatureID,Profiles)
{

	curProfs<-Profiles[[curFeatureID]]
	colnames(curProfs)[2]<-curFeatureID
	curProfs
	
}

## Checking span is even or odd number
## modified Sep 8, 2011, yan
peaks<-function (x, span = 3) 
{
	z <- embed(as.vector(x), span)
	s <- span%/%2
	
	if (span%%2==1) {
		result <- max.col(z,ties.method="first") == 1 + s
		c(rep(FALSE, s), result, rep(FALSE, s))
	}else if (span%%2==0) {
		result <- max.col(z,ties.method="first") == s
		c(rep(FALSE, s), result, rep(FALSE, s-1))
	}
}


## gaussian curve fitting
GaussianDisimilarity<-function(yy) {
	res<-90
	xx<-1:length(yy)
	tryCatch({fit1<-nls(yy ~ (w/sqrt(2*pi*sigma1^2)) * exp(-(xx - mu)^2/(2 * sigma1^2)),control = list(maxiter = 50),  
						start = list(w = max(yy)*sqrt(2*pi*1^2), sigma1=1,mu = xx[which.max(yy)]),trace = F)
				res<-dotProduct(fitted(fit1),yy)},
			error=function(er)
			{
				
				#cat(mz,er);
			}
	)
			#ifelse(is.na(res),90,res
		res
}

## Produdce two mirror image: left/right half
## The cutting point is apex
GaussianDisimilarity1<-function(yy,apex) {
	
	maxInd<-apex
	endInd<-length(yy)
	yy1<-yy[c(1:(maxInd-1),maxInd:1)]
	yy2<-yy[c(endInd:(maxInd+1),maxInd:endInd)]
	
	dist1<-GaussianDisimilarity(yy1)
	dist2<-GaussianDisimilarity(yy2)
	c(dist1,dist2)
	
}

## Find the apex and make a mirror image based on
## the gss calculation that indicates which half is 
## more closed to guassian peak

copyShape<-function(p,apex,from) {
	xx<-p$ET
	yy<-p$int
	maxInd<-which(xx==apex)
	endInd<-length(yy)	
	if(from==1) 
	{
		yy<-yy[c(1:(maxInd-1),maxInd:1)]
		xx<-xx[1]:(xx[1]+length(yy)-1)
		
	}else
	{
		yy<-yy[c(endInd:(maxInd+1),maxInd:endInd)]	
		
		# in case the first several peaks starting almost at the 1st scan
		# if the mirror image comes from the right part, the resulting profile
		# would have negative left boundary
		if ((xx[endInd]-length(yy)+1)<0) {
			yy<-yy[c(1:(maxInd-1),maxInd:1)]
			xx<-xx[1]:(xx[1]+length(yy)-1)
		}else {
			xx<-(xx[endInd]-length(yy)+1):xx[endInd]
		}	
	}	
	data.frame(ET=xx,int=yy)
}

Sharpness<-function(yy) {
	res<-0
	localPeakInd<-which.max(yy)
	
	nLength<-length(yy)
	shrpVec<-c((yy[2:localPeakInd]-yy[1:(localPeakInd-1)])/yy[1:(localPeakInd-1)],(yy[localPeakInd:(nLength-1)]-yy[(localPeakInd+1):nLength])/yy[(localPeakInd+1):nLength])
	res<-sum(shrpVec[is.finite(shrpVec)])
	res
}


Sharpness1<-function(yy,apex) {
	maxInd<-apex
	endInd<-length(yy)
	yy1<-yy[c(1:(maxInd-1),maxInd:1)]
	yy2<-yy[c(endInd:(maxInd+1),maxInd:endInd)]
	
	shrp1<-Sharpness(yy1)
	shrp2<-Sharpness(yy2)
	c(shrp1,shrp2)
}


modelPkDist <- function(modelPkList,profiles) {
	if (!is.null(nrow(modelPkList))) {
	modelPkprofiles<-profiles[as.character(rownames(modelPkList))]
	assign("Global.curProfs", value=modelPkprofiles, envir = .GlobalEnv)
	r<-DistCal4(isUnion=F)
	#distance <- as.dist(r)
	#distance
	r
	}
}

matrix.index <- function(a, value) {
	idx <- which(data.frame(a)==value)
	col.num <- ceiling(idx/nrow(a))
	row.num <- idx - (col.num-1) * nrow(a)
	return(c(row.num, col.num))
}


TICDenoising<-function(inFilePath,WorkDir,codeDir,isBslnCrt=F) 
{
	
	library("ncdf")
	source(paste(codeDir,"pipeline.r",sep="/"))
	###############
	#read EIC data
	###############
	#inFilePath <- DataFilelist[[fileindex]]
	fileName<-parseFileName(inFilePath)
	cat(fileName,"reading TIC data...\n")
	
	if(file.exists(inFilePath)==FALSE)next
	ncid <- open.ncdf(inFilePath,write=TRUE)
	vecInt <- get.var.ncdf(ncid, varid="total_intensity")
	vecET<-get.var.ncdf(ncid, varid="scan_acquisition_time")
	
	totalscan<-length(vecInt)
		
	###########
	#smoothing
	###########
	cat(fileName,"smoothing TIC data...\n")
	
	smoothingWindow<-10#50
	ma = rep(1, smoothingWindow)/smoothingWindow
	sn.ma=filter(vecInt,ma)
	sn.ma[which(is.na(sn.ma))]<-0
	vecInt<-sn.ma
	
	#####################
	#baseline correction
	#####################
	if(isBslnCrt==T)
	{
		nbinSize<-50#240
		isMin <- peaks(-vecInt, span=nbinSize)
		minInd<-which(isMin==TRUE)	
		
		intervalOfScans<-c(minInd[1],minInd[-1]-minInd[-length(minInd)])
		
		bgs<-c(rep(vecInt[minInd],intervalOfScans),rep(0,totalscan-minInd[length(minInd)]))
		f.lo <- loess(bgs ~ c(1:totalscan), span =0.05, degree = 2)
		bsln <- f.lo$fitted
		
		## baseline subtraction
		vecInt<-vecInt-bsln
		vecInt[vecInt<0]<-0
	}
	
	###################
	#write to new  CDF
	###################
	cat(fileName,"write denoised TIC data...\n")
	
	ncid1<-create.ncdf(filename=paste(WorkDir,"/output/TIC/denoised_",fileName,"_TIC.cdf",sep=""),vars=ncid$var)
	put.var.ncdf(ncid1,"total_intensity",vecInt)
	close.ncdf(ncid1)
	remove(ncid1)	
	
	close.ncdf(ncid)
	remove(ncid)	
}

Denoising<-function(inFilePath,WorkDir,codeDir,isBslnCrt=F,isSm=F,params) 
{
	library("ncdf")
	source(paste(codeDir,"pipeline.r",sep="/"))
	
	###############
	#read EIC data
	###############
	#inFilePath <- DataFilelist[[fileindex]]
	fileName<-parseFileName(inFilePath)
	EICFile<-paste(WorkDir,"/output/EIC/",fileName,"EIC.cdf",sep="")
	cat(fileName,"reading EIC data...\n")
	
	ncid <- open.ncdf(EICFile)
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	mzVec<-get.var.ncdf(ncid, varid="mzVec")
	totalscan<-length(vecInt)/length(mzVec)
	
	###########
	#smoothing
	###########
	cat(fileName,"smoothing EIC data...\n")
	
	DenoisedTIC <- 0
	
	time1<-Sys.time()
	for(mzInd in 1:length(mzVec))
	{
		mz<-mzVec[mzInd]	
		startInd<-(mzInd-1)*totalscan+1
		endInd<-startInd+totalscan-1
		curVecInt <- vecInt[startInd:endInd]
		
		if(isSm==T)
		{
			smoothingWindow<-as.integer(params$smoothing_span)# default is 20
			ma = rep(1, smoothingWindow)/smoothingWindow
			sn.ma=filter(curVecInt,ma)
			sn.ma[which(is.na(sn.ma))]<-0
			curVecInt<-sn.ma
		}	
		#####################
		#baseline correction
		#####################
		if(isBslnCrt==T)
		{
			nbinSize<-as.integer(params$baseline_span) #default is 240
			isMin <- peaks(-curVecInt, span=nbinSize)
			minInd<-which(isMin==TRUE)	
			
			if (length(minInd)==0) {bsln=0} # avoiding no minInd found - yan
			else {
				intervalOfScans<-c(minInd[1],minInd[-1]-minInd[-length(minInd)])
				
				bgs<-c(rep(curVecInt[minInd],intervalOfScans),rep(0,totalscan-minInd[length(minInd)]))
				f.lo <- loess(bgs ~ c(1:totalscan), span =0.05, degree = 2)
				bsln <- f.lo$fitted
			}
			
			## baseline subtraction
			curVecInt<-curVecInt-bsln
			curVecInt[curVecInt<0]<-0
		}
		
		## update vector of intensity
		vecInt[startInd:endInd]<-curVecInt
		DenoisedTIC <- DenoisedTIC + curVecInt
	}
	
	time2<-Sys.time()
	time2-time1

	######################
	#update EIC and TIC CDF
	#######################
	cat(fileName,"write denoised EIC data...\n")
	ncid1<-create.ncdf(filename=paste(WorkDir,"/output/EIC/denoised_",fileName,"EIC.cdf",sep=""),vars=ncid$var)
	put.var.ncdf(ncid1,"intVec",vecInt)
	put.var.ncdf(ncid1,"mzVec",mzVec)
	close.ncdf(ncid1)
	remove(ncid1)	
		
	cat(fileName,"write denoised TIC data...\n")
	dimTIC<-dim.def.ncdf(name="totalInt", units="", vals=1:length(DenoisedTIC), unlim=FALSE, create_dimvar=TRUE )
	varTIC<-var.def.ncdf("total_intensity","" ,dimTIC,0 )
	ncid1<-create.ncdf(filename=paste(WorkDir,"/output/TIC/denoised_",fileName,"_TIC.cdf",sep=""),vars=varTIC)
	put.var.ncdf(ncid1,"total_intensity",DenoisedTIC)
	close.ncdf(ncid1)
	remove(ncid1)	
			
	close.ncdf(ncid)
	remove(ncid)	
}


parTICDenoising<-function(params) 
{
	time1<-Sys.time()
	library("snow")
	###############
	#get parameters
	################
	WorkDir <- params$WorkDir
	DataFilelist <- params$DataFiles
	codeDir<-params$codeDir
	nNode<-as.integer(params$nNode)
	clustType<-params$clustType
	
	c1<-makeCluster(nNode,type=clustType)
	
	clusterExport(c1,"TICDenoising")#assign the current function to all nodes
	clusterApply(c1,DataFilelist,TICDenoising,WorkDir,codeDir,isBslnCrt=T)#parallel version
	
	stopCluster(c1)
	time2<-Sys.time()
	time2-time1
}

parDenoising<-function(params) 
{
	nNode<-params$nNode
	clustType<-params$clustType
	codeDir<-params$codeDir
	WorkDir<-params$WorkDir
	DataFilelist<-params$DataFiles

	cat("start denoising EIC data...\n")
	time1<-Sys.time()
	library("snow")

	c1<-makeCluster(nNode,type=clustType)

	## assign the current function to all nodes
	clusterExport(c1,"Denoising")
	## parallel version
	clusterApply(c1,DataFilelist,Denoising,WorkDir,codeDir,isBslnCrt=T,isSm=T,params)
	
	stopCluster(c1)
	time2<-Sys.time()
	cat("end denoising EIC data...\n")
	
	time2-time1
}

getPeaks<-function(vecInt,params,mz=NULL,mzVec=NULL)
{
	library("wmtsa")
	########################################
	#EIC- single mass
	#Get int vector from the mz position
	########################################
	if(!is.null(mz))
	{
		totalscan<-length(vecInt)/length(mzVec)
		mzInd<-which(mzVec==mz)
		startInd<-(mzInd-1)*totalscan+1
		endInd<-startInd+totalscan-1
		curVecInt <- vecInt[startInd:endInd]
		BHR<-as.numeric(params$BHR_EIC)#0.3##boundary/height raio
		edgeHightDiffRatio<-as.numeric(params$EHR_EIC)#0.2
		offset1<-startInd-1
		
	}else
	{
		########################################
		#TIC
		#int vector is directly from parameter
		########################################
		curVecInt<-vecInt
		totalscan<-length(curVecInt)
		BHR<-as.numeric(params$BHR_TIC)#0.4
		edgeHightDiffRatio<-as.numeric(params$EHR_TIC)#0.3
	}
	
	
	WorkDir<-params$WorkDir
	nPoints<-length(curVecInt)
	delaytime<-params$delaytime
	ScanInterval<-params$ScanInterval
	
	#################
	#apex detection
	#################
	
	#peakSpan<-min(9,ifelse(nPoints%%2==1,nPoints,nPoints-1))#peak detection
	peakSpan<-as.integer(params$Peak_span) #9
	isPeak <- peaks(x=curVecInt, span=peakSpan)
	peakInd<-which(isPeak==TRUE)#				
	
	###################
	#valley detection
	###################
	intCutoff<-0
	valleySpan<-as.integer(params$Valley_span)#5
	isMin <- peaks(-curVecInt, valleySpan)
	
	##those low-int points as valleys too
	##since sometime the consecutive low points may
	##not be detected by local minima if they have the same values
	isZeroThrsh<-curVecInt<=intCutoff
	isValley <- isMin|isZeroThrsh
	valleyInd<-which(isValley==TRUE)	
	
	##create a data frame with peak index 
 	##and intensity, ordered by peak index
	curPeakList<-data.frame(Intensity=curVecInt[peakInd],pkInd=peakInd)
	curPeakList<-curPeakList[order(peakInd),]
	
	######################
	#boudnary detection
	######################
	curPeakList$lboundInd<-0
	curPeakList$rboundInd<-0
	curPeakList$isApex<-0
	curPeakList$isShared<-0
	maxWindowlength<-as.integer(params$MaxWindow_length)#350
	maxlocalPkInd<-0#init the local peak apex position
	lbound<-0
	rbound<-0#init the rbound 
	preRbound<-0
	StN_Th <- as.integer(params$StN_Th1)
	
	for(i in 1:nrow(curPeakList))
	{
		curPeakInd<-curPeakList$pkInd[i]
		curPeakHight<-curPeakList$Intensity[i]
		
		##check if current peak apex is within the previous right boundary
		##if so, skip it since it has already been merged
		if(curPeakInd>preRbound)
		{
			########################################
			#lbound detection search lbound between 
			#previous rbound and current peak apex
			########################################
			leftValleyIndVec<-valleyInd[valleyInd<curPeakInd&valleyInd>=preRbound]
			if(length(leftValleyIndVec)==0)
			{
				##no left bound
				curPeakList[i,]$isApex=-2
			}else
			{
				##lboudnary has to satisfy a certain ratio of left boundary/height
				leftValleyIndVec1<-leftValleyIndVec[curVecInt[leftValleyIndVec]/curPeakHight<=BHR]
				lbound<-ifelse(length(leftValleyIndVec1)>0,min(leftValleyIndVec1),min(leftValleyIndVec))
				lbound<-max(lbound,0)
				lboundHeight<-curVecInt[lbound]
				curPeakList[i,]$lboundInd<-lbound
				
				##IF it's the last peak apex,
				##OR simply no value on the right side of peak
				##THEN directly set rbound as the last scan
				if(i==nrow(curPeakList)||length(valleyInd[valleyInd>curPeakInd])==0)
				{
					rbound=nPoints
					localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
					maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
					curPeakList[i,]$rboundInd<-rbound
					if(length(localPkInd)==1)
					{
						curPeakList[i,]$isApex=1
					}else
					{
						##merge the peak apexes within current boundarys.
						curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-2
						curPeakList[curPeakList$pkInd%in%localPkInd,]$lboundInd=lbound
						curPeakList[curPeakList$pkInd%in%localPkInd,]$rboundInd=rbound
						curPeakList[curPeakList$pkInd%in%localPkInd,]$isShared<-1
						##define the 1st and last peaks as -1
						curPeakList[curPeakList$pkInd%in%localPkInd,][c(1,length(localPkInd)),]$isApex=-1						
						pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
						curPeakList[pkID,]$isApex=1
					}
					preRbound<-rbound
					
				}else
				{
					########################################
					#rbound detection within the range of 
					#from current peakind to maxWindowlength
					########################################
					rightValleyIndVec<-valleyInd[valleyInd>curPeakInd&valleyInd<lbound+maxWindowlength]
					if(length(rightValleyIndVec)>0)
					{
						##filtered by peak rbound/height ratio
						rightValleyIndVec1<-rightValleyIndVec[curVecInt[rightValleyIndVec]/curPeakHight<=BHR]
						if(length(rightValleyIndVec1)>0)
						{
							##filtered by lbound-rbound/height ratio
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
								}else
								{
									##merge the peak apexes within current boundarys
									curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-2
									curPeakList[curPeakList$pkInd%in%localPkInd,]$lboundInd=lbound
									curPeakList[curPeakList$pkInd%in%localPkInd,]$rboundInd=rbound
									curPeakList[curPeakList$pkInd%in%localPkInd,]$isShared<-1
									# define the 1st and last peaks as -1
									curPeakList[curPeakList$pkInd%in%localPkInd,][c(1,length(localPkInd)),]$isApex=-1				
									pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
									curPeakList[pkID,]$isApex=1
								}
								
							}else
							{	
								###############################################
								#if edge difference is significant,try to merge
								###############################################
								rbound<-rightValleyIndVec1[1]
								
								##if lbound higher than rbound, then merge to previous peak
								if(maxlocalPkInd>0&&(lboundHeight>curVecInt[rbound]))
								{
									previousLbound<-curPeakList[curPeakList$pkInd==maxlocalPkInd,]$lboundInd	
									
									## within maximal window width
									if((rbound-previousLbound)<=maxWindowlength)
									{
										localPkInd<-peakInd[peakInd>previousLbound&peakInd<rbound]									
										curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-2
										curPeakList[curPeakList$pkInd%in%localPkInd,]$lboundInd=previousLbound
										curPeakList[curPeakList$pkInd%in%localPkInd,]$rboundInd=rbound
										curPeakList[curPeakList$pkInd%in%localPkInd,]$isShared<-1	
										# define the 1st and last peaks as -1
										curPeakList[curPeakList$pkInd%in%localPkInd,][c(1,length(localPkInd)),]$isApex=-1								
										maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
										pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
										curPeakList[pkID,]$isApex=1
									}else
									{
										##window exceed max size, then simply 
										##keep the closed boundary as the final cut
										localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
										maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
										curPeakList[i,]$rboundInd<-rbound
										if(length(localPkInd)==1)
										{
											curPeakList[i,]$isApex=1
										}else
										{
											##merge the peak apexes within current boundarys
											curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-2
											curPeakList[curPeakList$pkInd%in%localPkInd,]$lboundInd=lbound
											curPeakList[curPeakList$pkInd%in%localPkInd,]$rboundInd=rbound
											curPeakList[curPeakList$pkInd%in%localPkInd,]$isShared<-1											
											# define the 1st and last peaks as -1
											curPeakList[curPeakList$pkInd%in%localPkInd,][c(1,length(localPkInd)),]$isApex=-1											
											pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
											curPeakList[pkID,]$isApex=1
										}
									}
									preRbound<-rbound
									
								}else
								{
									##if lbound lower than rbound, then merge to next peak
									if((rbound-lbound)<=maxWindowlength)
									{
										curPeakList[i,]$isApex<--1
										
									}else
									{
										##if window exceeds max size, then force 
										##the current rbound as the final boundary
										localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
										maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
										curPeakList[i,]$rboundInd<-rbound
										if(length(localPkInd)==1)
										{											
											curPeakList[i,]$isApex=1
										}else
										{
											##merge the peak apexes within current boundarys
											curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-2
											curPeakList[curPeakList$pkInd%in%localPkInd,]$lboundInd=lbound
											curPeakList[curPeakList$pkInd%in%localPkInd,]$rboundInd=rbound
											curPeakList[curPeakList$pkInd%in%localPkInd,]$isShared<-1
											# define the 1st and last peaks as -1
											curPeakList[curPeakList$pkInd%in%localPkInd,][c(1,length(localPkInd)),]$isApex=-1											
											pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
											curPeakList[pkID,]$isApex=1
										}								
										preRbound<-rbound
									}
								}
							}
						}else
						{
							#######################################
							#if no rbound pass the pkheight ratio 
							#filter, then check the closed rbound
							#######################################
							rbound<-rightValleyIndVec[1]
							
							##if edge difference is not significant
							if(abs(curVecInt[rbound]-lboundHeight)/curPeakHight<=edgeHightDiffRatio)
							{
								curPeakList[i,]$rboundInd<-rbound
								preRbound<-rbound							
								localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
								maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
								if(length(localPkInd)==1)
								{									
									curPeakList[i,]$isApex=1
								}else
								{
									##merge the peak apexes within current boundarys
									curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-2
									curPeakList[curPeakList$pkInd%in%localPkInd,]$lboundInd=lbound
									curPeakList[curPeakList$pkInd%in%localPkInd,]$rboundInd=rbound
									curPeakList[curPeakList$pkInd%in%localPkInd,]$isShared<-1								
									# define the 1st and last peaks as -1
									curPeakList[curPeakList$pkInd%in%localPkInd,][c(1,length(localPkInd)),]$isApex=-1									
									pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
									curPeakList[pkID,]$isApex=1
								}
							}else
							{
								##if edge diff is significant,then try to merge
								##if lbound higher than rbound, merge to previous peak
								if(maxlocalPkInd>0&&(lboundHeight>curVecInt[rbound]))
								{
									previousLbound<-curPeakList[curPeakList$pkInd==maxlocalPkInd,]$lboundInd		
									if((rbound-previousLbound)<=maxWindowlength)
									{
										localPkInd<-peakInd[peakInd>previousLbound&peakInd<rbound]										
										curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-2
										curPeakList[curPeakList$pkInd%in%localPkInd,]$lboundInd=previousLbound
										curPeakList[curPeakList$pkInd%in%localPkInd,]$rboundInd=rbound
										curPeakList[curPeakList$pkInd%in%localPkInd,]$isShared<-1										
										# define the 1st and last peaks as -1
										curPeakList[curPeakList$pkInd%in%localPkInd,][c(1,length(localPkInd)),]$isApex=-1										
										maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
										pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
										curPeakList[pkID,]$isApex=1
									}else
									{
										##if window exceeds max size, simply keep 
										##the closed boundary as the final cut
										localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
										maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
										curPeakList[i,]$rboundInd<-rbound
										if(length(localPkInd)==1)
										{										
											curPeakList[i,]$isApex=1
										}else
										{
											##merge the peak apexes within current boundarys
											curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-2
											curPeakList[curPeakList$pkInd%in%localPkInd,]$lboundInd=lbound
											curPeakList[curPeakList$pkInd%in%localPkInd,]$rboundInd=rbound
											curPeakList[curPeakList$pkInd%in%localPkInd,]$isShared<-1
											# define the 1st and last peaks as -1
											curPeakList[curPeakList$pkInd%in%localPkInd,][c(1,length(localPkInd)),]$isApex=-1											
											pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
											curPeakList[pkID,]$isApex=1
										}
									}
									preRbound<-rbound
								}else
								{
									##if lbound lower than rbound, then merge to next peak
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
											curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-2
											curPeakList[curPeakList$pkInd%in%localPkInd,]$lboundInd=lbound
											curPeakList[curPeakList$pkInd%in%localPkInd,]$rboundInd=rbound
											curPeakList[curPeakList$pkInd%in%localPkInd,]$isShared<-1
											# define the 1st and last peaks as -1
											curPeakList[curPeakList$pkInd%in%localPkInd,][c(1,length(localPkInd)),]$isApex=-1
											
											pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
											curPeakList[pkID,]$isApex=1											
										}								
										preRbound<-rbound
									}
								}
							}
						}
						
					}else
					{
						##when the nearest rbound exceed the max 
						#window length, just pick the nearest one
						rbound<-min(valleyInd[valleyInd>curPeakInd])
						curPeakList[i,]$rboundInd<-rbound
						
						localPkInd<-peakInd[peakInd>lbound&peakInd<rbound]
						maxlocalPkInd<-localPkInd[which.max(curVecInt[localPkInd])]
						if(length(localPkInd)==1)
						{
							
							curPeakList[i,]$isApex=1
						}else
						{
							curPeakList[curPeakList$pkInd%in%localPkInd,]$isApex=-2
							curPeakList[curPeakList$pkInd%in%localPkInd,]$lboundInd=lbound
							curPeakList[curPeakList$pkInd%in%localPkInd,]$rboundInd=rbound
							curPeakList[curPeakList$pkInd%in%localPkInd,]$isShared<-1							
							# define the 1st and last peaks as -1
							curPeakList[curPeakList$pkInd%in%localPkInd,][c(1,length(localPkInd)),]$isApex=-1						
							pkID<-rownames(curPeakList[curPeakList$pkInd==maxlocalPkInd,])
							curPeakList[pkID,]$isApex=1
						}										
						preRbound<-rbound					
					}
				}		
			}
		}
	} # for loop
	
	#cat(paste("finish peak picking\n",sep=""))
	
	#########################
	# Outputing final results
	#########################
	
	if(!is.null(mz))
	{
		##--------------------------------------------------------------------------------
		## Small feature based calculation
		
		#########################
		#EIC peak picking output
		#########################
		
#		curPeakList$curMz=mz
#		
#		##############################################
#		#save offset of the current EIC block position 
#		#in the long intensity vector (vecInt) for 
#		#each EIC peak,so that peak profile could  
#		#be easily indexed from vecInt in decon
#		################################################
#		
#		curPeakList$offset<-offset1
#		curPeakList<-subset(curPeakList,rboundInd!=0)
#		curPeakList$StN<-0
#		
#		## identify unique peak feature by searching unique
#		## left boundary
#		unique.length <- length(unique(curPeakList$lboundInd))
#		index <- NULL
#		for (i in curPeakList$lboundInd) {
#			index <- c(index,which(unique(curPeakList$lboundInd)==i))
#		}
#		## index - unique peak feature ID
#		curPeakList$index <- index
#		
#		## to calcuate StN
#		for(i in c(1:unique.length)) {
#			cur.rows <- which(curPeakList$index==i)
#			curPk<-curPeakList[cur.rows[1],]
#			startInd<-curPk$lboundInd
#			endInd<-curPk$rboundInd
#			curProfile<-vecInt[(startInd+offset1):(endInd+offset1)]
#			## to identify the apex position (idex) of a peak feature
#			apex.vector <- which(c(startInd:endInd)%in%curPeakList[cur.rows,]$pkInd==TRUE)
#			apex.length <- length(apex.vector)
#			
#			#if (length(curProfile[curProfile!=0])>=5&length(curProfile)>=6)
#		if (length(curProfile[curProfile>=20])!=0&(length(curProfile[curProfile!=0])>=5&length(curProfile)>=15))
#			{
#				cwt <- wavCWT_Modify(curProfile)
#				tree <- wavCWTTree_Modify(cwt,type = "maxima")
#				
#				if (!is.null(tree)) {
#					WT_result <- wavCWTPeaks_Modify(tree)
#					if (length(WT_result$x)!=0) {
#						## extract identified peak apex pos and corresponding Signal to Noise
#						Time <- WT_result$x
#						sigToNoise <- as.numeric(attr(WT_result,"snr"))
#						Time.length <- length(Time)
#						
#						## 2 condition
#						## 1st: apex candidates are more than identified peaks from CentWave
#						if (apex.length <= Time.length) {
#							j =1
#							while (j <= apex.length) {
#								apex.temp <- apex.vector[j]
#								if (min(abs(Time-apex.temp)) <= 15) {
#									pos <- which.min(abs(Time-apex.temp))
#									curPeakList[cur.rows[j],]$StN <- sigToNoise[pos]
#									Time <- Time[-pos]
#									sigToNoise <- sigToNoise[-pos]
#								}
#								j = j+1				
#							}
#						}else {
#							## 2nd: identified peaks from CentWave are more than apex candidates 
#							z =1
#							while (z <= Time.length) {
#								time.temp <- Time[z]
#								if (min(abs(apex.vector-time.temp)) <= 15) {
#									pos <- which.min(abs(apex.vector-time.temp))
#									curPeakList[cur.rows[pos],]$StN <- sigToNoise[z]
#									apex.vector <- apex.vector[-pos]
#									cur.rows <- cur.rows[-pos]
#								}
#								z = z+1				
#							}
#						}	
#					}
#				}	
#			}	
#		}#for loop to cal StN
		

		##--------------------------------------------------------------------------------
		## whole feature based calculation

		curPeakList$curMz=mz
		curPeakList$offset<-offset1
		curPeakList<-subset(curPeakList,rboundInd!=0)
		curPeakList$StN<-0
		
#		startInd<-min(curPeakList$lboundInd)
#		endInd<-max(curPeakList$rboundInd)
#		curProfile<-vecInt[(startInd+offset1):(endInd+offset1)]
		cwt <- wavCWT_Modify(curVecInt)
		tree <- wavCWTTree_Modify(cwt,type = "maxima") 
		WT_result <- wavCWTPeaks_Modify(tree)
		RTtime <- WT_result$x
		sigToNoise <- as.numeric(attr(WT_result,"snr"))
		
		apex.length <- nrow(curPeakList)
		Time.length <- length(RTtime)
		apex.vector <- curPeakList$pkInd
		
		if (apex.length <= Time.length) {
			j =1
			while (j <= apex.length) {
				apex.temp <- apex.vector[j]
				if (min(abs(RTtime-apex.temp)) <= 20) {
					pos <- which.min(abs(RTtime-apex.temp))
					curPeakList[j,]$StN <- sigToNoise[pos]
					RTtime <- RTtime[-pos]
					sigToNoise <- sigToNoise[-pos]
				}
				j = j+1				
			}
		}else {
			z =1
			while (z <= Time.length) {
				time.temp <- RTtime[z]
				if (min(abs(apex.vector-time.temp)) <= 20) {
					pos <- which.min(abs(apex.vector-time.temp))
					curPeakList[which(curPeakList$pkInd == apex.vector[pos]),]$StN <- sigToNoise[z]
					apex.vector <- apex.vector[-pos]
				}
				z = z+1				
			}
		}#else

#		# *Modification 3*
#		# I notice there're very big peak assigned StN=0, which made wrong identification
#		# new algorithm: peak window based, locate each merged peak window, compare the wavelet
#		# result with ADAP, and assign StN to each identified peak
#		
#		unique_window_id <- unique(curPeakList$lboundInd)
#		updatePeakList <- NULL
#		
#		for (rowindex in 1:length(unique_window_id)) {		
#			cur_peak <- subset(curPeakList,lboundInd==unique_window_id[rowindex])
#			wave_index <- which(RTtime>=cur_peak$lboundInd[1]&RTtime<cur_peak$rboundInd[1])
#			
#			if (length(wave_index)==0) {
#				updatePeakList <- rbind(updatePeakList,cur_peak)
#				next
#			}
#			
#			else {
#				RTtime_wav <- RTtime[wave_index]
#				StN_wav <- sigToNoise[wave_index]
#				apex.length <- nrow(cur_peak)
#				StN_wav.length <- length(RTtime_wav)
#				apex.vector <- cur_peak$pkInd
#				
#				if (apex.length==1&StN_wav.length==1) {
#					cur_peak$StN <- StN_wav
#				} else if (apex.length <= StN_wav.length) {
#					j =1
#					while (j <= apex.length) {
#						apex.temp <- apex.vector[j]
#						pos <- which.min(abs(RTtime_wav-apex.temp))						
#						if (min(abs(RTtime_wav[pos]-apex.vector)) < min(abs(RTtime_wav-apex.temp))) {
#							j = j+1
#							next
#						}else {
#							cur_peak[j,]$StN <- StN_wav[pos]
#							RTtime_wav <- RTtime_wav[-pos]
#							StN_wav <- StN_wav[-pos]	
#							j = j+1	
#						}
#						
#					}
#				} else if (apex.length > StN_wav.length) {
#					z =1
#					while (z <= StN_wav.length) {
#						time.temp <- RTtime_wav[z]					
#						pos <- which.min(abs(apex.vector-time.temp))
#						if (min(abs(RTtime_wav-apex.vector[pos])) < min(abs(apex.vector-time.temp))) {
#							z = z+1	
#							next
#						}					
#						else {
#							cur_peak[which(cur_peak$pkInd == apex.vector[pos]),]$StN <- StN_wav[z]
#							apex.vector <- apex.vector[-pos]
#							z = z+1
#						}									
#					}
#				}#else	
#			}
#			nonzero_stn_peak <- which(cur_peak$StN!=0)
#			zero_stn_peak <- which(cur_peak$StN==0)
#			if (length(nonzero_stn_peak) >0&length(nonzero_stn_peak)>0) {
#				for (this.pk.id in zero_stn_peak) {
#					close_dis <- min(abs(cur_peak[-this.pk.id,]$pkInd)-cur_peak[this.pk.id,]$pkInd)
#					if (close_dis<=10){
#						cur_peak[this.pk.id,]$StN <- cur_peak[which.min(abs(cur_peak[-this.pk.id,]$pkInd)-cur_peak[this.pk.id,]$pkInd),]$StN
#					}
#				}
#			}
#			updatePeakList <- rbind(updatePeakList,cur_peak)
#		}
		
		# As we couldn't promise all StN assignment is right, so just remove those peaks
		# with zero temporirally

		##--------------------------------------------------------------------------------
		## calculate the gss and sharpness for each EIC peak profile
		#curPeakList<-subset(updatePeakList,isApex!=-2&StN>=10,select=c("curMz","pkInd","Intensity","lboundInd","rboundInd","isShared","isApex","offset","StN"))
		curPeakList<-subset(curPeakList,isApex!=-2&StN>=StN_Th,select=c("curMz","pkInd","Intensity","lboundInd","rboundInd","isShared","isApex","offset","StN"))
		if (nrow(curPeakList) >0) { # to check if it's empty when one EIC is pretty noisy, modified May 14, 2012
			curPeakList$gss<-90
			curPeakList$gssPos<-0
			curPeakList$shrp<-0
			
			for(i in 1:nrow(curPeakList))
			{
				curPk<-curPeakList[i,]
				startInd<-curPk$lboundInd
				endInd<-curPk$rboundInd
				curProfile<-vecInt[(startInd+offset1):(endInd+offset1)]			
				##to identify the apex position (index)
				apex <- which(c(startInd:endInd)==curPk$pkInd)
				##gaussian fitting of two mirror image profiles
				gss<-GaussianDisimilarity1(curProfile,apex)
				gssPos<-which.min(gss)
				curPeakList[i,]$gssPos<-gssPos
				curPeakList[i,]$gss<-gss[gssPos]
				##get sharpness value of selected image profile
				shrp<-Sharpness1(yy=curProfile,apex)	
				curPeakList[i,]$shrp<-shrp[gssPos]			
			}
			# remove those mass with no guassian shape
			curPeakList<-subset(curPeakList,gss!=90) 
			curPeakList
		}
	}else
	{
		#########################
		#TIC peak picking output
		#########################
		subset(curPeakList,select=c("pkInd","lboundInd","rboundInd","isShared","isApex"))
	}
	
}

parTICPeakpicking<-function(params,denoised=T)
{
	time1<-Sys.time()
	
	DataFilelist <- params$DataFiles
	nNode<-as.integer(params$nNode)
	clustType<-params$clustType
	
	library("snow")
	cl<-makeCluster(nNode,type=clustType)

	clusterExport(cl,"peaks")#assign the current profiles to all nodes
	clusterExport(cl,"getPeaks")#assign the current profiles to all nodes
	clusterExport(cl,"parseFileName")#assign the current profiles to all nodes

	##calculate the distance between each pair of profiles
	clusterApply(cl,DataFilelist,TICPeakpicking,params,denoised)#parallel version
	stopCluster(cl)
	time2<-Sys.time()
	time2-time1
}

TICPeakpicking<-function(inFilePath,params,denoised=T)
{
	
	library("ncdf")
	WorkDir <- params$WorkDir
		
	###############
	#read TIC data
	################
	fileName<-parseFileName(inFilePath)
	denoisedTICFile<-paste(WorkDir,"/output/TIC/denoised_",fileName,"_TIC.cdf",sep="")
	
	if(file.exists(inFilePath)==FALSE)next
	if(denoised)
	{
		cat(fileName,"reading denoised TIC data...\n")
		ncid <- open.ncdf(denoisedTICFile)
	}else
	{
		cat(fileName,"reading raw TIC data...\n")
		ncid <- open.ncdf(inFilePath)
	}
	
	vecInt <- get.var.ncdf(ncid, varid="total_intensity")
	close.ncdf(ncid)
	remove(ncid)	
	
	###############
	#peak picking
	###############
	cat(fileName,"TIC peak picking...")
	PeakList<-getPeaks(vecInt,params,mz=NULL,mzVec=NULL)
	## save peak picking results to cdf
	write.csv(PeakList,file=paste(WorkDir,"/output/peakpicking/",fileName,"_TIC_PeakList.csv",sep=""),row.names = F)
}


getPeaksGroup<-function(mzGroup,vecInt,mzVec,params,fileName)
{
	nMz<-length(mzGroup)
	groupResult<-vector("list",nMz)
	sink(paste(params$WorkDir,fileName,".txt",sep=""),append=T)	
	
	for(i in 1:nMz)
	{
		mz<-mzGroup[i]	
		print(mz)
		if (mz==0) next
		groupResult[[i]]<-getPeaks(vecInt,params,mz,mzVec)
	}
	sink()
	groupResult
}

parEICpeakpicking<-function(params)
{
	time1<-Sys.time()
	###############
	#get parameters
	###############
	#library(tcltk)
	library("snow")
	library("ncdf")
	library("wmtsa")
	
	WorkDir <- params$WorkDir
	DataFilelist <- params$DataFiles
	nNode<-as.integer(params$nNode)
	clustType<-params$clustType
	
	cl<-makeCluster(nNode,type=clustType)
	
	###############
	#peak picking
	###############
		
	for(fileindex in 1:length(DataFilelist))
	{
		
		###############
		#read EIC data
		################
		inFilePath <- DataFilelist[[fileindex]]
		fileName<-parseFileName(inFilePath)
		
		cat(fileName,"reading EIC data...\n")
		#EICFile<-paste(WorkDir,"/output/EIC/",fileName,"EIC.cdf",sep="")
	
		#check if it exist denoised EIC file
		rawEICFile <-paste(WorkDir,"output/EIC/",fileName,"EIC.cdf",sep="")
		denoiseEICFile<-paste(WorkDir,"output/EIC/denoised_",fileName,"EIC.cdf",sep="")
		EICFile<-ifelse(file.exists(denoiseEICFile),denoiseEICFile,rawEICFile)

		if(file.exists(inFilePath)==FALSE)next
		ncid <- open.ncdf(EICFile)
		vecInt <- get.var.ncdf(ncid, varid="intVec")
		mzVec<-get.var.ncdf(ncid, varid="mzVec")
	
		close.ncdf(ncid)
		remove(ncid)	
		
		###############
		#peak picking
		################
		cat(fileName,"peak picking...")
		#######peak picking
		mzGroups<-clusterSplit(cl,mzVec)
		
		clusterExport(cl,"peaks")#assign the current profiles to all nodes
		clusterExport(cl,"getPeaks")#assign the current profiles to all nodes
		clusterExport(cl,"GaussianDisimilarity")#assign the current profiles to all nodes
		clusterExport(cl,"GaussianDisimilarity1")#assign the current profiles to all nodes
		clusterExport(cl,"Sharpness")#assign the current profiles to all nodes
		clusterExport(cl,"Sharpness1")#assign the current profiles to all nodes
		clusterExport(cl,"dotProduct")#assign the current profiles to all nodes
		clusterExport(cl,"wavCWT_Modify")#assign the current profiles to all nodes
		clusterExport(cl,"wavCWTTree_Modify")#assign the current profiles to all nodes
		clusterExport(cl,"wavCWTPeaks_Modify")#assign the current profiles to all nodes
		time1<-Sys.time()
		PeakList<-clusterApply(cl,mzGroups,getPeaksGroup,vecInt,mzVec,params,fileName)#parallel version
		time2<-Sys.time()
		time2-time1
		PeakList<-unlist(PeakList,recursive=F)
		PeakList<-do.call(rbind,PeakList)

		write.csv(PeakList,file=paste(WorkDir,"/output/peakpicking/",fileName,"_EIC_PeakList.csv",sep=""),row.names=F)
		
	}
	stopCluster(cl)
	
	time2<-Sys.time()
	time2-time1
}

deconvolution<-function(params,isDistParallel=T,clustingType="h",isDataParalell)
{
	###############
	#get parameters
	################
	library(playwith)
	library("snow")
	DataFilelist<-params$DataFiles
	
	# Decide the optimal number of nodes to be used
	# yan modified on Mar 13
	params$nNode <- min(as.integer(params$nNode),length(DataFilelist))
	
	#######
	#decon
	#######
	cl<-makeCluster(params$nNode,type=params$clustType)
	
	if(isDataParalell==T)
	{
		time1<-Sys.time()
		clusterExport(cl,"decomposition")#assign the current profiles to all nodes
		clusterApply(cl,DataFilelist,decomposition,params,cl,isDistParallel=F,clustingType)#parallel version
		time2<-Sys.time()
		time2-time1
	}else
	{
		time1<-Sys.time()
#		fileindex<-4
#		isDistParallel=T
		clustingType="h"
		for(fileindex in 1:length(DataFilelist))
		{
			inFilePath <- DataFilelist[[fileindex]]
			
			decomposition(inFilePath,params,cl,isDistParallel=T,clustingType)
		}
		
		time2<-Sys.time()
		time2-time1
		
	}	
	
	stopCluster(cl)
	
}

decomposition_old<-function(inFilePath,params,cl,isDistParallel,clustingType)
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
	lib<-readMSP2Spec(filename="/home/ptpham/StdsLibrary/KQC.txt",withRT=T)
	
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

decomposition<-function(inFilePath,params,cl,isDistParallel,clustingType)
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
	libFile <- params$libFile_decon
	StN_Th <- as.integer(params$StN_Th2)
	Shrp_Th <- as.numeric(params$Sharpness_Th)
	Gss_Th <- as.numeric(params$Gss_Th)	
	TotalScore_factor <- as.numeric(params$TotalScore_factor)
	ClusterDist_Th <- as.numeric(params$ClusteringDist_Th)
	SpecSimilarity_Th <- as.numeric(params$IntroSpecSimilarity_Th)
	libMatch_Th <- as.numeric(params$libMatching_Th1)
	NonUmass1 <- as.numeric(unlist(strsplit(params$Non_Umass_above_equal_8_min,",")))
	NonUmass2 <- as.numeric(unlist(strsplit(params$Non_Umass_below_8_min,",")))
			
	##load library
	#lib<-readMSP2Spec(filename="../StdsLibrary/QcCompound.txt",withRT=T)
	lib<-readMSP2Spec(filename=libFile,withRT=T)
	
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
	
	#To check those potential compounds only
#	RTVector<-read.csv(paste(WorkDir,"output/RTlist.csv",sep=""))$ET
#	TICRTvec<-(TICpeaklist$pkInd-1)*ScanInterval+delaytime
#	winIDs<-NULL		
#	for(i in RTVector)
#	{
#		winIDs<-c(winIDs,which(abs(TICRTvec-i)<=0.2))
#	}
#	winIDs<-sort(unique(winIDs))
	
	#############################
	#read EIC peak picking result
	#############################
	cat(fileName,"reading EIC peak data...\n")
	EICPeakFile<-paste(WorkDir,"/output/peakpicking/",fileName,"_EIC_PeakList.csv",sep="")
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
	outputPath<-paste(paste(params$WorkDir,"output/decon/",sep=""),paste(fileName,"_ModelPkProfiles.txt",sep=""),sep="")
	strOutput <- NULL
	compNum <- 1
#	modelPkList3s <-NULL
#	modelPkList2s <-NULL
#	modelPkList1s <-NULL
	
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
			CurNonUmassVec<-NonUmass1
		}else
		{
			CurNonUmassVec<-NonUmass2
		}
		
		##################################################
		#get all the EIC peaks within current decon window
		##################################################
		localTICApexList<-subset(TICApexList,pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		localEICPeakList<-subset(EICpeaklist,flag==0&pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		
		if (nrow(localEICPeakList)!= 0) {
			
			if (max(localEICPeakList$Intensity)>=100) {
				##keep the unique mass profile of local EIC list
				##original local EIC list consists of merged peaks
				localEICPeakList <- subset(localEICPeakList,isApex==1&Intensity>=100)
				
				if(nrow(localEICPeakList)>1)
				{
					
					#localEICPeakList$Sig = 0
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
#						if (mean(curProfile$int)!=0) {
#							Sig <- var(curProfile$int)/mean(curProfile$int)
#							localEICPeakList[i,]$Sig <- Sig
#						}
						allprofiles[[i]]<-curProfile
						names(allprofiles)[i]<-rownames(curEICpk)	
						
					}		
					
					#CandidatePeaklist<-subset(localEICPeakList,StN>=300&shrp>=5&Sig>=300&gss<=3)
					CandidatePeaklist<-subset(localEICPeakList,StN>=StN_Th&shrp>=Shrp_Th&gss<=Gss_Th)
					############################
					##scoring system
					############################	
					
					if (nrow(CandidatePeaklist)>0) {
						#mass score
						CandidatePeaklist$f1<-scale(CandidatePeaklist$mz,min(CandidatePeaklist$mz),diff(range(CandidatePeaklist$mz)))	
						##gassian similarity score
						CandidatePeaklist$f2<-scale(cos(CandidatePeaklist$gss*pi/180),min(cos(CandidatePeaklist$gss*pi/180)),diff(range(cos(CandidatePeaklist$gss*pi/180))))
						## Signal to Noise
						CandidatePeaklist$f3<-scale(log(CandidatePeaklist$StN),min(log(CandidatePeaklist$StN)),diff(range(log(CandidatePeaklist$StN))))
						CandidatePeaklist$f4<-scale(log(CandidatePeaklist$Intensity),min(log(CandidatePeaklist$Intensity)),diff(range(log(CandidatePeaklist$Intensity))))
						## Singnificance level
						#CandidatePeaklist$f4<-scale(log(CandidatePeaklist$Sig),min(log(CandidatePeaklist$Sig)),diff(range(log(CandidatePeaklist$Sig))))
						CandidatePeaklist$f1[is.na(CandidatePeaklist$f1)]<-0
						CandidatePeaklist$f2[is.na(CandidatePeaklist$f2)]<-0
						CandidatePeaklist$f3[is.na(CandidatePeaklist$f3)]<-0
						CandidatePeaklist$f4[is.na(CandidatePeaklist$f4)]<-0
						##linear combination of 3 scores
						#scores<-(10/7)*(0.1*CandidatePeaklist$f1+0.3*CandidatePeaklist$f2+0.2*CandidatePeaklist$f3+0.1*CandidatePeaklist$f4)
						scores<-(10/8)*(0.1*CandidatePeaklist$f1+0.3*CandidatePeaklist$f2+0.2*CandidatePeaklist$f3+0.2*CandidatePeaklist$f4)
						#scores<-(10/7)*(0.1*CandidatePeaklist$f1+0.3*CandidatePeaklist$f2+0.3*CandidatePeaklist$f3)
						CandidatePeaklist$score<-scale(scores,F,0.001)
						
#						a <- hist(CandidatePeaklist$score)
#						plot(a$breaks[-1],a$counts/sum(a$counts),type="h")
#						
						
						## Create good peaks for clustering
						goodShapePeaklist<-subset(CandidatePeaklist,score>=min(CandidatePeaklist$score)+diff(range(CandidatePeaklist$score))*TotalScore_factor)
						goodShapePeaklist<-subset(goodShapePeaklist,!mz%in%CurNonUmassVec)
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
								maxIntraDist<- ClusterDist_Th #cutoff to decide if split
								
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
									modelPkList<- goodShapePeaklist[which.max(goodShapePeaklist$score),]
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
#						modelPkList3 <- NULL
#						modelPkList2 <- NULL
#						modelPkList1 <- NULL
							##No co-eluting compounds (only one or no model peak detected)
							##Then extract spectrum directly
							if(is.null(modelPkList))
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
											area=as.integer(),windowID=as.integer(),compoundID=as.integer(),isModel=as.integer(),isUnique=as.integer(),factor=as.numeric(),compID=as.integer())
									
									##name spec list with time and mz
#									clustID<-format((maxInd-1)*ScanInterval+delaytime,digit=4)
#									names(specList)[1]<-paste(clustID,curMz)
									
									names(specList)[1]<-paste(maxInd,windowID,1,curMz,compNum)
									
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
															pkInd=maxInd,Intensity=curInt,area=curArea,windowID=windowID,compoundID=1,isModel=0,isUnique=0,factor=0,compID=compNum))
										}
									}
									
									##Mark the ion with maximal intensity as model
									componentList[[1]][which.max(componentList[[1]]$Intensity),]$isModel=1
#								modelPkList1 <- curTICPk$pkInd
									
								}					
								componentList[[1]][which.max(componentList[[1]]$Intensity),]$isUnique=1
								##replace decom as 0, stop while loop
								decom=0
								compNum <- compNum+1
								
							}else
							{
#								## >= one model peaks
#								if (nrow(modelPkList)>1) {
#									rmmodelPks<-data.frame(row=as.integer(),col=as.integer())						
#									modelPkDistance <- as.matrix(modelPkDist(modelPkList,profiles))
#									
#									for (i in c(1:nrow(modelPkDistance))) {
#										subDistance <- modelPkDistance[i,]
#										par_name <- rownames(modelPkDistance)[i]
#										par_index <- modelPkList[par_name,]$pkInd
#										
#										for (j in c(1:length(subDistance))) {
#											dau_name <- names(subDistance)[j]
#											dau_index <- modelPkList[dau_name,]$pkInd
#											
#											comb_factor <- abs(par_index - dau_index)*modelPkDistance[i,j]
#											
#											if (comb_factor<= 12||(abs(par_index - dau_index)<=6&modelPkDistance[i,j]<=10)) {
#												rmmodelPks <- rbind(rmmodelPks,data.frame(row=i,col=j))
#											}	
#										}
#									}
#									
#									rmmodelPks <- subset(rmmodelPks,row!=col)
#									#######################################
#									#Step 2: grouping pairs
#									#######################################
#									if (nrow(rmmodelPks)>=1) {
#										m=c(1:nrow(rmmodelPks))
#										l = 1
#										modelGroup <- list()
#										
#										while (length(m)>0) {
#											temp <- NULL
#											for (i in m[2:length(m)]) {
#												if (length(rmmodelPks[i,][(rmmodelPks[i,]%in%rmmodelPks[m[1],])])>=1) {
#													temp <- c(temp, i) 
#												}
#											}
#											temp <- c(temp, m[1]) 
#											m <- m[!m%in%temp]
#											modelGroup[[l]] <- temp
#											l <- l+1
#										}
#										
#										###########################################
#										#Step 3: get representative for each group
#										###########################################
#										
#										modelPkList_temp <- NULL
#										rmmodelPkList_temp <- NULL
#										for (i in c(1:length(modelGroup))) {
#											rmmodelPkList <- modelPkList[unique(unlist(rmmodelPks[modelGroup[[i]],])),]
#											modelPkList_temp <- rbind(modelPkList_temp,rmmodelPkList[which.max(rmmodelPkList$score),])
#											rmmodelPkList_temp <- rbind(rmmodelPkList_temp,modelPkList[rownames(modelPkList)%in%rownames(rmmodelPkList),])
#											##Only keep the one with highest score
#											#modelPkList <- rbind(modelPkList,modelPkList_temp)	
#										}
#										
#										## get removed peaks and update model Peak list
#										rmmodelPkList_temp <- rmmodelPkList_temp[!rownames(rmmodelPkList_temp)%in%rownames(modelPkList_temp),]
#										modelPkList <- modelPkList[!rownames(modelPkList)%in%rownames(rmmodelPkList_temp),]
#										
#									}
#								}
								
								mdlPkIDs<-rownames(modelPkList)
								
								lbound<-min(localEICPeakList$lboundInd)
								rbound<-max(localEICPeakList$rboundInd)
								
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
											area=as.integer(),windowID=as.integer(),compoundID=as.integer(),isModel=as.integer(),isUnique=as.integer(),factor=as.numeric(),compID=as.integer())
									##add model peak to the spectrum
									specList[[i]]<-data.frame(mz=as.integer(),int=as.integer())
									#curRT<-format((modelPkList$pkInd[i]-1)*ScanInterval+delaytime,digit=4)
									#names(specList)[i]<-paste(curRT,modelPkList$mz[i])
									names(specList)[i]<-paste(modelPkList$pkInd[i],windowID,i,modelPkList$mz[i],compNum+i-1)
									
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
										
										if (A[m]==0&modelPkList[curmdlPkID,]$mz==curMz) {
											componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,
															rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=10,
															area=10,windowID=windowID,compoundID=j,isModel=1,isUnique=0,factor=A[m],compID=compNum+m-1))
										}
										
										if(A[m]>0)
										{
											specList[[j]]<-rbind(specList[[j]],data.frame(mz=as.integer(curMz),int=A[m,]))###add the curent restored mixing factor to spec
											if (modelPkList[curmdlPkID,]$mz==curMz) {
												componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,
																rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=as.integer(modelPkList[curmdlPkID,]$Intensity*A[m,]),
																area=as.integer(modelPkList[curmdlPkID,]$area*A[m,]),windowID=windowID,compoundID=j,isModel=1,isUnique=0,factor=A[m],compID=compNum+m-1))
											}
											
											else {
												componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,
																rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=as.integer(modelPkList[curmdlPkID,]$Intensity*A[m,]),
																area=as.integer(modelPkList[curmdlPkID,]$area*A[m,]),windowID=windowID,compoundID=j,isModel=0,isUnique=0,factor=A[m],compID=compNum+m-1))
											}
										}							
									}
								}
								
								###############################################
								#Makeup 2 for splitting issue
								#Checking pairwise resovlved spectra similarity
								###############################################
								Num <- nrow(modelPkList)									
								if (Num>1) {
									RepeatModelPk <- NULL									
									indexPairs<-combinations(Num,2)
									for(i in 1:nrow(indexPairs))
									{
										index <- indexPairs[i,]
										index1 <- index[1]
										index2 <- index[2]
										
										##calculate time difference between two model peaks
										modelPKET <- modelPkList$pkInd[c(index1,index2)]	
										
										#fragments with the maximal time difference 10 scans are considered a compound
										if (abs(modelPKET[1]-modelPKET[2]) <=10&nrow(specList[[index1]])!=0&nrow(specList[[index2]])!=0)
										{
											score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
											if (0.9*score+(1-(abs(modelPKET[1]-modelPKET[2])/10))*100 >=SpecSimilarity_Th) {
												selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index2],mdlPkIDs[index1])
												RepeatModelPk <- c(RepeatModelPk,selectedModelPkID)
											}
										}
									}						
									##select unique filtered model peaks 
									RepeatModelPk <- unique(RepeatModelPk)
									modelPkList <- modelPkList[!rownames(modelPkList)%in%RepeatModelPk,]
								}
								
								if (nrow(modelPkList) == Num) {
									decom = 0
									
									
									if (nrow(modelPkList) == 1) {
										componentList[[1]][which.max(componentList[[1]]$Intensity),]$isUnique=1
									}else if (nrow(modelPkList) == 2) {
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
									
									## output model pk profiles
									for (model in c(1:nrow(modelPkList))) {
										curModel <- modelPkList[model,]
										curModelInt <- S[,model]
										strOutput <- paste(strOutput,getModelPkProfile(curModel,curModelInt,windowID,lbound,rbound,model,compNum),sep="\n\n")
									}
									
									compNum <- compNum+Num
									
								} # Model Pk keep the same number
							} # model PK num is >2 at first
						} # while statement
						
						##flag analyzed EIC peaks 
						EICpeaklist[rownames(localEICPeakList),]$flag<-1
						componentResults<-c(componentResults,componentList)
						specResults<-c(specResults,specList)
#					modelPkList3s <- c(modelPkList3s,modelPkList3)
#					modelPkList2s <- c(modelPkList2s,modelPkList2)
#					modelPkList1s <- c(modelPkList1s,modelPkList1)
						
						##model peak id
						mdlPkIDVec<-c(mdlPkIDVec,mdlPkIDs)	
					} #candidate checking 				
				} # local Peaks checking3
			} # local peaks checking2	
		} # local Peaks checking1		
	}
	
#	tempCol <- NULL
#	for (i in c(1:length(componentResults))) {
#		tempCol <- c(tempCol,ncol(componentResults[[i]]))
#	}
#	#componentResults[[which(tempCol!=12)]] <- NULL
#	
#	if (length(which(tempCol!=12))>0) {
#		componentResults[[which(tempCol!=12)]] <- NULL
#	}
	componentResults<-do.call(rbind,componentResults)	
	cat(fileName,"finished deconvolution...\n")
	
	##################################
	#output lib matching results
	###################################
	
	##similarity cutoff=600
	libmatch <- libMatching_Top10(lib,specResults,libMatch_Th)
	write.csv(libmatch,file=paste(params$WorkDir,"/output/decon/",fileName,"_libmatch.csv",sep=""))
	write(strOutput,file=outputPath)
#	write.csv(modelPkList3s,file=paste(params$WorkDir,"/output/decon/",fileName,"_modelPkList3s.csv",sep=""))
#	write.csv(modelPkList2s,file=paste(params$WorkDir,"/output/decon/",fileName,"_modelPkList2s.csv",sep=""))
#	write.csv(modelPkList1s,file=paste(params$WorkDir,"/output/decon/",fileName,"_modelPkList1s.csv",sep=""))
	##############################################################
	#collect spectrum result in dataframe for reverse lib matching
	##############################################################
	
	##exp spectra as lib, get top 1 hit for each std
	#evaluateIdentification(refSpec=lib,inputSpec=specResults,RT_Tolerance=70,minSpecSimilarity=650,fileName,withRT=T,params)
	cat(fileName,"finish matching...\n")
	#write.csv(mdlPkIDVec,file=paste(params$WorkDir,"/output/decon/",fileName,"mdlpklist.csv",sep=""))
	
	################################################
	#collect component result in dataframe for alignment
	###################################################
	
	##remove those model peaks with 0 intensity
	Model_pos <- which(componentResults$isModel==1&componentResults$Intensity==0)
	if (length(Model_pos)>0) {
		for (modelpos in Model_pos) {
			#Remove the old model peak
			componentResults[modelpos,]$isModel=0
			
			#Make the unique as model peak
			compID <- componentResults[modelpos,]$compoundID
			Unique_pos <- which(componentResults$compoundID==compID&componentResults$isUnique==1)
			componentResults[Unique_pos,]$isModel=1
		}
	}
	
	##output and filter out zero-int fragments
	componentResults<-subset(componentResults,Intensity>0,select=c("mz","lboundInd","rboundInd","pkInd","Intensity","area","windowID","compoundID","isModel","isUnique","factor","compID"))
	
#	componentResults$compID <- 0
#	windowVScompound<-subset(componentResults,select=c("windowID","compoundID"))#get common mother Ions which occur in 80% of datasets
#	windowVScompound<-unique(windowVScompound)
#	for (WinComId in 1:nrow(windowVScompound)) {
#		uniqueID <- windowVScompound[WinComId,]
#		componentResults[which(componentResults$windowID==uniqueID$windowID&componentResults$compoundID==uniqueID$compoundID),]$compID=WinComId
#	}
	
	write.csv(componentResults,row.names=FALSE,file=paste(params$WorkDir,"/output/decon/",fileName,"Decon.csv",sep=""))
	cat(fileName,"writing component results...\n")
}


decomposition_temp<-function(inFilePath,params,cl,isDistParallel,clustingType)
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
	
	#To check those potential compounds only
#	RTVector<-read.csv(paste(WorkDir,"output/RTlist.csv",sep=""))$ET
#	TICRTvec<-(TICpeaklist$pkInd-1)*ScanInterval+delaytime
#	winIDs<-NULL		
#	for(i in RTVector)
#	{
#		winIDs<-c(winIDs,which(abs(TICRTvec-i)<=0.2))
#	}
#	winIDs<-sort(unique(winIDs))
	
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
	outputPath<-paste(paste(params$WorkDir,"output/decon/",sep=""),paste(fileName,"_ModelPkProfiles.txt",sep=""),sep="")
	strOutput <- NULL
	
#	modelPkList3s <-NULL
#	modelPkList2s <-NULL
#	modelPkList1s <-NULL
	
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
		localEICPeakList<-subset(EICpeaklist,flag==0&pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		
		if (nrow(localEICPeakList)!= 0) {
			
			if (max(localEICPeakList$Intensity)>=100) {
				##keep the unique mass profile of local EIC list
				##original local EIC list consists of merged peaks
				localEICPeakList <- subset(localEICPeakList,isApex==1&Intensity>=100)
				
				if(nrow(localEICPeakList)>1)
				{
					
					#localEICPeakList$Sig = 0
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
#						if (mean(curProfile$int)!=0) {
#							Sig <- var(curProfile$int)/mean(curProfile$int)
#							localEICPeakList[i,]$Sig <- Sig
#						}
						allprofiles[[i]]<-curProfile
						names(allprofiles)[i]<-rownames(curEICpk)	
						
					}		
					
					#CandidatePeaklist<-subset(localEICPeakList,StN>=300&shrp>=5&Sig>=300&gss<=3)
					CandidatePeaklist<-subset(localEICPeakList,StN>=500&shrp>=5&gss<=2.5)
					############################
					##scoring system
					############################	
					
					if (nrow(CandidatePeaklist)>0) {
						#mass score
						CandidatePeaklist$f1<-scale(CandidatePeaklist$mz,min(CandidatePeaklist$mz),diff(range(CandidatePeaklist$mz)))	
						##gassian similarity score
						CandidatePeaklist$f2<-scale(cos(CandidatePeaklist$gss*pi/180),min(cos(CandidatePeaklist$gss*pi/180)),diff(range(cos(CandidatePeaklist$gss*pi/180))))
						## Signal to Noise
						CandidatePeaklist$f3<-scale(log(CandidatePeaklist$StN),min(log(CandidatePeaklist$StN)),diff(range(log(CandidatePeaklist$StN))))
						CandidatePeaklist$f4<-scale(log(CandidatePeaklist$Intensity),min(log(CandidatePeaklist$Intensity)),diff(range(log(CandidatePeaklist$Intensity))))
						## Singnificance level
						#CandidatePeaklist$f4<-scale(log(CandidatePeaklist$Sig),min(log(CandidatePeaklist$Sig)),diff(range(log(CandidatePeaklist$Sig))))
						CandidatePeaklist$f1[is.na(CandidatePeaklist$f1)]<-0
						CandidatePeaklist$f2[is.na(CandidatePeaklist$f2)]<-0
						CandidatePeaklist$f3[is.na(CandidatePeaklist$f3)]<-0
						CandidatePeaklist$f4[is.na(CandidatePeaklist$f4)]<-0
						##linear combination of 3 scores
						#scores<-(10/7)*(0.1*CandidatePeaklist$f1+0.3*CandidatePeaklist$f2+0.2*CandidatePeaklist$f3+0.1*CandidatePeaklist$f4)
						scores<-(10/8)*(0.1*CandidatePeaklist$f1+0.3*CandidatePeaklist$f2+0.2*CandidatePeaklist$f3+0.2*CandidatePeaklist$f4)
						#scores<-(10/7)*(0.1*CandidatePeaklist$f1+0.3*CandidatePeaklist$f2+0.3*CandidatePeaklist$f3)
						CandidatePeaklist$score<-scale(scores,F,0.001)
						
#						a <- hist(CandidatePeaklist$score)
#						plot(a$breaks[-1],a$counts/sum(a$counts),type="h")
#						
						
						## Create good peaks for clustering
						goodShapePeaklist<-subset(CandidatePeaklist,score>=min(CandidatePeaklist$score)+diff(range(CandidatePeaklist$score))*0.25)
						goodShapePeaklist<-subset(goodShapePeaklist,!mz%in%CurNonUmassVec)
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
									modelPkList<- goodShapePeaklist[which.max(goodShapePeaklist$score),]
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
#						modelPkList3 <- NULL
#						modelPkList2 <- NULL
#						modelPkList1 <- NULL
							##No co-eluting compounds (only one or no model peak detected)
							##Then extract spectrum directly
							if(is.null(modelPkList))
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
#									clustID<-format((maxInd-1)*ScanInterval+delaytime,digit=4)
#									names(specList)[1]<-paste(clustID,curMz)
								
									names(specList)[1]<-paste(maxInd,windowID,1,curMz)
									
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
#								modelPkList1 <- curTICPk$pkInd
									
								}					
								componentList[[1]][which.max(componentList[[1]]$Intensity),]$isUnique=1
								##replace decom as 0, stop while loop
								decom=0

							}else
							{
								## >= one model peaks
								if (nrow(modelPkList)>1) {
									rmmodelPks<-data.frame(row=as.integer(),col=as.integer())						
									modelPkDistance <- as.matrix(modelPkDist(modelPkList,profiles))
									
									for (i in c(1:nrow(modelPkDistance))) {
										subDistance <- modelPkDistance[i,]
										par_name <- rownames(modelPkDistance)[i]
										par_index <- modelPkList[par_name,]$pkInd
										
										for (j in c(1:length(subDistance))) {
											dau_name <- names(subDistance)[j]
											dau_index <- modelPkList[dau_name,]$pkInd
											
											comb_factor <- abs(par_index - dau_index)*modelPkDistance[i,j]
											
											if (comb_factor<= 12||(abs(par_index - dau_index)<=6&modelPkDistance[i,j]<=10)) {
												rmmodelPks <- rbind(rmmodelPks,data.frame(row=i,col=j))
											}	
										}
									}
									
									rmmodelPks <- subset(rmmodelPks,row!=col)
									#######################################
									#Step 2: grouping pairs
									#######################################
									if (nrow(rmmodelPks)>=1) {
										m=c(1:nrow(rmmodelPks))
										l = 1
										modelGroup <- list()
										
										while (length(m)>0) {
											temp <- NULL
											for (i in m[2:length(m)]) {
												if (length(rmmodelPks[i,][(rmmodelPks[i,]%in%rmmodelPks[m[1],])])>=1) {
													temp <- c(temp, i) 
												}
											}
											temp <- c(temp, m[1]) 
											m <- m[!m%in%temp]
											modelGroup[[l]] <- temp
											l <- l+1
										}
										
										###########################################
										#Step 3: get representative for each group
										###########################################
										
										modelPkList_temp <- NULL
										rmmodelPkList_temp <- NULL
										for (i in c(1:length(modelGroup))) {
											rmmodelPkList <- modelPkList[unique(unlist(rmmodelPks[modelGroup[[i]],])),]
											modelPkList_temp <- rbind(modelPkList_temp,rmmodelPkList[which.max(rmmodelPkList$score),])
											rmmodelPkList_temp <- rbind(rmmodelPkList_temp,modelPkList[rownames(modelPkList)%in%rownames(rmmodelPkList),])
											##Only keep the one with highest score
											#modelPkList <- rbind(modelPkList,modelPkList_temp)	
										}
										
										## get removed peaks and update model Peak list
										rmmodelPkList_temp <- rmmodelPkList_temp[!rownames(rmmodelPkList_temp)%in%rownames(modelPkList_temp),]
										modelPkList <- modelPkList[!rownames(modelPkList)%in%rownames(rmmodelPkList_temp),]
										
									}
								}
								
								mdlPkIDs<-rownames(modelPkList)

								lbound<-min(localEICPeakList$lboundInd)
								rbound<-max(localEICPeakList$rboundInd)
								
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
									#curRT<-format((modelPkList$pkInd[i]-1)*ScanInterval+delaytime,digit=4)
									#names(specList)[i]<-paste(curRT,modelPkList$mz[i])
									names(specList)[i]<-paste(modelPkList$pkInd[i],windowID,i,modelPkList$mz[i])
									
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
										
										if (A[m]==0&modelPkList[curmdlPkID,]$mz==curMz) {
											componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,
															rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=10,
															area=10,windowID=windowID,compoundID=j,isModel=1,isUnique=0,factor=A[m]))
										}
										
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
								Num <- nrow(modelPkList)									
								if (Num>1) {
									RepeatModelPk <- NULL									
									indexPairs<-combinations(Num,2)
									for(i in 1:nrow(indexPairs))
									{
										index <- indexPairs[i,]
										index1 <- index[1]
										index2 <- index[2]
										
										##calculate time difference between two model peaks
										modelPKET <- modelPkList$pkInd[c(index1,index2)]	
										
										#fragments with the maximal time difference 10 scans are considered a compound
										if (abs(modelPKET[1]-modelPKET[2]) <=10&nrow(specList[[index1]])!=0&nrow(specList[[index2]])!=0)
										{
											score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
											if (score>=790|0.9*score+(1-(abs(modelPKET[1]-modelPKET[2])/10))*100 >=850) {
												selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index2],mdlPkIDs[index1])
												RepeatModelPk <- c(RepeatModelPk,selectedModelPkID)
											}
										}
									}						
									##select unique filtered model peaks 
									RepeatModelPk <- unique(RepeatModelPk)
									modelPkList <- modelPkList[!rownames(modelPkList)%in%RepeatModelPk,]
								}
										
								if (nrow(modelPkList) == Num) {
									decom = 0
									
									if (nrow(modelPkList) == 1) {
										componentList[[1]][which.max(componentList[[1]]$Intensity),]$isUnique=1
									}else if (nrow(modelPkList) == 2) {
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
									
									## output model pk profiles
									for (model in c(1:nrow(modelPkList))) {
										curModel <- modelPkList[model,]
										curModelInt <- S[,model]
										strOutput <- paste(strOutput,getModelPkProfile(curModel,curModelInt,windowID,lbound,rbound,model),sep="\n\n")
									}
								} # Model Pk keep the same number
							} # model PK num is >2 at first
						} # while statement
						
						##flag analyzed EIC peaks 
						EICpeaklist[rownames(localEICPeakList),]$flag<-1
						componentResults<-c(componentResults,componentList)
						specResults<-c(specResults,specList)
#					modelPkList3s <- c(modelPkList3s,modelPkList3)
#					modelPkList2s <- c(modelPkList2s,modelPkList2)
#					modelPkList1s <- c(modelPkList1s,modelPkList1)
						
						##model peak id
						mdlPkIDVec<-c(mdlPkIDVec,mdlPkIDs)	
					} #candidate checking 				
				} # local Peaks checking3
			} # local peaks checking2	
		} # local Peaks checking1		
	}
	
	componentResults<-do.call(rbind,componentResults)	
	cat(fileName,"finished deconvolution...\n")
	
	##################################
	#output lib matching results
	###################################
	
	##similarity cutoff=600
	libmatch <- libMatching(lib,specResults,600)
	write.csv(libmatch,file=paste(params$WorkDir,"/output/decon/",fileName,"_libmatch.csv",sep=""))
	write(strOutput,file=outputPath)
#	write.csv(modelPkList3s,file=paste(params$WorkDir,"/output/decon/",fileName,"_modelPkList3s.csv",sep=""))
#	write.csv(modelPkList2s,file=paste(params$WorkDir,"/output/decon/",fileName,"_modelPkList2s.csv",sep=""))
#	write.csv(modelPkList1s,file=paste(params$WorkDir,"/output/decon/",fileName,"_modelPkList1s.csv",sep=""))
	##############################################################
	#collect spectrum result in dataframe for reverse lib matching
	##############################################################
	
	##exp spectra as lib, get top 1 hit for each std
	#evaluateIdentification(refSpec=lib,inputSpec=specResults,RT_Tolerance=70,minSpecSimilarity=650,fileName,withRT=T,params)
	cat(fileName,"finish matching...\n")
	#write.csv(mdlPkIDVec,file=paste(params$WorkDir,"/output/decon/",fileName,"mdlpklist.csv",sep=""))
	
	################################################
	#collect component result in dataframe for alignment
	###################################################
	
	##remove those model peaks with 0 intensity
	Model_pos <- which(componentResults$isModel==1&componentResults$Intensity==0)
	if (length(Model_pos)>0) {
		for (modelpos in Model_pos) {
			#Remove the old model peak
			componentResults[modelpos,]$isModel=0
			
			#Make the unique as model peak
			compID <- componentResults[modelpos,]$compoundID
			Unique_pos <- which(componentResults$compoundID==compID&componentResults$isUnique==1)
			componentResults[Unique_pos,]$isModel=1
		}
	}
	
	##output and filter out zero-int fragments
	componentResults<-subset(componentResults,Intensity>0,select=c("mz","lboundInd","rboundInd","pkInd","Intensity","area","windowID","compoundID","isModel","isUnique","factor"))
	
	componentResults$compID <- 0
	windowVScompound<-subset(componentResults,select=c("windowID","compoundID"))#get common mother Ions which occur in 80% of datasets
	windowVScompound<-unique(windowVScompound)
	for (WinComId in 1:nrow(windowVScompound)) {
		uniqueID <- windowVScompound[WinComId,]
		componentResults[which(componentResults$windowID==uniqueID$windowID&componentResults$compoundID==uniqueID$compoundID),]$compID=WinComId
	}
	
	write.csv(componentResults,row.names=FALSE,file=paste(params$WorkDir,"/output/decon/",fileName,"Decon.csv",sep=""))
	cat(fileName,"writing component results...\n")
}

##################################################
# old version of decomposition"
# consider single model peak without decomposition 
##################################################

decomposition_singleModelPk_preserved<-function(inFilePath,params,cl,isDistParallel,clustingType)
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

	#To check those potential compounds only
#	RTVector<-read.csv(paste(WorkDir,"output/RTlist.csv",sep=""))$ET
#	TICRTvec<-(TICpeaklist$pkInd-1)*ScanInterval+delaytime
#	winIDs<-NULL		
#	for(i in RTVector)
#	{
#		winIDs<-c(winIDs,which(abs(TICRTvec-i)<=0.2))
#	}
#	winIDs<-sort(unique(winIDs))
	
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
	outputPath<-paste(paste(params$WorkDir,"output/decon/",sep=""),paste(fileName,"_ModelPkProfiles.txt",sep=""),sep="")
	strOutput <- NULL
	
#	modelPkList3s <-NULL
#	modelPkList2s <-NULL
#	modelPkList1s <-NULL
	
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
		localEICPeakList<-subset(EICpeaklist,flag==0&pkInd>=curTICPk$lboundInd&pkInd<=curTICPk$rboundInd)
		
		if (nrow(localEICPeakList)!= 0) {
			
			if (max(localEICPeakList$Intensity)>=100) {
				##keep the unique mass profile of local EIC list
				##original local EIC list consists of merged peaks
				localEICPeakList <- subset(localEICPeakList,isApex==1&Intensity>=100)
				
				if(nrow(localEICPeakList)>1)
				{
					
					#localEICPeakList$Sig = 0
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
#						if (mean(curProfile$int)!=0) {
#							Sig <- var(curProfile$int)/mean(curProfile$int)
#							localEICPeakList[i,]$Sig <- Sig
#						}
						allprofiles[[i]]<-curProfile
						names(allprofiles)[i]<-rownames(curEICpk)	
						
					}		
					
					#CandidatePeaklist<-subset(localEICPeakList,StN>=300&shrp>=5&Sig>=300&gss<=3)
					CandidatePeaklist<-subset(localEICPeakList,StN>=500&shrp>=5&gss<=2.5)
					############################
					##scoring system
					############################	
					
					if (nrow(CandidatePeaklist)>0) {
						#mass score
						CandidatePeaklist$f1<-scale(CandidatePeaklist$mz,min(CandidatePeaklist$mz),diff(range(CandidatePeaklist$mz)))	
						##gassian similarity score
						CandidatePeaklist$f2<-scale(cos(CandidatePeaklist$gss*pi/180),min(cos(CandidatePeaklist$gss*pi/180)),diff(range(cos(CandidatePeaklist$gss*pi/180))))
						## Signal to Noise
						CandidatePeaklist$f3<-scale(log(CandidatePeaklist$StN),min(log(CandidatePeaklist$StN)),diff(range(log(CandidatePeaklist$StN))))
						CandidatePeaklist$f4<-scale(log(CandidatePeaklist$Intensity),min(log(CandidatePeaklist$Intensity)),diff(range(log(CandidatePeaklist$Intensity))))
						## Singnificance level
						#CandidatePeaklist$f4<-scale(log(CandidatePeaklist$Sig),min(log(CandidatePeaklist$Sig)),diff(range(log(CandidatePeaklist$Sig))))
						CandidatePeaklist$f1[is.na(CandidatePeaklist$f1)]<-0
						CandidatePeaklist$f2[is.na(CandidatePeaklist$f2)]<-0
						CandidatePeaklist$f3[is.na(CandidatePeaklist$f3)]<-0
						CandidatePeaklist$f4[is.na(CandidatePeaklist$f4)]<-0
						##linear combination of 3 scores
						#scores<-(10/7)*(0.1*CandidatePeaklist$f1+0.3*CandidatePeaklist$f2+0.2*CandidatePeaklist$f3+0.1*CandidatePeaklist$f4)
						scores<-(10/8)*(0.1*CandidatePeaklist$f1+0.3*CandidatePeaklist$f2+0.2*CandidatePeaklist$f3+0.2*CandidatePeaklist$f4)
						#scores<-(10/7)*(0.1*CandidatePeaklist$f1+0.3*CandidatePeaklist$f2+0.3*CandidatePeaklist$f3)
						CandidatePeaklist$score<-scale(scores,F,0.001)
						
#						a <- hist(CandidatePeaklist$score)
#						plot(a$breaks[-1],a$counts/sum(a$counts),type="h")
#						
						
						## Create good peaks for clustering
						goodShapePeaklist<-subset(CandidatePeaklist,score>=min(CandidatePeaklist$score)+diff(range(CandidatePeaklist$score))*0.25)
						goodShapePeaklist<-subset(goodShapePeaklist,!mz%in%CurNonUmassVec)
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
									modelPkList<- goodShapePeaklist[which.max(goodShapePeaklist$score),]
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
#						modelPkList3 <- NULL
#						modelPkList2 <- NULL
#						modelPkList1 <- NULL
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
#								modelPkList1 <- curTICPk$pkInd
									
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
#							modelPkList2 <- modelPkList$pkInd
							}else
							{
								##More than one model peaks
								##Decompostion
								
								#######################################
								#Makeup 1 for splitting issue
								#Checking pairwise model peaks distance
								#######################################
								
#							modelPkDistance <- as.matrix(modelPkDist(modelPkList,profiles))
#							rmmodelPk <- NULL
#							ModelPkNames <- rownames(modelPkDistance)
#							
#							##model same matrix as modelPKDistance with index
#							psmatrix <- matrix(c(1:length(modelPkDistance)),ncol=nrow(modelPkDistance),byrow=T)
#							psmodel <- which(modelPkDistance<=4&modelPkDistance>0)
#							
#							##get index of model peaks with smaller distance
#							for (i in psmodel) {
#								pos <- matrix.index(psmatrix, i)
#								rmmodelPk <- c(rmmodelPk,ModelPkNames[pos[1]],ModelPkNames[pos[2]])
#							}
#							rmmodelPks <- unique(rmmodelPk)
#							
#							rmmodelPkList <- modelPkList[rmmodelPks,]
#							modelPkList_temp <- rmmodelPkList[which.max(rmmodelPkList$score),]
#							modelPkList <- modelPkList[!rownames(modelPkList)%in%rmmodelPks,]
#							
#							##Only keep the one with highest score
#							modelPkList <- rbind(modelPkList,modelPkList_temp)				
#							mdlPkIDs<-rownames(modelPkList)
#							
#							##only one model peak after testing
#							##GO back to the upper step
#							if (nrow(modelPkList)==1) decom=1
								
								#######################################
								#Makeup 1 for splitting issue
								#Checking pairwise model peaks distance
								#Step 1: finding pairs
								#######################################
								
								rmmodelPks<-data.frame(row=as.integer(),col=as.integer())						
								modelPkDistance <- as.matrix(modelPkDist(modelPkList,profiles))
								
								for (i in c(1:nrow(modelPkDistance))) {
									subDistance <- modelPkDistance[i,]
									par_name <- rownames(modelPkDistance)[i]
									par_index <- modelPkList[par_name,]$pkInd
									
									for (j in c(1:length(subDistance))) {
										dau_name <- names(subDistance)[j]
										dau_index <- modelPkList[dau_name,]$pkInd
										
										comb_factor <- abs(par_index - dau_index)*modelPkDistance[i,j]
										
										if (comb_factor<= 12||(abs(par_index - dau_index)<=6&modelPkDistance[i,j]<=10)) {
											rmmodelPks <- rbind(rmmodelPks,data.frame(row=i,col=j))
										}	
									}
								}
								
								rmmodelPks <- subset(rmmodelPks,row!=col)
								#######################################
								#Step 2: grouping pairs
								#######################################
								if (nrow(rmmodelPks)>=1) {
									m=c(1:nrow(rmmodelPks))
									l = 1
									modelGroup <- list()
									
									while (length(m)>0) {
										temp <- NULL
										for (i in m[2:length(m)]) {
											if (length(rmmodelPks[i,][(rmmodelPks[i,]%in%rmmodelPks[m[1],])])>=1) {
												temp <- c(temp, i) 
											}
										}
										temp <- c(temp, m[1]) 
										m <- m[!m%in%temp]
										modelGroup[[l]] <- temp
										l <- l+1
									}
									
									###########################################
									#Step 3: get representative for each group
									###########################################
									
									modelPkList_temp <- NULL
									rmmodelPkList_temp <- NULL
									for (i in c(1:length(modelGroup))) {
										rmmodelPkList <- modelPkList[unique(unlist(rmmodelPks[modelGroup[[i]],])),]
										modelPkList_temp <- rbind(modelPkList_temp,rmmodelPkList[which.max(rmmodelPkList$score),])
										rmmodelPkList_temp <- rbind(rmmodelPkList_temp,modelPkList[rownames(modelPkList)%in%rownames(rmmodelPkList),])
										##Only keep the one with highest score
										#modelPkList <- rbind(modelPkList,modelPkList_temp)	
									}
									
									## get removed peaks and update model Peak list
									rmmodelPkList_temp <- rmmodelPkList_temp[!rownames(rmmodelPkList_temp)%in%rownames(modelPkList_temp),]
									modelPkList <- modelPkList[!rownames(modelPkList)%in%rownames(rmmodelPkList_temp),]
									
								}
								##only one model peak after testing
								##GO back to the upper step
								mdlPkIDs<-rownames(modelPkList)
								if (nrow(modelPkList)==1) decom=1
								
								else {
									##still multiple models
									lbound<-min(localEICPeakList$lboundInd)
									rbound<-max(localEICPeakList$rboundInd)
																	
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
											
											if (A[m]==0&modelPkList[curmdlPkID,]$mz==curMz) {
												componentList[[j]]<-rbind(componentList[[j]],data.frame(mz=as.integer(curMz),lboundInd=modelPkList[curmdlPkID,]$lboundInd,
																rboundInd=modelPkList[curmdlPkID,]$rboundInd,pkInd=modelPkList[curmdlPkID,]$pkInd,Intensity=10,
																area=10,windowID=windowID,compoundID=j,isModel=1,isUnique=0,factor=A[m]))
											}
											
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
										
										#score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
#										if (score >=970||(abs(modelPKET[1]-modelPKET[2]) <=6&&score >=900))
#										{
#											selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index2],mdlPkIDs[index1])
#											RepeatModelPk <- c(RepeatModelPk,selectedModelPkID)
#										}
										
										
										#fragments with the maximal time difference 10 scans are considered a compound
										if (abs(modelPKET[1]-modelPKET[2]) <=10&nrow(specList[[index1]])!=0&nrow(specList[[index2]])!=0)
										{
											score <- specDistCal(specList[[index1]],specList[[index2]],isWeight=T,isNist=T)
											if (0.9*score+(1-(abs(modelPKET[1]-modelPKET[2])/10))*100 >=850) {
												selectedModelPkID <- ifelse(modelPkList$Intensity[index1]>=modelPkList$Intensity[index2],mdlPkIDs[index2],mdlPkIDs[index1])
												RepeatModelPk <- c(RepeatModelPk,selectedModelPkID)
											}
										}
									}
									
									##select unique filtered model peaks 
									RepeatModelPk <- unique(RepeatModelPk)
									modelPkList <- modelPkList[!rownames(modelPkList)%in%RepeatModelPk,]
#								modelPkList3 <- modelPkList$pkInd
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
									
									for (model in c(1:nrow(modelPkList))) {
										curModel <- modelPkList[model,]
										curModelInt <- S[,model]
										strOutput <- paste(strOutput,getModelPkProfile(curModel,curModelInt,windowID,lbound,rbound,model),sep="\n\n")
									}
									
								} # else model PK num is >2 after testing distance	
							} # model PK num is >2 at first
						} # while statement
						
						##flag analyzed EIC peaks 
						EICpeaklist[rownames(localEICPeakList),]$flag<-1
						componentResults<-c(componentResults,componentList)
						specResults<-c(specResults,specList)
#					modelPkList3s <- c(modelPkList3s,modelPkList3)
#					modelPkList2s <- c(modelPkList2s,modelPkList2)
#					modelPkList1s <- c(modelPkList1s,modelPkList1)
						
						##model peak id
						mdlPkIDVec<-c(mdlPkIDVec,mdlPkIDs)	
					} #candidate checking 				
				} # local Peaks checking3
			} # local peaks checking2	
		} # local Peaks checking1		
	}
	
	componentResults<-do.call(rbind,componentResults)	
	cat(fileName,"finished deconvolution...\n")
	
	##################################
	#output lib matching results
	###################################
	
	##similarity cutoff=600
	libmatch <- libMatching(lib,specResults,600)
	write.csv(libmatch,file=paste(params$WorkDir,"/output/decon/",fileName,"_libmatch.csv",sep=""))
	write(strOutput,file=outputPath)
#	write.csv(modelPkList3s,file=paste(params$WorkDir,"/output/decon/",fileName,"_modelPkList3s.csv",sep=""))
#	write.csv(modelPkList2s,file=paste(params$WorkDir,"/output/decon/",fileName,"_modelPkList2s.csv",sep=""))
#	write.csv(modelPkList1s,file=paste(params$WorkDir,"/output/decon/",fileName,"_modelPkList1s.csv",sep=""))
	##############################################################
	#collect spectrum result in dataframe for reverse lib matching
	##############################################################
	
	##exp spectra as lib, get top 1 hit for each std
	#evaluateIdentification(refSpec=lib,inputSpec=specResults,RT_Tolerance=70,minSpecSimilarity=650,fileName,withRT=T,params)
	cat(fileName,"finish matching...\n")
	#write.csv(mdlPkIDVec,file=paste(params$WorkDir,"/output/decon/",fileName,"mdlpklist.csv",sep=""))
	
	################################################
	#collect component result in dataframe for alignment
	###################################################
	
	##remove those model peaks with 0 intensity
	Model_pos <- which(componentResults$isModel==1&componentResults$Intensity==0)
	if (length(Model_pos)>0) {
		for (modelpos in Model_pos) {
			#Remove the old model peak
			componentResults[modelpos,]$isModel=0
			
			#Make the unique as model peak
			compID <- componentResults[modelpos,]$compoundID
			Unique_pos <- which(componentResults$compoundID==compID&componentResults$isUnique==1)
			componentResults[Unique_pos,]$isModel=1
		}
	}

	##output and filter out zero-int fragments
	componentResults<-subset(componentResults,Intensity>0,select=c("mz","lboundInd","rboundInd","pkInd","Intensity","area","windowID","compoundID","isModel","isUnique","factor"))
	write.csv(componentResults,row.names=FALSE,file=paste(params$WorkDir,"/output/decon/",fileName,"Decon.csv",sep=""))
	cat(fileName,"writing component results...\n")
}

fr <- function(x,M,S) {
	m<-0
	for(i in 1:ncol(S))
	{
		m<-m+x[i]*S[,i]
	}
	sum((M-m)^2)
	
}

#the main routine to call all the steps
parGcProprocess <-function(params)
{	
	#get the first scan and scan interval from *_RT.csv file which is produced by SPlitEIC in C code
	fileInfo<-read.csv(paste(params$WorkDir,"/output/",params$JobName,"_RT.csv",sep=""))
	params$delaytime<-fileInfo$firstRT[1]
	params$ScanInterval<-fileInfo$ScanInterval[1]
	#denoise TIC before TIC peakpicking
#	parTICDenoising(params) 
	#TIC peak picking and boundary detection in order to set the preliminary decon window	
	parTICPeakpicking(params,denoised=F)
	#denoising EIC before EIC peak picking
	parDenoising(params,codeDir="code")
	
	#EIC peak picking and boundary detection in order to get EIC profile for each fragment ion
	parEICpeakpicking(params)
	#set the exluded ion list
	params$NonUmassVec<-c(1:50,73,147,221)
	#deconvolute the unresolved the peaks
	deconvolution(params,denoised=F)
	
}

##change the file name in batch manner
renameFiles<-function()
{
	filelist<-list.files(full.names = F, paste(params$WorkDir, "/raw",sep=""))
	for(i in 1:length(filelist))
	{
		curFile<-filelist[i]
		
		fileName<-unlist(strsplit(curFile,split="\\."))[1]
		EICFile<-paste(fileName,"EIC.cdf",sep="")
		
		#######remove denoised_ prefix
		file.rename(from=paste(params$WorkDir, "/output/EIC/","denoised_",EICFile,sep=""),to=paste(params$WorkDir, "/output/EIC/",EICFile,sep=""))	
		#########rename .CDF to .cdf
#		file.rename(from=paste(params$WorkDir, "/raw/",curFile,sep=""),to=paste(params$WorkDir, "/raw/",fileName,".cdf",sep=""))
	}
}

wavCWT_Modify<-function (x, scale.range = deltat(x) * c(1, length(x)), n.scale = 100, 
		wavelet = "gaussian2", shift = 5, variance = 1) 
{
	checkVectorType(scale.range, "numeric")
	checkScalarType(n.scale, "integer")
	checkScalarType(wavelet, "character")
	checkScalarType(shift, "numeric")
	checkScalarType(variance, "numeric")
	checkRange(n.scale, c(1, Inf))
	series.name <- deparse(substitute(x))
	if (length(scale.range) != 2) 
		stop("scale.range must be a two-element numeric vector")
	if (variance <= 0) 
		stop("variance input must be positive")
	sampling.interval <- deltat(x)
	octave <- logb(scale.range, 2)
	scale <- ifelse1(n.scale > 1, 2^c(octave[1] + seq(0, n.scale - 
									2) * diff(octave)/(floor(n.scale) - 1), octave[2]), scale.range[1])
	scale <- unique(round(scale/sampling.interval) * sampling.interval)
	n.scale <- length(scale)
	if (abs(min(scale) - sampling.interval) > .Machine$double.eps) 
		stop("Minimum scale must be greater than or equal to sampling interval ", 
				"of the time series")
	if (inherits(x, "signalSeries")) 
		times <- as(x@positions, "numeric")
	else times <- time(x)
	x <- as.vector(x)
	storage.mode(x) <- "double"
	gauss1 <- c("gaussian1", "gauss1")
	gauss2 <- c("gaussian2", "gauss2", "mexican hat", "sombrero")
	supported.wavelets <- c("haar", gauss1, gauss2, "morlet")
	wavelet <- match.arg(lowerCase(wavelet), supported.wavelets)
	filter <- mutilsFilterTypeContinuous(wavelet)
	if (filter == 4) {
		filter.arg <- sqrt(variance)
		wavelet <- "gaussian1"
	}
	else if (filter == 5) {
		filter.arg <- sqrt(variance)
		wavelet <- "gaussian2"
	}
	else if (filter == 6) {
		filter.arg <- shift
		wavelet <- "morlet"
	}
	else if (filter == 7) {
		filter.arg <- 0
		wavelet <- "haar"
		scale <- sampling.interval * unique(round(scale/sampling.interval))
	}
	else stop("Unsupported filter type")
	z <- .Call("RS_wavelets_transform_continuous_wavelet", as.numeric(x), 
			as.numeric(sampling.interval), as.integer(filter), as.numeric(filter.arg), 
			as.numeric(scale), COPY = rep(FALSE, 5), CLASSES = c("numeric", 
					"numeric", "integer", "numeric", "numeric"), PACKAGE = "ifultools")
	if (wavelet != "morlet") 
		z <- Re(z)
	attr(z, "scale") <- scale
	attr(z, "time") <- as.vector(times)
	attr(z, "wavelet") <- wavelet
	attr(z, "series") <- x
	attr(z, "sampling.interval") <- sampling.interval
	attr(z, "series.name") <- series.name
	attr(z, "n.sample") <- length(x)
	attr(z, "n.scale") <- n.scale
	attr(z, "filter.arg") <- filter.arg
	oldClass(z) <- "wavCWT"
	z
}

wavCWTTree_Modify<-function (x, n.octave.min = 1, tolerance = 0, type = "maxima") 
{
	
	"WTMM" <- function(x, tolerance = NULL, type = "maxima") {
		if (!is(x, "wavCWT")) 
			stop("Input object must be of class wavCWT")
		x.attr <- attributes(x)
		times <- x.attr$time
		scales <- x.attr$scale
		n.sample <- x.attr$n.sample
		series <- x.attr$series
		if (is.null(tolerance)) {
			tolerance <- mad(Mod(x[, 1]))/scales
		}
		if (length(tolerance) < length(scales)) 
			tolerance <- tolerance[1]/sqrt(scales)
		wtmmz <- .Call("RS_wavelets_transform_continuous_wavelet_modulus_maxima", 
				as.matrix(x) + (0+0i), tolerance, mutilsTransformPeakType(type), 
				CLASSES = c("matrix", "numeric", "integer"), COPY = rep(FALSE, 
						3), PACKAGE = "ifultools")
		z <- matrix(0, nrow = nrow(x), ncol = ncol(x))
		z[matrix(unlist(wtmmz), ncol = 2) + 1] <- 1
		z
	}# WTMM func
	
	"wtmmBranches" <- function(wtmm, extrema.mask, times, scales, 
			span.min = 5, gap.max = 3, skip = NULL, sampling.interval = 1) {
		scales <- as.integer(scales/sampling.interval)
		n.scale <- ncol(extrema.mask)
		n.sample <- nrow(extrema.mask)
		if (is.null(scales)) 
			scales <- 1:n.scale
		iwtmm <- which(extrema.mask[, n.scale] > 0)		
		if (length(iwtmm)!=0) {
			
			iscale <- seq(n.scale - 1, 1, -1)
			tree <- as.list(iwtmm)
			names(tree) <- iwtmm
			peakStatus <- as.list(rep(0, length(iwtmm)))
			names(peakStatus) <- iwtmm
			orphanRidgeList <- NULL
			orphanRidgeName <- NULL
			n.level <- length(iscale)
			for (j in seq(n.level)) {
				iscale.j <- iscale[j]
				scale.j <- scales[iscale.j]
				if (length(iwtmm) == 0) {
					iwtmm <- which(extrema.mask[, iscale.j] > 0)
					next
				}
				span <- scale.j * 2 + 1
				if (span < span.min) 
					span <- span.min
				remove.j <- selPeak.j <- NULL
				for (k in seq(along = iwtmm)) {
					itime <- iwtmm[k]
					itime.start <- itime - span
					if (itime.start < 1) 
						itime.start <- 1
					itime.end <- itime + span
					if (itime.end > n.sample) 
						itime.end <- n.sample
					itime.candidates <- which(extrema.mask[itime.start:itime.end, 
									iscale.j] > 0) + itime.start - 1
					if (length(itime.candidates) == 0) {
						status.k <- peakStatus[[as.character(itime)]]
						
						if (length(status.k)>0&length(scale.j)>0) {
							if (status.k > gap.max & scale.j >= 2) {
								temp <- tree[[as.character(itime)]]
								orphanRidgeList <- c(orphanRidgeList, list(temp[1:(length(temp) - 
																	status.k)]))
								orphanRidgeName <- c(orphanRidgeName, paste(iscale.j + 
														status.k + 1, itime, sep = "_"))
								remove.j <- c(remove.j, as.character(itime))
								next
							}
							else {
								itime.candidates <- itime
								peakStatus[[as.character(itime)]] <- status.k + 
										1
							}
						}
					}
					else {
						peakStatus[[as.character(itime)]] <- 0
						if (length(itime.candidates) >= 2) 
							itime.candidates <- itime.candidates[which.min(abs(itime.candidates - 
															itime))]
					}
					tree[[as.character(itime)]] <- c(tree[[as.character(itime)]], 
							itime.candidates)
					selPeak.j <- c(selPeak.j, itime.candidates)
				}
				if (length(remove.j) > 0) {
					bad.tree <- which(is.element(names(tree), remove.j))
					tree <- tree[-bad.tree]
					peakStatus <- peakStatus[-bad.tree]
				}
				dupPeak.j <- unique(selPeak.j[duplicated(selPeak.j)])
				if (length(dupPeak.j) > 0) {
					bad.tree <- NULL
					for (dupPeak.jk in dupPeak.j) {
						selInd <- which(selPeak.j == dupPeak.jk)
						selLen <- sapply(tree[selInd], length)
						bad.tree.jk <- which.max(selLen)
						bad.tree <- c(bad.tree, selInd[-bad.tree.jk])
						orphanRidgeList <- c(orphanRidgeList, tree[bad.tree.jk])
						orphanRidgeName <- c(orphanRidgeName, paste(iscale.j, 
										selPeak.j[bad.tree.jk], sep = "_"))
					}
					selPeak.j <- selPeak.j[-bad.tree]
					tree <- tree[-bad.tree]
					peakStatus <- peakStatus[-bad.tree]
				}
				names(tree) <- selPeak.j
				names(peakStatus) <- selPeak.j
				if (scale.j >= 2) {
					maxInd.next <- which(extrema.mask[, iscale.j] > 
									0)
					unSelPeak.j <- maxInd.next[!is.element(maxInd.next, 
									selPeak.j)]
					newPeak.j <- as.list(unSelPeak.j)
					names(newPeak.j) <- unSelPeak.j
					tree <- c(tree, newPeak.j)
					iwtmm <- c(selPeak.j, unSelPeak.j)
					newPeakStatus <- as.list(rep(0, length(newPeak.j)))
					names(newPeakStatus) <- newPeak.j
					peakStatus <- c(peakStatus, newPeakStatus)
				}
				else {
					iwtmm <- selPeak.j
				}
			}
			
			if (length(tree)!=0) {
				names(tree) <- paste(1, names(tree), sep = "_")
				names(orphanRidgeList) <- orphanRidgeName
				tree <- c(tree, orphanRidgeList)
				tree <- lapply(tree, rev)
				tree <- tree[unique(names(tree))]
				tree <- lapply(seq(along = tree), function(i, tree, iscale.min, 
								times, scales, wtmm) {
							itime <- tree[[i]]
							iscale <- seq(iscale.min[i], length = length(itime))
							list(itime = itime, iscale = iscale, time = times[itime], 
									scale = scales[iscale], extrema = wtmm[cbind(itime, 
													iscale)])
						}, tree = tree, iscale.min = as.integer(gsub("_.*", "", 
										names(tree))), times = times, scales = scales * sampling.interval, 
						wtmm = wtmm)
				iflat <- lapply(tree, function(x, nr) (x$iscale - 1) * 
									nr + x$itime, nr = nrow(wtmm))
				flatset <- iflat[[1]]
				bad <- NULL
				if(length(iflat)>1)
				{
					for (i in seq(2, length(iflat))) {
						if (any(is.element(iflat[[i]], flatset))) 
							bad <- c(bad, i)
						else flatset <- c(flatset, iflat[[i]])
					}
				}
				
				if (length(bad) > 0) 
					tree <- tree[-bad]
				tree
			}# tree length checking		
		}
	}# wtmmBranches function
	
	
	x.attr <- attributes(x)
	times <- x.attr$time
	scales <- x.attr$scale
	n.sample <- x.attr$n.sample
	sampling.interval <- x.attr$sampling.interval
	border.times <- range(times) + sampling.interval * c(1, -1)
	extrema.mask <- WTMM(x, tolerance = tolerance, type = type)
	if (!identical(dim(x), dim(extrema.mask))) 
		stop("Input WTMM dimenions do not match those of the input CWT matrix")
	
	z <- wtmmBranches(ifelse1(is.complex(x), Mod(as.matrix(x)), 
					as.matrix(x)), extrema.mask, times, scales, sampling.interval = sampling.interval)
	
	if (!is.null(z)) {
		n.scale <- length(scales)
		n.octave <- log2(max(scales)/min(scales))
		n.voice <- (n.scale - 1)/n.octave
		n.scale.min <- as.integer(n.voice * n.octave.min)
		
		good <- which(unlist(lapply(z, function(x, n.scale.min) length(x[[1]]) > 
											n.scale.min, n.scale.min = n.scale.min)))
		if (length(good)!=0) {
			z <- z[good]
			endtime <- unlist(lapply(z, function(x, iscale) x$itime[iscale], 
							iscale = which.min(scales)))
			isort <- order(endtime)
			z <- z[isort]
			names(z) <- seq(z)
			attr(z, "iendtime") <- endtime[isort]
			attr(z, "endtime") <- times[endtime[isort]]
			attr(z, "time") <- times
			attr(z, "scale") <- scales
			attr(z, "extrema.mask") <- extrema.mask
			attr(z, "noise") <- x[, 1]
			attr(z, "branch.hist") <- colSums(extrema.mask * abs(x))
			attr(z, "wavelet") <- attr(x, "wavelet")
			attr(z, "filter.arg") <- attr(x, "filter.arg")
			attr(z, "series.name") <- attr(x, "series.name")
			attr(z, "series") <- attr(x, "series")
			attr(z, "sampling.interval") <- attr(x, "sampling.interval")
			oldClass(z) <- "wavCWTTree"
			z
		}
	}
}


wavCWTPeaks_Modify <- function( x, snr.min = 3, scale.range = NULL, length.min = 10, 
		noise.span = NULL, noise.fun = "quantile", noise.min = NULL)
{
	
	if (!is(x, "wavCWTTree")) 
		stop("Input must be an object of class wavCWTTree")
	xatt <- attributes(x)
	endtimes <- attr(x, "endtime")
	times <- attr(x, "time")
	scale <- attr(x, "scale")
	noise <- attr(x, "noise")
	wavelet <- attr(x, "wavelet")
	series <- attr(x, "series")
	branch.hist <- attr(x, "branch.hist")
	sampling.interval <- abs(diff(times[1:2]))
	if (!is.element(wavelet, "gaussian2")) 
		stop("Only CWT developed using the Mexican hat (gaussian2) filter are supported")
	if (is.null(noise.min)) 
		noise.min <- quantile(abs(attr(x, "noise")), prob = 0.05)
	if (is.null(scale.range)) 
		scale.range <- scale[range(which(branch.hist > quantile(branch.hist, 
										prob = 0.8)))]
	if (is.null(noise.span)) 
		noise.span <- max(0.01 * diff(range(times)), 5 * sampling.interval)
	noise.levels <- unlist(lapply(endtimes, function(x, noise.fun, 
							times, times.range, noise, noise.min, noise.span) {
						time.start <- x - noise.span
						if (time.start < times.range[1]) 
							time.start <- times.range[1]
						time.end <- x + noise.span
						if (time.end < times.range[2]) 
							time.end <- times.range[2]
						ix <- which(times >= time.start & times <= time.end)
						noise.local <- noise.fun(abs(noise[ix]))
						if (noise.local < noise.min) 
							noise.local <- noise.min
						noise.local
					}, noise.fun = switch(noise.fun, quantile = function(x) {
								quantile(x, probs = 0.95)
							}, sd = sd, mad = function(x) {
								mad(x, center = 0)
							}), times = times, times.range = range(times), noise = noise, 
					noise.min = noise.min, noise.span = noise.span))
	tmpargs <- lapply(x, function(x) unlist(lapply(x, function(x, 
										imax) x[imax], imax = which.max(x$extrema))))
	peaks <- data.frame(do.call("rbind", tmpargs))
	peaks <- cbind(data.frame(branch = row.names(peaks)), peaks, 
			data.frame(iendtime = attr(x, "iendtime")))
	peak.snr <- round(peaks[["extrema"]]/noise.levels)
	peak.scale <- peaks[["scale"]]
	branch.lengths <- unlist(lapply(x, function(x, scale.range) length(which(x$scale >= 
												scale.range[1] & x$scale <= scale.range[2])), scale.range = scale.range))
	#good.snr <- peak.snr >= snr.min
	good.scale <- peak.scale >= scale.range[1]      
#	good.scale <- ((peak.scale >= max(scale.range[1],2))& (peak.scale <= min(scale.range[2], 25)))   #Modifed by Wenchao Zhang, Use lower scale boundary and Upper scale boundary to filter  
	#good.length <- branch.lengths >= length.min
	iendtime.min <- max(as.integer(noise.span/sampling.interval/4), 3)
	iendtime.max <- length(times) - iendtime.min + 1
	good.end <- peaks[["iendtime"]] > iendtime.min & peaks[["iendtime"]] < iendtime.max
	
	# Added some codes to filter out some very very small peaks. Added by wenchao zhang
	#peak_extrema <- peaks[["extrema"]]
	
	#good.peak_extrema <- peak_extrema> max(peak_extrema)/10.0
	
	#peaks <- peaks[which(good.snr & good.scale & good.length & good.end & good.peak_extrema), ]	
#	peaks <- peaks[which(good.scale  & good.end), ]	
#	snr <- peak.snr[which(good.scale & good.end)]
	snr <- peak.snr
	
	#Debug for none good peak cases, 
	if(nrow(peaks)>0) row.names(peaks) <- as.character(seq(nrow(peaks)))
	z <- list(x = times[peaks$iendtime], y = series[peaks$iendtime])
	attr(z, "peaks") <- peaks
	#attr(z, "snr.min") <- snr.min
	attr(z, "scale.range") <- scale.range
	attr(z, "length.min") <- length.min
	attr(z, "noise.span") <- noise.span
	attr(z, "noise.fun") <- noise.fun
	attr(z, "noise.min") <- noise.min
	attr(z, "snr") <- snr
	z
}


