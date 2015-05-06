vizTIC<- function(inFilePath) {
	
	
	ncid <- open.ncdf(inFilePath);
    
	# get data from ncdf
	scan_acquisition_time <- get.var.ncdf(ncid, varid="scan_acquisition_time")
	total_intensity <- get.var.ncdf(ncid, varid="total_intensity")
	
	# chromatogram
	fileName<-parseFileName(inFilePath)
	playwith({
	plot(0,type="n", main=fileName, xlab="Retension Time", ylab="intensity", xlim=c(min(scan_acquisition_time)/60, max(scan_acquisition_time)/60) ,ylim=c(0,100))
	points(scan_acquisition_time / 60, total_intensity * 100 / max(total_intensity), type="l",col="blue")
	})
}

#viz EIC peak picking results

vizEICPeak<-function(strMz,curFeatureset,ncid,delaytime,ScanInterval)
{	
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	vectorMz<-get.var.ncdf(ncid, varid="mzVec")
	totalscan<-length(vecInt)/length(vectorMz)
	vecET<-((1:totalscan)-1)*ScanInterval+delaytime
	
	mzVec<-strsplit(strMz,",")[[1]]
	nMz<-length(mzVec)

	nColor<-8

	
	#if no mz specified, then all mz are selected within the windows
	if(nMz==0)
	{
		mzVec<-unique(curFeatureset$curMz)
		nMz<-length(mzVec)
	}else
	{

		curFeatureset<-subset(curFeatureset,curMz%in%mzVec)
	}
	
	minET<-min((curFeatureset$lboundInd-1)*ScanInterval+delaytime)
	maxET<-	max((curFeatureset$rboundInd-1)*ScanInterval+delaytime)
	maxInt<-max(curFeatureset$Intensity)
	
	#get RT vec
	curFeatureset$RT<-(curFeatureset$pkInd-1)*ScanInterval+delaytime
	#normalize int
	curFeatureset$Intensity<-100*curFeatureset$Intensity/maxInt
	datapoints<-subset(curFeatureset,select=c("RT","Intensity"))
	datalabels<-curFeatureset$curMz
	
	
	
	playwith(
	{
		plot(x=0,type="n",main="peak picking results",xlim=c(minET,maxET),ylim=c(0,110),xlab="Retention Time(minutes)",ylab="Intensity",cex=0.03)
		xy<-par("usr")
		xlength<-xy[2]-xy[1] 
		ylength<-xy[4]-xy[3] 
		for(i in 1:nMz)
		{
			mz<-mzVec[i]
			#get EIC
			mzInd<-which(vectorMz==mz)
			startInd<-(mzInd-1)*totalscan+1
			endInd<-startInd+totalscan-1
			curVecInt <- vecInt[startInd:endInd]
			EICind<-which(vecET>=minET&vecET<=maxET)
			#get peak apex
			curFeatures<-subset(curFeatureset,curMz==as.numeric(mz))
			#plot mz legend
			mtext(las = 1,at=i*ylength/30,side=4,text=mz,col=rainbow(nColor)[i%%nColor])
			#plot EIC
			
			points(vecET[EICind],100*curVecInt[EICind]/maxInt,type="l",col=rainbow(nColor)[i%%nColor])
			#plot peak apex
			points(curFeatures$RT,curFeatures$Intensity, col=rainbow(nColor)[i%%nColor],cex=0.7)
			
		}
	},data.points=datapoints,labels=datalabels,new=TRUE)
}


	
#viz deconvolution results by 
vizDecon<-function(strMz,strWindowID,strComponentID,LocaldeconFeatures,ncid,delaytime,ScanInterval)
{	
	vecInt <- get.var.ncdf(ncid, varid="intVec")
	vectorMz<-get.var.ncdf(ncid, varid="mzVec")
	totalscan<-length(vecInt)/length(vectorMz)
	vecET<-((1:totalscan)-1)*ScanInterval+delaytime
	#filter by mz
	mzVec<-unlist(strsplit(strMz,","))
	nMz<-length(mzVec)
	ifelse(nMz>0,curFeatureset<-subset(LocaldeconFeatures,mz%in%mzVec),curFeatureset<-LocaldeconFeatures)
	
	#filter by windowID
	windowVec<-unlist(strsplit(strWindowID,","))
	nWindow<-length(windowVec)
	if(nWindow>0)
	{
		curFeatureset<-subset(curFeatureset,windowID%in%windowVec)
		#filter by compoundID
		compVec<-unlist(strsplit(strComponentID,","))
		nComp<-length(compVec)
		if(nComp>0)
			curFeatureset<-subset(curFeatureset,compoundID%in%compVec)
	
	}
	
	#if no mz specified, then all mz are selected within the windows
	if(nMz==0)
	{
		mzVec<-unique(curFeatureset$mz)
		nMz<-length(mzVec)
	}
	
	curFeatureset$key<-paste(curFeatureset$windowID,curFeatureset$compoundID)
	
	#get the boundary of EIC
	minLbound<-min(curFeatureset$lboundInd)
	maxRbound<-max(curFeatureset$rboundInd)
	minET<-(minLbound-1)*ScanInterval+delaytime
	maxET<-	(maxRbound-1)*ScanInterval+delaytime
	maxInt<-max(curFeatureset$Intensity)
	
	#get RT vec
	curFeatureset$RT<-(curFeatureset$pkInd-1)*ScanInterval+delaytime
	#normalize int
	curFeatureset$Intensity<-100*curFeatureset$Intensity/maxInt
	componentFeatures<-subset(curFeatureset,compoundID>0)
	#get data points
	datapoints<-subset(componentFeatures,select=c("RT","Intensity"))
	datalabels<-paste("window:",componentFeatures$windowID," component:",componentFeatures$compoundID," mz:",componentFeatures$mz," scan:",componentFeatures$pkInd,sep="")
	
	
	playwith(
	{
		plot(x=0,type="n",main="deconvolution results",xlim=c(minET,maxET),ylim=c(0,110),xlab="Retention Time(minutes)",ylab="Intensity",cex=0.03)
		xy<-par("usr")
		xlength<-xy[2]-xy[1] 
		ylength<-xy[4]-xy[3] 
		#plot EIC
		for(i in 1:nMz)
		{
			mz<-mzVec[i]
			#get EIC
			#vecInt <- get.var.ncdf(ncid, varid=as.character(mz))
			#vecET<-((1:length(vecInt))-1)*ScanInterval+delaytime
			mzInd<-which(vectorMz==mz)
			startInd<-(mzInd-1)*totalscan+1
			endInd<-startInd+totalscan-1
			curVecInt <- vecInt[startInd:endInd]
			
			#EIC<-data.frame(RT=vecET,int=curVecInt)
#			EIC<-subset(EIC,RT>=minET&RT<=maxET)
			EICind<-which(vecET>=minET&vecET<=maxET)
			#get peak apex
			#curFeatures<-subset(curFeatureset,mz==as.numeric(mz))
			#plot mz legend
			mtext(las = 1,at=i*ylength/30,side=4,text=mz,col=rainbow(nMz+1)[i+1])
			#plot EIC
			points(vecET[EICind],100*curVecInt[EICind]/maxInt,type="l",col=rainbow(nMz+1)[i+1])

			
		}
		#plot deconvoluted peak apex
		pchVec<-20:30
		npchVec<-length(pchVec)

		
		ComponentID<-unique(componentFeatures$key)
		nComponent<-length(ComponentID)
		for(i in 1:nComponent)
		{

			curComponentID<-ComponentID[i]
			fInd<-which(componentFeatures$key==curComponentID)
			#curFeatures<-subset(componentFeatures,key==curComponentID)			
			points(componentFeatures[fInd,]$RT,componentFeatures[fInd,]$Intensity,pch=pchVec[i%%npchVec],cex=0.7)#col=rainbow(5)[i%%5]
		}

		#plot denoised peak apex
		
		#DenoisedFeatureset<-subset(curFeatureset,compoundID<0)
#		points(DenoisedFeatureset$RT,DenoisedFeatureset$Intensity,pch=4,col="red",cex=0.7)#col=rainbow(5)[i%%5]
	}
	,data.points=datapoints,labels=datalabels,new=TRUE)
	
}

#viz component spectrum
vizSpec<-function(curWindowID,curCompoundID,deconFeatures,delaytime,ScanInterval)
{	

	curComponent<-subset(deconFeatures,windowID==curWindowID&compoundID==curCompoundID)
	RT<-mean((curComponent$pkInd)-1)*ScanInterval+delaytime
	
	colnames(curComponent)[which(colnames(curComponent)=="Intensity")]<-"int"
	componentName<-paste("win",curWindowID," component",curCompoundID,sep="")
	strOutput<-getMSPBlock(curComponent,RT*60,componentName,dbind=1)
	cat(strOutput,"\n")

	spec<-subset(curComponent,select=c("mz","int"))
	
	spec$int<-100*spec$int/max(spec$int)

	playwith(
	{
		plot(x=0,type="n",xlim=c(min(spec$mz),max(spec$mz)),ylim=c(1,100),main="component spectrum",xlab="m/z",ylab="Intensity",cex=0.03)
		points(spec,type="h")
	},data.points=spec,labels=spec$mz,new=TRUE)
}

vizLinearQuantitation<-function(FileIDVsConcertration,alignedFeature,selectedCompoundRow,strMzs,quantitationType)
{
	#plot linear quatitation
		QMassVec<-strsplit(strMzs,",")[[1]]
		measure<-ifelse(quantitationType=="Height","int","pkArea")
		curAlignedFeature<-subset(alignedFeature,mz%in%QMassVec&clustID==selectedCompoundRow$clustID);
		
		curClustID<-selectedCompoundRow$clustID
		curCompName<-selectedCompoundRow$Name
		RT<-selectedCompoundRow$RT

		#curUmass<-subset(curAlignedFeature,clustID==curClustID&mz==curQmz)
		#select Umass
		nrows<-ceiling(sqrt(length(QMassVec)))
		ncols<-ceiling(length(QMassVec)/nrows)

		#x11()
		par(mfrow=c(nrows,ncols))
		for(i in 1:length(QMassVec))
		{
			curQmz<-QMassVec[i]
			curUmass<-subset(curAlignedFeature,mz==curQmz)
			curUmass<-merge(curUmass,FileIDVsConcertration,by.y="fileid",by.x="FileID")
			curUmass<-subset(curUmass,select=c("concerntration",measure))
			LinearRegression<-lm(curUmass$int ~ curUmass$concerntration)
			
			plot(curUmass,type="p",main=paste(curCompName," ",round(RT,3),"min",sep=""),xlab=paste("QMass:",curQmz," R-squared:",round(summary(LinearRegression)$r.square,3),sep=""),ylab="")
			points(curUmass$concerntration,LinearRegression$fitted,type="l",col="red")
		}
		print(curUmass)

		

}

vizLinearQuantitationALL<-function(FileIDVsConcertration,alignedFeature,selectedCompoundRow,quantitationType)
{
	#plot linear quatitation
		#QMassVec<-strsplit(strMzs,",")[[1]]
		
		measure<-ifelse(quantitationType=="Height","int","pkArea")
		
		#curUmass<-subset(curAlignedFeature,clustID==curClustID&mz==curQmz)
		#select Umass
		nrows<-floor(sqrt(nrow(selectedCompoundRow)))
		ncols<-ceiling(nrow(selectedCompoundRow)/nrows)

		#x11()
		#png("~/Janus/cc.png",w=1024,h=768)
		par(mfrow=c(nrows,ncols))
		r.squared.tot<-0
		nComp<-nrow(selectedCompoundRow)
		for(i in 1:nComp)
		{
			curClustID<-selectedCompoundRow[i,]$clustID
			curCompName<-selectedCompoundRow[i,]$compoundName
			RT<-selectedCompoundRow[i,]$RT

			curQmz<-selectedCompoundRow[i,]$Qmass
			
			curAlignedFeature<-subset(alignedFeature,mz==curQmz&clustID==curClustID);
		
			#curUmass<-subset(curAlignedFeature,mz==curQmz)
			curUmass<-curAlignedFeature
			curUmass<-merge(curUmass,FileIDVsConcertration,by.y="fileid",by.x="FileID")
			curUmass<-subset(curUmass,select=c("concerntration",measure))
			LinearRegression<-lm(curUmass$int ~ curUmass$concerntration)
			r.squared<-summary(LinearRegression)$r.square
			r.squared.tot<-r.squared.tot+r.squared
			plot(curUmass,type="p",main=paste(curCompName," ",round(RT,3),"min",sep=""),xlab=paste("QMass:",curQmz," R-squared:",round(r.squared,3),sep=""),ylab="")
			points(curUmass$concerntration,LinearRegression$fitted,type="l",col="red")
		}
		print("average R.Squared:")
		print(r.squared.tot/nComp)

		

}
#vizAlignedEICs<-function(strMz,strClustID,alignedFeatureA,delaytime,RTstart,RTend)
#{	
#
#	plot(x=0,type="n",main="alignment results",xlim=c(minET,maxET),ylim=c(0,110),xlab="Retention Time(minutes)",ylab="Intensity",cex=0.03)
#
#	curFileID<-1
#	plotAlignedEIC(curFileID)
#	LocalFeatures=subset(alignedFeatureA,compoundET>(RTstart)&compoundET<(RTend))
#	
#	plotAlignedEIC<-function(curFileID)
#	{
#		fileName<-parseFileName(vizFileNames[curFileID])
#		dir1<-params$workingDir
#		EICFileName<-paste(dir1,"/output/EIC/",fileName,"EIC.cdf",sep="")
#		curNcid <- open.ncdf(EICFileName)			#read raw spectra
#		
#		vecInt <- get.var.ncdf(curNcid, varid="intVec")
#		vectorMz<-get.var.ncdf(curNcid, varid="mzVec")
#		totalscan<-length(vecInt)/length(vectorMz)
#		vecET<-((1:totalscan)-1)*ScanInterval+delaytime
#		filter by mz
#		mzVec<-unlist(strsplit(strMz,","))
#		nMz<-length(mzVec)
#		ifelse(nMz>0,curFeatureset<-subset(LocalFeatures,mz%in%mzVec&FileID==curFileID),curFeatureset<-subset(LocalFeatures,FileID==curFileID))
#		
#		curFeatureset<-subset(LocalFeatures,mz%in%mzVec)
#		write.csv(
#		filter by windowID
#		windowVec<-unlist(strsplit(strWindowID,","))
#		nWindow<-length(windowVec)
#		if(nWindow>0)
#		{
#			curFeatureset<-subset(curFeatureset,windowID%in%windowVec)
#			filter by compoundID
#			compVec<-unlist(strsplit(strComponentID,","))
#			nComp<-length(compVec)
#			if(nComp>0)
#				curFeatureset<-subset(curFeatureset,compoundID%in%compVec)
#		
#		}
#		
#		if no mz specified, then all mz are selected within the windows
#		if(nMz==0)
#		{
#			mzVec<-unique(curFeatureset$mz)
#			nMz<-length(mzVec)
#		}
#		curFeatureset$key<-paste(curFeatureset$windowID,curFeatureset$compoundID)
#	
#		curFeatureset
#		get the boundary of EIC
#		minLbound<-min(curFeatureset$lboundInd)
#		maxRbound<-max(curFeatureset$rboundInd)
#		minET<-(minLbound-1)*ScanInterval+delaytime
#		maxET<-	(maxRbound-1)*ScanInterval+delaytime
#		maxInt<-max(curFeatureset$Intensity)
#		
#		get RT vec
#		curFeatureset$RT<-(curFeatureset$pkInd-1)*ScanInterval+delaytime
#		normalize int
#		curFeatureset$Intensity<-100*curFeatureset$Intensity/maxInt
#		componentFeatures<-subset(curFeatureset,compoundID>0)
#		get data points
#		datapoints<-subset(componentFeatures,select=c("RT","Intensity"))
#		datalabels<-paste("window:",componentFeatures$windowID," component:",componentFeatures$compoundID," mz:",componentFeatures$mz," scan:",componentFeatures$pkInd,sep="")
#		
#	
#		xy<-par("usr")
#		xlength<-xy[2]-xy[1] 
#		ylength<-xy[4]-xy[3] 
#		plot EIC
#		for(i in 1:nMz)
#		{
#			mz<-mzVec[i]
#			get EIC
#			vecInt <- get.var.ncdf(ncid, varid=as.character(mz))
#			vecET<-((1:length(vecInt))-1)*ScanInterval+delaytime
#			mzInd<-which(vectorMz==mz)
#			startInd<-(mzInd-1)*totalscan+1
#			endInd<-startInd+totalscan-1
#			curVecInt <- vecInt[startInd:endInd]
#			
#			EIC<-data.frame(RT=vecET,int=curVecInt)
#			EIC<-subset(EIC,RT>=minET&RT<=maxET)
#			EICind<-which(vecET>=minET&vecET<=maxET)
#			get peak apex
#			curFeatures<-subset(curFeatureset,mz==as.numeric(mz))
#			plot mz legend
#			mtext(las = 1,at=i*ylength/30,side=4,text=mz,col=rainbow(nMz+1)[i+1])
#			plot EIC
#			points(vecET[EICind],100*curVecInt[EICind]/maxInt,type="l",col=rainbow(nMz+1)[i+1])
#
#			
#		}
#		plot deconvoluted peak apex
#		pchVec<-19:25
#		npchVec<-length(pchVec)
#
#		
#		ComponentID<-unique(componentFeatures$key)
#		nComponent<-length(ComponentID)
#		for(i in 1:nComponent)
#		{
#
#			curComponentID<-ComponentID[i]
#			fInd<-which(componentFeatures$key==curComponentID)
#			curFeatures<-subset(componentFeatures,key==curComponentID)			
#			points(componentFeatures[fInd,]$RT,componentFeatures[fInd,]$Intensity,pch=pchVec[i%%npchVec],cex=0.7)#col=rainbow(5)[i%%5]
#		}
#
#	}
#	
#	
#	
#	
#	playwith(
#	{
#
#	
#		plot denoised peak apex
#		
#		DenoisedFeatureset<-subset(curFeatureset,compoundID<0)
#		points(DenoisedFeatureset$RT,DenoisedFeatureset$Intensity,pch=4,col="red",cex=0.7)#col=rainbow(5)[i%%5]
#	}
#	,data.points=datapoints,labels=datalabels,new=TRUE)
#	
#}
plotAlignedEIC_singleCompound<-function(curFeature,delaytimeTable,type,curColor)
{


		fileInfo<-subset(delaytimeTable,FileID==curFeature$FileID)
		fileName<-fileInfo$FileName
		delaytime<-fileInfo$firstRT
		ScanInterval<-fileInfo$ScanInterval
		dir1<-params$WorkDir
		EICFileName<-paste(dir1,"/output/EIC/",fileName,"EIC.cdf",sep="")
		curNcid <- open.ncdf(EICFileName)			#read raw spectra
		
		vecInt <- get.var.ncdf(curNcid, varid="intVec")
		vectorMz<-get.var.ncdf(curNcid, varid="mzVec")
		totalscan<-length(vecInt)/length(vectorMz)
		vecET<-((1:totalscan)-1)*ScanInterval+delaytime
		
	
		#get the boundary of EIC
		maxInt<-max(curFeature$int)
		minET<-(curFeature$lboundInd-1)*ScanInterval+delaytime
		maxET<-	(curFeature$rboundInd-1)*ScanInterval+delaytime
		#get data points
#		datapoints<-subset(curFeatureset,select=c("RT","Intensity"))
#		datalabels<-paste("window:",curFeatureset$windowID," component:",curFeatureset$compoundID," mz:",curFeatureset$mz," scan:",curFeatureset$pkInd,sep="")
#		
	
		xy<-par("usr")
		xlength<-xy[2]-xy[1] 
		ylength<-xy[4]-xy[3] 
		#plot EIC
		
		
		mzInd<-which(vectorMz==curFeature$mz)
		startInd<-(mzInd-1)*totalscan+1
		endInd<-startInd+totalscan-1
		curVecInt <- vecInt[startInd:endInd]
		EICind<-which(vecET>=minET&vecET<=maxET)
		#mtext(las = 1,at=i*ylength/30,side=4,text=curFeature$FileID,col=rainbow(nFeatures)[i])
		if(type=="Before Alignment")
		{
			points(vecET[EICind],100*curVecInt[EICind]/maxInt,type="l",col=curColor)
		}else{			
			RTShift<-curFeature$alignedET-curFeature$ET
			points(vecET[EICind]+RTShift,100*curVecInt[EICind]/maxInt,type="l",col=curColor)
			}
}

plotAlignedEIC_multipleCompounds<-function(clustIDVec,LocalFeatures,type,mzVec)
{
	minLbound<-min(LocalFeatures$lboundInd)
	maxRbound<-max(LocalFeatures$rboundInd)
	minET<-(minLbound-1)*ScanInterval+max(delaytimeTable$firstRT)
	maxET<-	(maxRbound-1)*ScanInterval+min(delaytimeTable$firstRT)
	plot(x=0,type="n",main=type,xlim=c(minET,maxET),ylim=c(0,110),xlab="Retention Time(min)",ylab="Intensity",cex=0.03)
	for(j in 1:length(clustIDVec))
	{
		curClustID<-clustIDVec[j]
		curMz<-mzVec[j]
		curFeatures<-subset(LocalFeatures,clustID==curClustID&mz==curMz)
		curFeatures<-curFeatures[order(curFeatures$FileID),]

		for(i in 1:nrow(curFeatures))
		{
			curFeature<-curFeatures[i,]
			nFeatures<-nrow(curFeatures)
			curColor<-rainbow(nFeatures)[i]
			plotAlignedEIC_singleCompound(curFeature,delaytimeTable,type=type,curColor)
		}
	}
}			
	
vizAlignedEICs<-function(alignedFeatureA,strClustIDs,strMzs)
{	
	clustIDVec<-unlist(strsplit(strClustIDs,split=","))
	mzVec<-unlist(strsplit(strMzs,split=","))
	LocalFeatures<-subset(alignedFeatureA,clustID%in%clustIDVec)
	par(mfrow=c(2,1))
	plotAlignedEIC_multipleCompounds(clustIDVec,LocalFeatures,type="Before Alignment",mzVec)
	plotAlignedEIC_multipleCompounds(clustIDVec,LocalFeatures,type="After Alignment",mzVec)
}

vizRTdeviation<-function(alignedFeatureA,strClustIDs,strMzs)
{
	clustIDVec<-unlist(strsplit(strClustIDs,split=","))
	mzVec<-unlist(strsplit(strMzs,split=","))
	if(length(clustIDVec)>0)
	{
		LocalFeatures<-subset(alignedFeatureA,clustID%in%clustIDVec)
	}else{
		#add all compound part
	}
	minLbound<-min(LocalFeatures$lboundInd)
	maxRbound<-max(LocalFeatures$rboundInd)
	minET<-(minLbound-1)*ScanInterval+max(delaytimeTable$firstRT)
	maxET<-	(maxRbound-1)*ScanInterval+min(delaytimeTable$firstRT)

	plot(x=0,type="n",main=type,xlim=c(minET,maxET),ylim=c(-60,60),xlab="Retention Time(min)",ylab="RT Deviation(scan)",cex=0.03)
	for(j in 1:length(clustIDVec))
	{
		curClustID<-clustIDVec[j]
		curMz<-mzVec[j]
		curFeatures<-subset(LocalFeatures,clustID==curClustID&mz==curMz)
		curFeatures<-curFeatures[order(curFeatures$FileID),]

		for(i in 1:nrow(curFeatures))
		{
			curFeature<-curFeatures[i,]
			nFeatures<-nrow(curFeatures)
			curColor<-rainbow(nFeatures)[i]
			points(x=curFeature$ET,y=(curFeature$alignedET-curFeature$ET)*1200,pch=j,col=curColor)
		}
	}
}
	
plotAlignedSpec<-function(componentFeatures)
{


	
	componentFeatures<-componentFeatures[order(componentFeatures$FileID),]
	componentFeatures$key<-paste(componentFeatures$FileID,componentFeatures$CompoundID)
	tt<-unique(subset(componentFeatures,select=c("clustID","FileID","isRef","compoundET","key")))#get common mother Ions which occur in 80% of datasets
	minET<-min(tt$compoundET)
	maxInt<-max(componentFeatures$int)
	FileIds<-sort(unique(componentFeatures$FileID))
	#CompoundIDs<-sort(unique(componentFeatures$CompoundID))
	#nCompound<-length(CompoundIDs)
	nfile<-length(FileIds)
	maxMz<-max(componentFeatures$mz)
	minMz<-min(componentFeatures$mz)
	keys<-unique(tt$key)
	nkey<-length(keys)
	
	#par(mfrow=c(3,ceiling(nfile/3)))
	strOutput<-NULL
	#x11()
	par(mfrow=c(3,ceiling(nkey/3)))

	isRedmarked<-FALSE

	for(i in 1:nkey)
	{
		#curComponents<-subset(componentFeatures,FileID==FileIds[i])
		curComponents<-subset(componentFeatures,key==keys[i])
		curCompoundET<-unique(curComponents$compoundET)
		curClustID<-unique(curComponents$clustID)
		componentName<-paste(curCompoundET,curClustID)
		strOutput<-paste(strOutput,getMSPBlock(curComponents,curCompoundET*60,componentName,i),sep="\n\n")

		curFileID<-unique(curComponents$FileID)
		maxInt<-max(curComponents$int)
		color<-ifelse(unique(curComponents$isRef)==1,"blue","black")
	
		if(!isRedmarked&&unique(curComponents$compoundET)==minET)
		{
			color<-"red"
			isRedmarked<-TRUE
		}
		spec<-subset(curComponents,select=c("mz","int"))
		spec$int<-100*spec$int/maxInt
		
	
		plot(x=0,type="n",xlim=c(minMz,maxMz),ylim=c(1,100),main=paste("sample_",curFileID,"RT:",curCompoundET,sep=""),xlab="m/z",ylab="Intensity",cex=0.03,new=TRUE)
		points(spec,col=color,type="h")
		
	
	}

	#cat(paste(strOutput,"\n\n"))
	#cal dist matrix for spectra

	library(gtools)
	#FileIds<-sort(unique(componentFeatures$FileID))
	#nfile<-length(FileIds)
	r<-specDistMatrix(componentFeatures,addRTScore=FALSE)
	
	print("average distance:\n")
	for(col_ind in 1:ncol(r))
	{
		print(paste(unlist(strsplit(colnames(r)[col_ind],split=" "))[1],format(mean(r[-col_ind,col_ind]),digit=3),sep=":"))

		}
		

	print("distance matrix:\n")
	print(r)



	#print(avgScore)
	
	w2<-gwindow()
	gtext(text=strOutput,cont=w2)
    #distance <- as.dist(r)#calculate distance among each pair
	#distance


	#write spec to MSP format
	#NistoutputName<-"paper/JPR_April/9.72.msp"

#	for(i in 1:nfile)
#	{
#		curComponent<-subset(componentFeatures,FileID==FileIds[i])
#		strOutput<-paste(strOutput,getMSPBlockFromAlignedFeatures(curComponent,i),sep="\n\n")
#	}
#	
	#write(strOutput,file=NistoutputName)

}
