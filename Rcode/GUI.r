
LoadEIC<-function(h,...)
{
	svalue(sb) <<- paste("start the loading EIC...")
	begTime <- Sys.time()
	
	#read EIC 
	fileName<-parseFileName(vizFileNames)
	dir1<-params$WorkDir
	EICFileName<-paste(dir1,"/output/EIC/",fileName,"EIC.cdf",sep="")
	ncid <<- open.ncdf(EICFileName)			#read raw spectra


	runTime <- Sys.time()-begTime

	svalue(sb) <<- paste("Loaded EIC of:",fileName,runTime)
	
}


LoadPeaklist<-function(h,...)
{
	svalue(sb) <<- paste("start the loading peak list...")
	begTime <- Sys.time()
	
	fileName<-parseFileName(vizFileNames)
	dir1<-params$WorkDir
	#read peak list
	FeatureFileName<-paste(dir1,"/output/peakpicking/",fileName,"PeakList.csv",sep="")
	Features<<-read.csv(FeatureFileName)	
	#read deconvluted peak list
	FeatureFileName<-paste(dir1,"/output/decon/",fileName,"Decon.csv",sep="")
	deconFeatures<<-read.csv(FeatureFileName)
	

	runTime <- Sys.time()-begTime

	svalue(sb) <<- paste("Loaded peak list of:",fileName,runTime)
	
}
startPreprocess<-function(h,...)
{
	params$clustType<<-svalue(clust.type)
	params$nNode<<-svalue(clust.number)
	params$EicPkRatio<<-svalue(EIC.pkRatio)
	params$EICWindow<<-svalue(EIC.window)
	params$startET<<-svalue(EIC.start)
	params$endET<<-svalue(EIC.end)
	params$nFragments<<-svalue(Frag.number)
	params$maxShift<<-svalue(Frag.shift)
	params$FragPkAreaRatio<<-svalue(Frag.pkRatio)
	

	#enabled(h$obj$start)<-FALSE
	svalue(sb) <<- paste("start the preprocessing...")
	begTime <- Sys.time()
	#parGcFeatureDetect_Snow(params)
	runTime <- Sys.time()-begTime
	#enabled(h$obj$start)<-TRUE
	svalue(sb) <<- paste("Done:",params$ProcessedFiles,runTime)
	
}
startAlignment<-function(h,...)
{
	params$clustType<<-svalue(clust.type)
	params$nNode<<-svalue(clust.number)
	params$Max_Compound_ET_shift<<-svalue(gspMaxCompshift)
	params$maxSpecAngle<<-svalue(gspMaxSpecAngle)
	params$JobName<<-svalue(gtJobName)
	params$nFragments<<-svalue(Frag.number)
	svalue(sb) <<- paste("start the alignment...")
	begTime <- Sys.time()
	GcCompoundAlignment(params)
	#Sys.sleep(1)
	runTime <- Sys.time()-begTime
	svalue(sb) <<- paste("Done:",params$ProcessedFiles,runTime)
	
}
searchUMass<-function(h,...)
{
	params$clustType<<-svalue(clust.type)
	params$nNode<<-svalue(clust.number)
	params$occuRatio<<-svalue(gspRatio)
	params$uMassWindow<<-svalue(gspUMassWindow)
	params$JobName<<-svalue(gtJobName)
	params$nFragments<<-svalue(Frag.number)
	svalue(sb) <<- paste("searching Unique mass...")
	begTime <- Sys.time()
	GcUniqueIonTable(params)
	#Sys.sleep(1)
	runTime <- Sys.time()-begTime
	svalue(sb) <<- paste("Done:",params$ProcessedFiles,runTime)
	
}

plotEICPeak<-function(h,...)
{
	strMz<-svalue(gtMz)
	#delaytime<-as.numeric(svalue(gtDelayTime))
	
	delaytime<-subset(delaytimeTable,FileName==parseFileName(vizFileNames))$firstRT
	ScanInterval<-subset(delaytimeTable,FileName==parseFileName(vizFileNames))$ScanInterval
	LocalFeatures=subset(Features,(pkInd-1)>((svalue(RTstart)-delaytime)/ScanInterval)&(pkInd-1)<((svalue(RTend)-delaytime)/ScanInterval))
	vizEICPeak(strMz,curFeatureset=LocalFeatures,ncid,delaytime,ScanInterval)

	
}
plotDeconEIC<-function(h,...)
{
	strMz<-svalue(gtMz)
	strWindowID<-svalue(gtWindowID)	
	strComponentID<-svalue(gtComponentID)
	#delaytime<-as.numeric(svalue(gtDelayTime))
	delaytime<-subset(delaytimeTable,FileName==parseFileName(vizFileNames))$firstRT
	ScanInterval<-subset(delaytimeTable,FileName==parseFileName(vizFileNames))$ScanInterval
	LocaldeconFeatures=subset(deconFeatures,(pkInd-1)>((svalue(RTstart)-delaytime)/ScanInterval)&(pkInd-1)<((svalue(RTend)-delaytime)/ScanInterval))
	vizDecon(strMz,strWindowID,strComponentID,LocaldeconFeatures,ncid,delaytime,ScanInterval)
}

#plotAlignedEIC<-function(h,...)
#{
#	strMz<-svalue(gtMz)
#	strWindowID<-svalue(gtWindowID)
#	delaytime<-subset(delaytimeTable,FileName==parseFileName(vizFileNames))$firstRT
#	LocalAlignedFeatures=subset(alignedFeature,pkInd>((svalue(RTstart)-delaytime)/ScanInterval)&pkInd<((svalue(RTend)-delaytime)/ScanInterval))
#	vizDecon(strMz,strWindowID,LocaldeconFeatures,ncid,delaytime)
#}

plotSpec<-function(h,...)
{
	#strMz<-svalue(gtMz)
	strWindowID<-svalue(gtWindowID)
	strComponentID<-svalue(gtComponentID)
	delaytime<-subset(delaytimeTable,FileName==parseFileName(vizFileNames))$firstRT
	ScanInterval<-subset(delaytimeTable,FileName==parseFileName(vizFileNames))$ScanInterval
	vizSpec(strWindowID,strComponentID,deconFeatures,delaytime,ScanInterval)
}
plotPCA<-function(h,...)
{

	PCA(trainingSet,groupVec,colVec)
}
plotPLSDA<-function(h,...)
{

	PLSDA(trainingSet,groupVec,colVec)
}
plotClusterA<-function(h,...)
{

	ClusterA(trainingSet,groupVec,colVec)
}

plotQuantitation<-function(h,...)
{
	
	vizLinearQuantitation(FileIDVsConcertration,alignedFeature,selectedCompoundRow,strMzs=svalue(gtQMz),svalue(quantitation.type));

}
hPlotQuantitationAll<-function(h,...)
{
	
	vizLinearQuantitationALL(FileIDVsConcertration,alignedFeature,selectedCompoundRow,quantitationType=svalue(quantitation.type));

}
selectEICPeakData<-function(h,...)
{
	visible(w,FALSE)
	w1 = gbasicdialog(title="select raw data files to process ",handler=function(h,...)params$DataFiles<<-tag(h$obj,1))
	w1 = gwindow(title="select sample data files to plot ")
	size(w1) <- c(500,500)
	g1 <- ggroup(horizontal=FALSE,cont = w1)
	gtData<-gtable(items=list.files(full.names=TRUE,paste(params$WorkDir,"/raw",sep="")),multiple=TRUE,cont=g1)
	size(gtData)<-c(500,400)
	g2 <- ggroup(cont = g1)
	gbOk<-gbutton(action=gtData,text="ok",border=TRUE,cont=g2,handler=function(h,...){vizFileNames<<-svalue(h$action);dispose(w1);svalue(sb) <<- paste("selected",params$DataFiles);visible(w,TRUE)})
	gbCancel<-gbutton(text="cancel",handler=function(h,...){dispose(w1);visible(w,set=TRUE)},border=TRUE,cont=g2)	
	
	svalue(sb) <<- paste("You selected",vizFileNames)
}


LoadConcerntrationData<-function(h,...)
{
	visible(w,FALSE)
	curRelativePath<-gsub(getwd(),"",params$WorkDir)
	concerntrationFile<-gfile(initialfilename=paste(curRelativePath,"output",sep="/"),text="select Concerntration data file")
	FileIDVsConcertration<-read.csv(concerntrationFile)
	FileIDVsConcertration<<-FileIDVsConcertration[order(FileIDVsConcertration$concerntration),]
	svalue(sb) <<- paste("loaded:",concerntrationFile)

	visible(w,TRUE)
}
LoadmatchResult<-function(h,...)
{
	visible(w,FALSE)
	curRelativePath<-gsub(getwd(),"",params$WorkDir)
	matchResult<-gfile(initialfilename=paste(curRelativePath,"output",sep="/"),text="select matchResult file")
	
	clustVsCompName<-read.csv(matchResult)
	#clustIDs<-NULL
#	for(i in 1:nrow(clustVsCompName))
#	{
#
#		clustID<-unlist(strsplit(as.character(clustVsCompName[i,]$componentName)," "))[2]
#		clustIDs<-c(clustIDs,substr(clustID,6,nchar(clustID)))
#		}
#
#
#	clustVsCompName$clustID<-clustIDs
	clustVsCompName<<-clustVsCompName
	svalue(sb) <<- paste("loaded:",matchResult)

	visible(w,TRUE)
}
LoadAlignedFeatureQ<-function(h,...)
{
	visible(w,FALSE)
	dir1<-params$WorkDir
	
	curRelativePath<-gsub(getwd(),"",dir1)
	alignedFeatureFile<-gfile(initialfilename=paste(curRelativePath,"output",sep="/"),text="select aligned data file")
	alignedFeature<<-read.csv(alignedFeatureFile)
	alignedFeature<<-subset(alignedFeature,clustID%in%clustVsCompName$clustID)
	
	svalue(sb) <<- paste("loaded:",alignedFeatureFile)

	visible(w,TRUE)
}
LoadAllEIC<-function(h,...)
{
	#visible(w,FALSE)
	dir1<-params$WorkDir
	ncidVec<-vector(mode="list",length=nrow(delaytimeTable))
	for(i in 1:nrow(delaytimeTable))
	{
		filename<-delaytimeTable[i,]$FileName
		fileID<-delaytimeTable[i,]$FileID
		EICFile<-paste(dir1,"/output/EIC/",filename,"EIC.cdf",sep="")
		ncidVec[[i]] <- open.ncdf(EICFile)			#read raw spectra
		print(paste(filename,"loaded.\n"))
	}
	ncidVec<<-ncidVec
	svalue(sb) <<- paste("All EIC loaded:")

	visible(w,TRUE)
}
LoadAlignedFeature<-function(h,...)
{
	visible(w,FALSE)
	dir1<-params$WorkDir
	
	curRelativePath<-gsub(getwd(),"",dir1)
	alignedFeatureFile<-gfile(initialfilename=paste(curRelativePath,"output",sep="/"),text="select aligned data file")
	alignedFeatureA<<-read.csv(alignedFeatureFile)
	
	svalue(sb) <<- paste("loaded:",alignedFeatureFile)

	visible(w,TRUE)
}
LoadbioMarker<-function(h,...)
{
	visible(w,FALSE)
	curRelativePath<-gsub(getwd(),"",params$WorkDir)
	biomarkerFile<-gfile(initialfilename=paste(curRelativePath,"output",sep="/"),text="select biomarker file")
	biomarker<<-read.csv(biomarkerFile)
	biomarker<<-subset(biomarker,clustID%in%clustVsCompName$clustID)
	
	svalue(sb) <<- paste("loaded:",biomarkerFile)

	visible(w,TRUE)
}


	
getRawData<-function(h,...)
{
	visible(w,FALSE)
	w1 = gwindow(title="select raw data files to process ")
	size(w1) <- c(500,500)
	g1 <- ggroup(horizontal=FALSE,cont = w1)
	gtData<-gtable(items=list.files(full.names=TRUE,paste(params$WorkDir,"/raw",sep="")),multiple=TRUE,cont=g1)
	size(gtData)<-c(500,400)
	g2 <- ggroup(cont = g1)
	gbOk<-gbutton(action=gtData,text="ok",border=TRUE,cont=g2,handler=function(h,...){params$DataFiles<<-svalue(h$action);dispose(w1);svalue(sb) <<- paste("selected",params$DataFiles);visible(w,TRUE)})
	gbCancel<-gbutton(text="cancel",handler=function(h,...){dispose(w1);visible(w,set=TRUE)},border=TRUE,cont=g2)
	
	svalue(sb) <<- paste("You selected",params$DataFiles)
}
getProcessedData<-function(h,...)
{
	visible(w,FALSE)
	
	w1 = gwindow(title="select processed data files to align ")
	size(w1) <- c(500,500)
	dir1<-params$WorkDir
	g1 <- ggroup(horizontal=FALSE,cont = w1)
	gtData<-gtable(items=list.files(pattern = ".csv",full.names=TRUE,paste(dir1,"/output/decon",sep="")),multiple=TRUE,cont=g1)
	size(gtData)<-c(500,400)
	g2 <- ggroup(cont = g1)
	gbOk<-gbutton(action=gtData,text="ok",border=TRUE,cont=g2,handler=function(h,...){params$ProcessedFiles<<-svalue(h$action);dispose(w1);svalue(sb) <<- paste("selected",params$ProcessedFiles);visible(w,TRUE)})
	gbCancel<-gbutton(text="cancel",handler=function(h,...){dispose(w1);visible(w,set=TRUE)},border=TRUE,cont=g2)
	
	#svalue(sb) <<- paste("You selected",params$DataFiles)
}
SelectCompoundtoQuantify<-function(h,...)
{
	visible(w,FALSE)
	
	w1 = gwindow(title="select compound: ")
	size(w1) <- c(500,500)
	dir1<-params$WorkDir
	g1 <- ggroup(horizontal=FALSE,cont = w1)
	gtData<-gtable(items=clustVsCompName,chosencol=c(1:ncol(clustVsCompName)),multiple=TRUE,drop=FALSE,index=FALSE,cont=g1)
	size(gtData)<-c(500,400)
	g2 <- ggroup(cont = g1)
	gbOk<-gbutton(action=gtData,text="ok",border=TRUE,cont=g2,
				handler=function(h,...)
					{
						selectedCompoundRow<<-svalue(h$action);
						#clustVsCompName[CompoundInd,]$clustID
						dispose(w1);
						visible(w,TRUE)
					}
				)
	gbCancel<-gbutton(text="cancel",handler=function(h,...){dispose(w1);visible(w,set=TRUE)},border=TRUE,cont=g2)
	
	#svalue(sb) <<- paste("You selected",params$DataFiles)
}

SelectQmass<-function(h,...)
{
	visible(w,FALSE)
	
	w1 = gwindow(title="select QMass: ")
	size(w1) <- c(500,500)

	g1 <- ggroup(horizontal=FALSE,cont = w1)
	

	curbiomarker<-subset(biomarker,clustID==selectedCompoundRow$clustID&isUnique==1);

	
	gtData<-gtable(items=sort(unique(curbiomarker$mz)),multiple=TRUE,cont=g1)
	size(gtData)<-c(500,400)
	g2 <- ggroup(cont = g1)
	gbOk<-gbutton(action=gtData,text="ok",border=TRUE,cont=g2,
					handler=function(h,...)
					{

						QMassVec<-svalue(h$action)
						strMzs<-QMassVec[1]
						for(i in QMassVec[-1])
						{
							strMzs<-paste(strMzs,i,sep=",")
						}
						svalue(gtQMz)<<-strMzs
						#vizLinearQuantitation(FileIDVsConcertration,alignedFeature,selectedCompoundRow,QMassVec,svalue(quantitation.type));
						dispose(w1);
						visible(w,TRUE);
					}
				)
	gbCancel<-gbutton(text="cancel",handler=function(h,...){dispose(w1);visible(w,set=TRUE)},border=TRUE,cont=g2)
	
	#svalue(sb) <<- paste("You selected",params$DataFiles)
}

SelectMz<-function(h,...)
{
	visible(w,FALSE)
	
	w1 = gwindow(title="select m/z: ")
	size(w1) <- c(500,500)

	g1 <- ggroup(horizontal=FALSE,cont = w1)
	



	
	gtData<-gtable(items=sort(unique(Features$curMz)),multiple=TRUE,cont=g1)
	size(gtData)<-c(500,400)
	g2 <- ggroup(cont = g1)
	gbOk<-gbutton(action=gtData,text="ok",border=TRUE,cont=g2,
					handler=function(h,...)
					{

						QMassVec<-svalue(h$action)
						strMzs<-QMassVec[1]
						for(i in QMassVec[-1])
						{
							strMzs<-paste(strMzs,i,sep=",")
						}
						svalue(gtMz)<<-strMzs
						#vizLinearQuantitation(FileIDVsConcertration,alignedFeature,selectedCompoundRow,QMassVec,svalue(quantitation.type));
						dispose(w1);
						visible(w,TRUE);
					}
				)
	gbCancel<-gbutton(text="cancel",handler=function(h,...){dispose(w1);visible(w,set=TRUE)},border=TRUE,cont=g2)
	
	#svalue(sb) <<- paste("You selected",params$DataFiles)
}


setWorkDir<-function(h,...)
{
	visible(w,FALSE)
	params$WorkDir<<-gfile(initialfilename="obesity",text="select a working directory",type="selectdir");
	params$JobName<<-svalue(gtJobName)
	svalue(sb) <<- paste("set working directory as:",params$WorkDir)
	svalue(gedtDir)<<-params$WorkDir
	delaytimeTable<<-read.csv(paste(params$WorkDir,"/output/",params$JobName,"_RT.csv",sep=""))
	visible(w,TRUE)
}

getAlignedFile<-function(h,...)
{
	visible(w,FALSE)
	dir1<-params$WorkDir
	curRelativePath<-gsub(getwd(),"",dir1)
	
	params$AlignedFile<<-gfile(initialfilename=paste(curRelativePath,"output",sep="/"),text="select a aligned data file");
	svalue(sb) <<- paste("selected:",params$AlignedFile)
	svalue(gedtAignedFile)<<-params$AlignedFile
	visible(w,TRUE)
}
#LoadAlignedFileViz<-function(h,...)
#{
#	visible(w,FALSE)
#	dir1<-params$WorkDir
#	curRelativePath<-gsub(getwd(),"",dir1)
#	
#	AlignedFile<-gfile(initialfilename=paste(curRelativePath,"output",sep="/"),text="select a aligned data file");
#	alignedFeature<-read.csv(AlignedFile)
#	svalue(sb) <<- paste("Loaded:",AlignedFile)
#	visible(w,TRUE)
#}
getRTFilterFile<-function(h,...)
{
	visible(w,FALSE)
	dir1<-params$WorkDir
	curRelativePath<-gsub(getwd(),"",dir1)
	RTFilterFile<-gfile(initialfilename=paste(curRelativePath,"output",sep="/"),text="select a RT Filter file");
	svalue(sb) <<- paste("selected:",params$RTFilterFile)
	svalue(gedtRTFilterFile)<<-RTFilterFile
	visible(w,TRUE)
}
getlibFile<-function(h,...)
{
	visible(w,FALSE)
	dir1<-params$WorkDir
	curRelativePath<-gsub(getwd(),"",dir1)
	libfile<-gfile(initialfilename=paste(curRelativePath,"output",sep="/"),text="select a lib file");
	svalue(sb) <<- paste("selected:",libfile)
	svalue(gedtlibFile)<<-libfile
	visible(w,TRUE)
}
getNISTFile<-function(h,...)
{
	visible(w,FALSE)
	dir1<-params$WorkDir
	curRelativePath<-gsub(getwd(),"",dir1)
	NISTFile<-gfile(initialfilename=paste(curRelativePath,"output",sep="/"),text="select a NIST output file");
	svalue(sb) <<- paste("selected:",NISTFile)
	svalue(gedtNISTFile)<<-NISTFile
	visible(w,TRUE)
}
getIDnistFile<-function(h,...)
{
	visible(w,FALSE)
	dir1<-params$WorkDir
	curRelativePath<-gsub(getwd(),"",dir1)
	IDnistFile<-gfile(initialfilename=paste(curRelativePath,"output",sep="/"),text="select a NIST output file");
	svalue(sb) <<- paste("selected:",IDnistFile)
	svalue(gedtIDFile)<<-IDnistFile
	visible(w,TRUE)
}
showProprocessGroup<-function(h,...)
{
	#visible(ggUMass,FALSE)
	visible(ggAlignment,FALSE)
	visible(ggViz,FALSE)
	visible(ggQuantitation,FALSE)
	#visible(ggAlignmentVIz,TRUE)	
	visible(ggPreprocess,TRUE)
	#visible(ggIdentification,FALSE)			
	visible(ggMVA,FALSE)
	size(w)<<-c(989,373)
}

showAlignmentGroup<-function(h,...)
{
	#visible(ggUMass,FALSE)
	visible(ggPreprocess,FALSE)
	visible(ggViz,FALSE)
	visible(ggQuantitation,FALSE)
	#visible(ggAlignmentVIz,TRUE)	
	visible(ggAlignment,TRUE)
	#visible(ggIdentification,FALSE)			
	visible(ggMVA,FALSE)
	size(w)<<-c(1200,600)
}
#showUMassGroup<-function(h,...)
#{
#	visible(ggPreprocess,FALSE)
#	visible(ggAlignment,FALSE)
#	visible(ggViz,FALSE)
#	visible(ggQuantitation,FALSE)
	#visible(ggAlignmentVIz,TRUE)
	#visible(ggUMass,TRUE)
	#visible(ggIdentification,FALSE)			
#	visible(ggMVA,FALSE)
#	size(w)<<-c(743,500)
#}
showMVAGroup<-function(h,...)
{
	visible(ggPreprocess,FALSE)
	visible(ggAlignment,FALSE)
	visible(ggViz,FALSE)
	visible(ggQuantitation,FALSE)
	#visible(ggAlignmentVIz,FALSE)
	#visible(ggUMass,FALSE)
	#visible(ggIdentification,FALSE)			
	visible(ggMVA,TRUE)
	size(w)<<-c(743,200)
}
showVizGroup<-function(h,...)
{
	visible(ggPreprocess,FALSE)
	visible(ggAlignment,FALSE)
	#visible(ggUMass,FALSE)
	visible(ggViz,TRUE)
	visible(ggQuantitation,TRUE)
	#visible(ggAlignmentVIz,TRUE)
	#visible(ggIdentification,FALSE)		
	visible(ggMVA,FALSE)	
	size(w)<<-c(1400,800)
}
#showIDGroup<-function(h,...)
#{
#	visible(ggPreprocess,FALSE)
#	visible(ggAlignment,FALSE)
	#visible(ggUMass,FALSE)
#	visible(ggViz,FALSE)
#	visible(ggQuantitation,FALSE)
	#visible(ggAlignmentVIz,TRUE)
#	visible(ggMVA,FALSE)
	#visible(ggIdentification,TRUE)		
#	
#	size(w)<<-c(1400,800)
#}




#startMainWindow<-function()	
#{
	setwd("/Users/yni1/workspace/")
	source("/Users/yni1/workspace/ADAP/pipeline.r")
	library(gWidgetsRGtk2)
	library("snow")
	library("playwith")
	library("ncdf")
	options(guiToolkit = "RGtk2")
	#global variable
	
	params<-NULL

	#ncid<-NULL
	#Features<-NULL
	#deconFeatures<-NULL
	#vizFileNames<-NULL

	w = gwindow("GC Data Analysis Pipeline")
	#menu
	mbl <- list()
	mbl$task$Preprocess$handler=showProprocessGroup
	mbl$task$Alignment$handler=showAlignmentGroup
	#mbl$task$UMassSearch$handler=showUMassGroup
	#mbl$task$Identification$handler=showIDGroup  
	mbl$task$MVA$handler=showMVAGroup 
	mbl$task$visualization$handler=showVizGroup 
	mbl$quit$handler = function(h,...)dispose(w)
	mbl$quit$icon="quit"
	mb <- gmenu(mbl, container=w)

	ggGeneral<-ggroup(cont=w)
		gf1 <- gframe(text="Cluster",con=ggGeneral,horizontal=FALSE)
			glabel(text="Type:",con=gf1,anchor=c(-1,1))
			clust.type<-gradio(horizontal=TRUE,c("SOCK","MPI"),selected=2,con=gf1,,anchor=c(-1,1))
			glabel(text="Node:",con=gf1,anchor=c(-1,1))
			clust.number<-gspinbutton(from=1,to=32,by=1,value=8,con=gf1)
		gf4 <- gframe(con=ggGeneral,horizontal=FALSE)
		glabel(text="job name:",con=gf4,anchor=c(-1,1))
		gtJobName<-gedit(text=ifelse(is.null(params),"",params$JobName),con=gf4)
#		glabel(text="RT of first scan:",con=gf4,anchor=c(-1,1))
#		gtDelayTime<-gedit(text="4.083",con=gf4)
		glabel(text="Minimum number of fragments for each Compound",con=gf4,anchor=c(-1,1))
		Frag.number<-gspinbutton(from=3,to=10,by=1,value=5,con=gf4)
		gbDir<-gbutton(text="directory...",cont=gf4,handler=setWorkDir)
		gedtDir<-gedit(text=ifelse(is.null(params),"",params$WorkDir),con=gf4)
		enabled(gedtDir)<-FALSE
		#Sys.sleep(20)
	ggPreprocess<-ggroup(horizontal=FALSE,cont=w)
		gg1<-ggroup(cont=ggPreprocess)
		gf2 <- gframe(text="EIC Peak Picking",con=gg1,horizontal=FALSE)
		glabel(text="EIC local denoising window size(s):",con=gf2,anchor=c(-1,1))	
		EIC.window<-gspinbutton(from=1,to=120,by=5,value=30,con=gf2)
		glabel(text="Minimum EIC peak ratio:",con=gf2,anchor=c(-1,1))	
		EIC.pkRatio<-gspinbutton(from=0,to=1000,by=50,value=5,con=gf2)
		glabel(text="Retention Time",con=gf2,anchor=c(-1,1))
		glabel(text="From",con=gf2,anchor=c(-1,1))
		EIC.start<-gspinbutton(from=0,to=60,by=1,value=0,con=gf2)
		glabel(text="To",con=gf2,anchor=c(-1,1))
		EIC.end<-gspinbutton(from=0,to=60,by=1,value=60,con=gf2)
		
		gf3 <- gframe(text="deconvolution",con=gg1,horizontal=FALSE)
		glabel(text="Minimum fragments peak ratio:",con=gf3,anchor=c(-1,1))
		Frag.pkRatio<-gspinbutton(from=0,to=100,by=10,value=10,con=gf3)
		glabel(text="Maximum Fragments shift(scans)",con=gf3,anchor=c(-1,1))
		Frag.shift<-gspinbutton(from=0,to=5,by=1,value=3,con=gf3)
		
		gg2<-ggroup(cont=ggPreprocess)
		gbData<-gbutton(text="open...",cont=gg2,handler=getRawData)
		gbStart<-gbutton(text="start",cont=gg2,handler=startPreprocess)
	#enabled(ggPreprocess)<-FALSE
	#visible(ggPreprocess,FALSE)
	#Sys.sleep(2)
	ggAlignment<-ggroup(horizontal=TRUE,cont=w)
		gfAlignment<-gframe(text="",con=ggAlignment,horizontal=FALSE)
			gfAlignmentProcess<-gframe(text="Alignment",con=gfAlignment,horizontal=FALSE)
				glabel(text="Max Compound ET shift(s):",con=gfAlignmentProcess,anchor=c(-1,1))
				gspMaxCompshift<-gspinbutton(from=1,to=60,by=5,value=10,con=gfAlignmentProcess)
				glabel(text="spectrum similarity angle:",con=gfAlignmentProcess,anchor=c(-1,1))
				gspMaxSpecAngle<-gspinbutton(from=1,to=60,by=2,value=30,con=gfAlignmentProcess)
				gbAData<-gbutton(text="open...",cont=gfAlignmentProcess,handler=getProcessedData)
				gbAStart<-gbutton(text="start",cont=gfAlignmentProcess,handler=startAlignment)
				enabled(gfAlignmentProcess)<-FALSE
				
#			gfAlignmentViz <- gframe(text="Alignment results",con=gfAlignment,horizontal=FALSE)
			gbutton(text="Load aligned Data...",cont=gfAlignment,handler=LoadAlignedFeature)

			gbldEICs<-gbutton(text="Load EIC data...",cont=gfAlignment,handler=LoadAllEIC)
			enabled(gbldEICs)<-FALSE
			glabel(text="clustID:",con=gfAlignment,anchor=c(-1,1))
			gtClustID<-gedit(text="",con=gfAlignment)
			gfAlignedSpecViz <- gframe(text="view aligned spectra",con=gfAlignment,horizontal=FALSE)
				gbutton(text="plot",cont=gfAlignedSpecViz,handler=function(h,...){plotAlignedSpec(componentFeatures=subset(alignedFeatureA,clustID%in%(unlist(strsplit(svalue(gtClustID),",")))))})
			gfAlignedEIC <- gframe(text="view aligned EICs",con=gfAlignment,horizontal=FALSE)
				glabel(text="mz:",con=gfAlignedEIC,anchor=c(-1,1))
				gtAmz<-gedit(text="",con=gfAlignedEIC)
				gbutton(text="plot",cont=gfAlignedEIC,handler=function(h,...)
				{vizAlignedEICs(alignedFeatureA,strClustIDs=svalue(gtClustID),strMzs=svalue(gtAmz))})
			gfRTDeviation <- gframe(text="view RT Deviation",con=gfAlignment,horizontal=FALSE)
				gbutton(text="plot ",cont=gfRTDeviation,handler=function(h,...)
				{vizRTdeviation(alignedFeatureA,strClustIDs=svalue(gtClustID),strMzs=svalue(gtAmz))})

		gfpostAlignment<- gframe(text="post Alignment",con=ggAlignment,horizontal=FALSE)
			gbAlignedFile<-gbutton(text="select Aligned File...",cont=gfpostAlignment,handler=getAlignedFile)
			gedtAignedFile<-gedit(text=,con=gfpostAlignment)
			enabled(gedtAignedFile)<-FALSE
			glabel(text="minimum sample ratio:",con=gfpostAlignment,anchor=c(-1,1))
			gspRatio<-gspinbutton(from=0,to=1.0,by=0.1,value=0.5,con=gfpostAlignment)

			gfUmass<-gframe(text="Search UMASS",con=gfpostAlignment,horizontal=FALSE)
				glabel(text="Unique Mass search window(s):",con=gfUmass,anchor=c(-1,1))
				gspUMassWindow<-gspinbutton(from=1,to=20,by=1,value=6,con=gfUmass)
				gbUStart<-gbutton(text="start",cont=gfUmass,handler=searchUMass)

			gfMSP<-gframe(text="Output",con=gfpostAlignment,horizontal=FALSE)
				glabel(text="RT tolerance(s)",con=gfMSP,anchor=c(-1,1))
				gspRTFilterWindow<-gspinbutton(from=1,to=20,by=1,value=6,con=gfMSP)
				gbFilterbyRT<-gbutton(text="Open RT Filter Data...",cont=gfMSP,handler=getRTFilterFile)
				gedtRTFilterFile<-gedit(text="",con=gfMSP)
				enabled(gedtRTFilterFile)<-FALSE
				#gbNistFile<-gbutton(text="Open NIST MSP File...",cont=gfMSP,handler=getNISTFile)
				#gedtNISTFile<-gedit(text="",con=gfMSP)
				#enabled(gedtNISTFile)<-FALSE
				#gbutton(text="Filter MSP By RT",cont=gfMSP,handler=function(h,...)
				#{FilterByRT(mspFileName=svalue(gedtNISTFile),filterFile=svalue(gedtRTFilterFile),RT_Tolerance=svalue(gspRTFilterWindow))})
				gtFilterClustID<-gedit(text="",con=gfMSP)
				#gbutton(text="Filter MSP By ClustID",cont=gfMSP,handler=function(h,...)
				#{FilterByClustID(mspFileName=svalue(gedtNISTFile),strClustID=svalue(gtFilterClustID))})
				glabel(text="Filter Type:",con=gfMSP,anchor=c(-1,1))
				Filter.type<-gradio(horizontal=TRUE,c("RT","clustID","None"),con=gfMSP,,anchor=c(-1,1))

				gbutton(text="write reference Spectra to MSP",cont=gfMSP,
				handler=function(h,...){refSpec2MSP(sampleRatio=svalue(gspRatio),JobName=svalue(gtJobName),params,
				filterFile=svalue(gedtRTFilterFile),RT_Tolerance=svalue(gspRTFilterWindow),strClustID=svalue(gtFilterClustID),filterType=svalue(Filter.type))})		
	
			gfIdentification<- gframe(text="Identification",con=ggAlignment,horizontal=FALSE)
				#ggIdentification<-ggroup(horizontal=FALSE,cont=w)
				glabel(text="RT tolerance(s)",con=gfIdentification,anchor=c(-1,1))
				gspIDRTWindow<-gspinbutton(from=1,to=20,by=1,value=1,con=gfIdentification)
				
	
				
				
				#glabel(text="similarity tolerance(s)",con=ggIdentification,anchor=c(-1,1))
				#gedtIDAngle<-gedit(text="22",con=ggIdentification)
				gbutton(text="select lib file...",cont=gfIdentification,handler=getlibFile)
				gedtlibFile<-gedit(text="",con=gfIdentification)
				enabled(gedtlibFile)<-FALSE
				gbutton(text="select File to identify...",cont=gfIdentification,handler=getIDnistFile)
				gedtIDFile<-gedit(text="",con=gfIdentification)
				enabled(gedtIDFile)<-FALSE
				
				gbutton(text="start to evaluate",cont=gfIdentification,
				handler=function(h,...){evaluateIdentification(refSpec=svalue(gedtlibFile),
				inputSpec=svalue(gedtIDFile),RT_Tolerance=svalue(gspIDRTWindow),minSpecSimilarity=svalue(gspSpecSimilarity),fileName=parseFileName(svalue(gedtIDFile)),withRT=T)})
				
				glabel(text="similarity tolerance",con=gfIdentification,anchor=c(-1,1))
				gspSpecSimilarity<-gspinbutton(from=700,to=999,by=50,value=800,con=gfIdentification)
				gbutton(text="start to match",cont=gfIdentification,
				handler=function(h,...){result<-libMatching(refSpec=svalue(gedtlibFile),
				inputSpec=svalue(gedtIDFile),minSpecSimilarity=svalue(gspSpecSimilarity));
				write.csv(result,paste(params$WorkDir,"/output/libMatchResult_",parseFileName(svalue(gedtIDFile)),".csv",sep=""))})
			

		#visible(ggIdentification,FALSE)
	visible(ggAlignment,FALSE)
	
	
	
	#Sys.sleep(2)
	#ggUMass<-ggroup(horizontal=FALSE,cont=w)
#		glabel(text="Unique Mass search window(s):",con=ggUMass,anchor=c(-1,1))
#		gspUMassWindow<-gspinbutton(from=1,to=20,by=1,value=6,con=ggUMass)
#		glabel(text="fragment occurency ratio:",con=ggUMass,anchor=c(-1,1))
#		gspRatio<-gspinbutton(from=0,to=1.0,by=0.1,value=0.5,con=ggUMass)
#		gbAlignedFile<-gbutton(text="Open Aligned Data...",cont=ggUMass,handler=getAlignedFile)
#		gedtAignedFile<-gedit(text=,con=ggUMass)
#		enabled(gedtAignedFile)<-FALSE
#		gbUStart<-gbutton(text="start",cont=ggUMass,handler=searchUMass)
#	visible(ggUMass,FALSE)
	
	ggViz<-ggroup(horizontal=FALSE,cont=w)
		#ggEIC<-ggroup(cont=ggViz)
	
		
		gbData<-gbutton(text="select Sample...",cont=ggViz,handler=selectEICPeakData)
		gfTIC<- gframe(text="TIC",con=ggViz,horizontal=FALSE)
			gbEICPlot<-gbutton(text="plot",cont=gfTIC,handler=function(h,...){vizTIC(inFilePath=vizFileNames)})

		gfEIC <- gframe(text="visualize Preprocessing results",con=ggViz,horizontal=FALSE)		
			gbLoadEIC<-gbutton(text="load EIC",cont=gfEIC,handler=LoadEIC)
			gbLoadPeak<-gbutton(text="load Peak Info",cont=gfEIC,handler=LoadPeaklist)
			
			gbutton(text="select M/z...",cont=gfEIC,handler=SelectMz)
			#glabel(text="M/z:(seperate by comma)",con=gfEIC,anchor=c(-1,1))	
			gtMz<-gedit(text="73",con=gfEIC)
			glabel(text="window ID:(seperate by comma)",con=gfEIC,anchor=c(-1,1))	
			gtWindowID<-gedit(text="",con=gfEIC)	
			glabel(text="Component ID:",con=gfEIC,anchor=c(-1,1))	
			gtComponentID<-gedit(text="",con=gfEIC)
					
			
			gfPeakPicking<- gframe(text="peak Picking results",con=gfEIC,horizontal=FALSE)
				gbEICPlot<-gbutton(text="plot",cont=gfPeakPicking,handler=plotEICPeak)
			
			gfDecon<- gframe(text="deconvolution results",con=gfEIC,horizontal=FALSE)
				glabel(text="Retention Time (min)",con=gfDecon,anchor=c(-1,1))
				glabel(text="From",con=gfDecon,anchor=c(-1,1))
				RTstart<-gspinbutton(from=0,to=60,by=0.01,value=0,con=gfDecon)
				glabel(text="To",con=gf2,anchor=c(-1,1))
				RTend<-gspinbutton(from=0,to=80,by=0.01,value=80,con=gfDecon)
				gbDeconEICPlot<-gbutton(text="plot",cont=gfDecon,handler=plotDeconEIC)
			gfSpec<- gframe(text="spectrum",con=gfEIC,horizontal=FALSE)
				gbutton(text="plot",cont=gfSpec,handler=plotSpec)
	visible(ggViz,FALSE)
	
	
	
	
	ggQuantitation<-ggroup(horizontal=FALSE,cont=w)
		ggQviz<-ggroup(cont=ggQuantitation)
	
		gfQ <- gframe(text="visualize Quantitation results",con=ggQviz,horizontal=FALSE)
		gbData<-gbutton(text="Load concerntration file...",cont=gfQ,handler=LoadConcerntrationData)
		gbData<-gbutton(text="Load compound list...",cont=gfQ,handler=LoadmatchResult)
		gbData<-gbutton(text="Load aligned Data...",cont=gfQ,handler=LoadAlignedFeatureQ)
		gbData<-gbutton(text="Load biomarker Data...",cont=gfQ,handler=LoadbioMarker)

		glabel(text="Quantitation Type:",con=gfQ,anchor=c(-1,1))
		quantitation.type<-gradio(horizontal=TRUE,c("Height","Area"),con=gfQ,,anchor=c(-1,1))
		
		gbutton(text="select compound...",cont=gfQ,handler=SelectCompoundtoQuantify)
		gbutton(text="select QMass...",cont=gfQ,handler=SelectQmass)
		gtQMz<-gedit(text="",con=gfQ)
		gbutton(text="plot",cont=gfQ,handler=plotQuantitation)
		gbutton(text="plot all selected compounds",cont=gfQ,handler=hPlotQuantitationAll)
				
	visible(ggQuantitation,FALSE)
	
	ggMVA<-ggroup(horizontal=FALSE,cont=w)
		ggMVAviz<-ggroup(cont=ggMVA)
	
		gfMVA <- gframe(text="MVA",con=ggMVAviz,horizontal=FALSE)
		gbutton(text="Load processed Data",cont=gfMVA,handler=function(h,...){params$JobName<<-svalue(gtJobName);preprocessAnalysisData(params)})
		
		gbutton(text="PCA",cont=gfMVA,handler=plotPCA)
		gbutton(text="PLSDA",cont=gfMVA,handler=plotPLSDA)
		gbutton(text="Clustering Analysis",cont=gfMVA,handler=plotClusterA)
#		gbData<-gbutton(text="Load compound list...",cont=gfQ,handler=LoadmatchResult)
#		gbData<-gbutton(text="Load aligned Data...",cont=gfQ,handler=LoadAlignedFeature)
#		gbData<-gbutton(text="Load biomarker Data...",cont=gfQ,handler=LoadbioMarker)
#
#		glabel(text="Quantitation Type:",con=gfQ,anchor=c(-1,1))
#		quantitation.type<-gradio(horizontal=TRUE,c("Height","Area"),con=gfQ,,anchor=c(-1,1))
#		
#		gbutton(text="select compound...",cont=gfQ,handler=SelectCompoundtoQuantify)
#		gbutton(text="select QMass...",cont=gfQ,handler=SelectQmass)
#		gtQMz<-gedit(text="",con=gfQ)
#		gbutton(text="plot",cont=gfQ,handler=plotQuantitation)
#		gbutton(text="plot all selected compounds",cont=gfQ,handler=hPlotQuantitationAll)
#				
	visible(ggMVA,FALSE)
	
	sb <- gstatusbar("ready to run", cont=w)

#}
