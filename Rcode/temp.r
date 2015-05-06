##########################################################
#plot old QC 12.22-12.30
####################################################
curTICPk$lboundInd<-(12.22-delaytime)/ScanInterval+1
curTICPk$rboundInd<-(12.30-delaytime)/ScanInterval+1

write.csv(curPeakList,file="~ktsuttle/peakpicking/pklist73.csv")

for(i in 1:nrow(modelPkList))
{
	curID<-rownames(modelPkList[i,])
	modelPkList[i,]$gss<-GaussianDisimilarity(allprofiles[[curID]]$int)
	modelPkList[i,]$shrp<-Sharpness(allprofiles[[curID]]$int)
	
}

############################################################
#check quantitation of the decon results
########################################################################
specList<-readMSP2Spec(filename="~/testing/KQC2/output/testNist.txt",T)
evaluateIdentification(refSpec=lib,inputSpec=specList,RT_Tolerance=20,minSpecSimilarity=650,"KQC",withRT=T,params)
matchresult<-read.csv("~/testing/KQC2/output/evaluationResult_KQC.csv")
alignedResult<-read.csv("~/testing/KQC2/output/KQC_aligned.csv")


qresult<-NULL
rsd<-NULL
for(i in 1:nrow(matchresult))
{
	
	curQmass<-matchresult$Qmass[i]
	curclustID<-as.integer(do.call(rbind,strsplit(do.call(rbind,strsplit(as.character(matchresult$componentName[i]),split=" "))[,2],"t"))[,2])
	
	localResult<-subset(alignedResult,clustID==curclustID&mz==curQmass,select=c("alignedET","mz","pkArea","FileID"))
	rsd<-c(rsd,sd(localResult$pkArea)/mean(localResult$pkArea))
	qresult<-rbind(qresult,localResult)
	
}
write.csv(rsd,file="rsd.csv")

write.csv(qresult,file="qresult.csv")







		
####################################################################################

Sharpness(allprofiles[["10591"]]$int)

oxoproline<-allprofiles[[rownames(subset(CandidatePeaklist,mz==236))[1]]]
xx<-allprofiles[[rownames(subset(CandidatePeaklist,mz==245))]]
plot((oxoproline$ET-1)*ScanInterval+delaytime,oxoproline$int,type="l")
 plot((xx$ET-1)*ScanInterval+delaytime,xx$int,type="l")

r["92835","40099"]


aa<-read.csv("~/testing/KQC2/output/decon/QC_2_1mdlpklist_bck.csv")$x
bb<-EICpeaklist[as.character(aa),]

bb$gss1<-90
bb$shrp1<-0
for(i in 1:nrow(bb))
{
	curPk<-bb[i,]
	
	yy<-vecInt[(curPk$lboundInd+curPk$offset):(curPk$rboundInd+curPk$offset)]
	bb[i,]$gss1<-GaussianSimilarity1(yy)
	bb[i,]$shrp1<-Sharpness1(yy)
	
}

colVec<-rainbow(nrow(bb))
names(colVec)<-bb$mz
labelVec<-paste("mz:",bb$mz," ,gss:",format(bb$gss,digits=2)," ,gss1:",format(bb$gss1,digits=2)," ,shrp:",format(bb$shrp,digits=2)," ,shrp1:",format(bb$shrp1,digits=2))
dataPoints<-subset(bb,select=c("pkInd","Intensity"))
dataPoints$pkInd<-(dataPoints$pkInd-1)*ScanInterval+delaytime
playwith({
		plot(x=0,type="n",xlim=(c(min(bb$lboundInd),max(bb$lboundInd))-1)*ScanInterval+delaytime,ylim=c(min(bb$Intensity),max(bb$Intensity)))
		for(i in 1:nrow(bb))
		{
			curPk<-bb[i,]
			points(((curPk$lboundInd:curPk$rboundInd)-1)*ScanInterval+delaytime,vecInt[(curPk$lboundInd+curPk$offset):(curPk$rboundInd+curPk$offset)],type="l",col=colVec[i])	
		}
		
		},data.points=dataPoints,labels=labelVec)

plot(allprofiles[["20105"]]$int,type="l")

yy<-xx$int
plot(yy)
GaussianDisimilarity1(yy)
GaussianDisimilarity(yy)


coneproj=function(y,amat){
	n=length(y);m=length(amat)/n
	sm=1e-8;h=1:m<0;obs=1:m;check=0
	delta=-amat;b2=delta%*%y
	if(max(b2)>sm){
		i=min(obs[b2==max(b2)])
		h[i]=TRUE
	}else{check=1;theta=1:n*0}
	while(check==0){
		xmat=matrix(delta[h,],ncol=n)
		a=solve(xmat%*%t(xmat))%*%xmat%*%y
		if(min(a)<(-sm)){
			avec=1:m*0;avec[h]=a
			i=min(obs[avec==min(avec)])
			h[i]=FALSE;check=0
		}else{
			check=1
			theta=t(xmat)%*%a
			b2=delta%*%(y-theta)
			if(max(b2)>sm){
				i=min(obs[b2==max(b2)])		
				h[i]=TRUE;check=0
			}
		}
	}
	bhat=y-theta
	bhat
}



write.csv(m1,row.names=FALSE,file="m.txt")
write.csv(S,row.names=FALSE,file="s.txt")

fr <- function(x,M,S) {
#	M<-read.csv("m.txt")
#	S<-read.csv("s.txt")
	m<-0
	for(i in 1:ncol(S))
	{
		m<-m+x[i]*S[,i]
	}
 	sum((M-m)^2)
	
}



m1<-S[,1]*0.5+S[,2]*1


pseudoinverse(S[,1:3])%*%M
optim(par=rep(0,ncol(S)),fn=fr, gr = NULL,M=m1,S=S,lower = 0, method="L-BFGS-B")$par



plot(M,type="l")
points(S[,1]*120,type="l",col="red")
points(S[,2]*2,type="l",col="blue")
points(S[,3]*2,type="l",col="green")
points(S[,4]*300,type="l",col="pink")

componentResults<-read.csv("~/testing/QC15/output/decon/control_1Decon.csv")
pkResults<-read.csv("~/testing/QC15/output/peakpicking/control_1PeakList.csv")

subset(pkResults,abs((pkInd-1)*ScanInterval+delaytime-12.255)<0.1&curMz==75)$Intensity


subset(componentResults,windowID==37&compoundID==2&mz==75)$Intensity
subset(componentResults,windowID==37&compoundID==4&mz==75)$Intensity


Rsd1<-read.csv("~/testing/QC15/output/1.0/RsdResult_S1.csv",sep="\t")
SbioMarker<-read.csv("~/testing/QC15/output/1.0/QC15_S_BiomarkerTable.csv")

matchedResult1<-read.csv("~/testing/QC15/output/1.0/evaluation_QcCompound.csv",sep=";")
matchedClustIDs<-unique(as.integer(as.character(subset(matchedResult1,score>0)$componentName)))

Rsd1<-subset(Rsd1,clustID%in%matchedClustIDs)
SbioMarker<-subset(SbioMarker,clustID%in%matchedClustIDs)
write.csv(SbioMarker,file="~/testing/QC15/output/1.0/selectedSbioMarker.csv")

mean(Rsd1$rsdVec)

Rsd2<-read.csv("~/testing/QC15/output/2.0/RsdResult_S2.csv",sep="\t")
SbioMarker<-read.csv("~/testing/QC15/output/2.0/QC15_S_BiomarkerTable.csv")

matchedResult2<-read.csv("~/testing/QC15/output/2.0/evaluation_QcCompound.csv",sep=";")
matchedClustIDs<-unique(as.integer(as.character(subset(matchedResult2,score>0)$componentName)))

Rsd2<-subset(Rsd2,clustID%in%matchedClustIDs)
SbioMarker<-subset(SbioMarker,clustID%in%matchedClustIDs)
write.csv(SbioMarker,file="~/testing/QC15/output/2.0/selectedSbioMarker.csv")

mean(Rsd2$rsdVec)
Rsd2$std/Rsd2$meanVec


inFilePath <- DataFilelist[[fileindex]]
fileName<-parseFileName(inFilePath)
denoisedTICFile<-paste(WorkDir,"/output/denoised_",fileName,"_TIC.cdf",sep="")
ncid <- open.ncdf(denoisedTICFile)
RawDenoisedTIC <- get.var.ncdf(ncid, varid="total_intensity")
close.ncdf(ncid)
remove(ncid)	

ncid <- open.ncdf(inFilePath,write=TRUE)
RawTIC <-  get.var.ncdf(ncid, varid="total_intensity")


fileName<-parseFileName(inFilePath)
EICFile<-paste(WorkDir,"/output/EIC/",fileName,"EIC.cdf",sep="")
cat(fileName,"reading EIC data...\n")



ncid <- open.ncdf(EICFile)
vecInt <- get.var.ncdf(ncid, varid="intVec")
mzVec<-get.var.ncdf(ncid, varid="mzVec")
totalscan<-length(vecInt)/length(mzVec)

##############
#smoothing
################
cat(fileName,"smoothing EIC data...\n")

DenoisedTIC <- 0
getEIC<-function(vecInt,mzVec,mz)
{
	mzInd <- which(mzVec==mz)
	startInd<-(mzInd-1)*totalscan+1
	endInd<-startInd+totalscan-1
	curVecInt <- vecInt[startInd:endInd]
	curVecInt
}




	

library(playwith)
playwith(
		{ 
			plot((1:length(DenoisedTIC))*params$ScanInterval+params$delaytime,vecInt,type = "l", col = "blue")
			points((1:length(DenoisedTIC))*params$ScanInterval+params$delaytime,RawDenoisedTIC,type="l",col = "red")
#			points((1:length(DenoisedTIC))*params$ScanInterval+params$delaytime,RawTIC,type="l")
#			points((1:length(DenoisedTIC))*params$ScanInterval+params$delaytime,getEIC(vecInt,mzVec,mz=133),type="l",col = "purple")
			
		})


time1<-Sys.time()

time2<-Sys.time()
time2-time1