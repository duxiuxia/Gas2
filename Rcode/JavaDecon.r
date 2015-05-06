codeDir<-"code"
source(paste(codeDir,"pipeline.r",sep="/"))
params<-readParaFromCmd(commandArgs(TRUE))
#params<-readParaFromCmd("~/testing/KQC2/KQC.db")
#params<-readParaFromCmd("~mike/testing/QC15/QC15.db")
#params<-readParaFromCmd("~yni1/Test/GCTestforADAP/UrineStds.db")
params$codeDir<-codeDir
fileInfo<-read.csv(paste(params$WorkDir,"output/",params$JobName,"_RT.csv",sep=""))
params$delaytime<-fileInfo$firstRT[1]
params$ScanInterval<-fileInfo$ScanInterval[1]
#params$NonUmassVec<-c(1:50,72,73,147,221)
#params$NonUmassVec<-c(1:50,73,147,221)
#parTICPeakpicking(params,denoised=F)
#parEICpeakpicking(params)
deconvolution(params,isDistParallel=F,clustingType="h",isDataParalell=T)



