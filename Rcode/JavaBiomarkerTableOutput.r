codeDir<-"code"
source(paste(codeDir,"pipeline.r",sep="/"))
Arg <- commandArgs(TRUE)
params<-readParaFromCmd(Arg)
# the arg[2] save the quantitation type from command
# the arg[1] is the dbfile
params$QuantitationType <- Arg[2]
data_select(params)