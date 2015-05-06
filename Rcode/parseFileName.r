parseFileName<-function(filePath,isPath=T)
{
	if(isPath)
	{
		fileName<-strsplit(filePath,"/")[[1]]#parse the path
		fileName<-fileName[length(fileName)]#get last element as filename
	}else
	{
		fileName<-filePath
	}
	
	sepFilename<-strsplit(fileName,".",fixed=TRUE)[[1]]
	extName<-sepFilename[length(sepFilename)]
	fileName<-substr(fileName,0,nchar(fileName)-nchar(extName)-1)#remove extention name
	return(fileName)
}
