########################
# Read db file
########################

readParaSQLite<-function(paraDB)
{
	
	# create a SQLite instance and create one connection.
	m <- dbDriver("SQLite")
	con <- dbConnect(m, dbname = paraDB)
	# query
	rs <- dbGetQuery(con, "select * from options where pname='File'")
	# clean up
	dbDisconnect(con)
	
	params<-NULL
	
	for(i in 1:nrow(rs))
	{
		paraName<-rs$pname[i]
		paraValue<-rs$pvalue[i]
		expr<-paste("params$DataFiles$",paraName,i,"=paraValue",sep="")
		eval(parse(text=expr))
	}
	
	con <- dbConnect(m, dbname = paraDB)
	# query
	rs <- dbGetQuery(con, "select * from options where pname<>'File'")
	# clean up
	dbDisconnect(con)
	for(i in 1:nrow(rs))
	{
		paraName<-rs$pname[i]
		paraValue<-rs$pvalue[i]
		expr<-paste("params$",paraName,"=paraValue",sep="")
		eval(parse(text=expr))
	}
	
	
	return(params)
	
}

getDataFileFullPath<-function(params)
{
	DataFilelist=params$DataFiles
	for(i in 1:length(DataFilelist))
		DataFilelist[[i]]<-paste(params$DataDir,DataFilelist[[i]],sep="")
	DataFilelist
}

readParaFromCmd<-function(args)
{
	library(RSQLite)
#	dbFile<-"~/testing/KQC2/KQC.db"
	dbFile<-args[1] # the first arg is dbFile name (includes path)
	params<-readParaSQLite(dbFile)
	params$DataFiles<-getDataFileFullPath(params)
	params$clustType="MPI"
	params$dbFile<-dbFile
	params
	
}