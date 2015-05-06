split.matrix<-function(x,f) {
  #print('processing matrix')
  v=lapply(
      1:dim(x)[[2]]
      , function(i) {
        base:::split.default(x[,i],f)#the difference is here
      }
      )

  w=lapply(
      seq(along=v[[1]])
      , function(i) {
        result=do.call(
            cbind
            , lapply(v,
                function(vj) {
                  vj[[i]]
                }
                )
            )
        colnames(result)=colnames(x)
        return(result)
      }
      )
  names(w)=names(v[[1]])
  return(w)
}