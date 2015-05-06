dotProduct<- function(x, y) 
	{
		if(length(x)<=10||length(y)<=10)
		{
			dotprod<-0
		}else
		{
			dotprod<-(x%*%y)/((x%*%x)^0.5*(y%*%y)^0.5)
			dotprod<-as.numeric(dotprod[1])
		}
		
		angle<-acos(dotprod)*180/pi
		if(is.na(angle))
			angle<-90
		angle
		#dotprod
	}
f<- function(x, y) 
	{
		(1-cor(x,y))*1000
	}
distCal<-function (x) 
{
    
       
            ncy <- ncx <- ncol(x)
            if (ncx == 0) 
                stop("'x' is empty")
            r <- matrix(0, nrow = ncx, ncol = ncy)
            for (i in seq_len(ncx)) {
                for (j in seq_len(i)) {
                  x2 <- x[, i]
                  y2 <- x[, j]
 				r[i, j] <- ifelse(i==j,0,dotProduct(x2, y2))
                }
            }
            r <- r + t(r) - diag(diag(r))
            rownames(r) <- colnames(x)
            colnames(r) <- colnames(x)
            r
        
       
   
}
