# convert a vector to a string seperated by comma
vector2string <- function(x) {
	len <- length(x)
	tmp <- paste("", x[1],sep="")
	for(i in 2:len) {
		tmp <- paste(tmp, x[i], sep=", ")
	}
	tmp <- paste(tmp, "", sep="")
}

tab2matrix <- function(tab=NULL){
	intv.mat <- c()
	for(i in 1:dim(tab)[1]) {
		for(j in 1:dim(tab)[2])  {
			if(tab[i,j] >0)  {
				tmp.row <- c(rownames(tab)[i], colnames(tab)[j], tab[i,j])
				intv.mat <- cbind(intv.mat, as.numeric(tmp.row))
			}
		}
	}
	intv.mat <- t(intv.mat)    
	intv.mat <- cbind(c(1:dim(intv.mat)[1]), intv.mat)
	intv.mat
}