#####################################################################
dir <- "D:/04_tomatoperpper/GS/pepper/fold"
setwd(dir)

f <- 5 #folds number
nn <- 162 #sample number
n <- floor(nn/f)
id <- rep(c(1:f),n)
if (nn > length(id)){
	id <- c(id,c(1:(nn-length(id))))
}

rr <- 50 #replicate number
foldall <- NULL
for (r in 1:rr){
	ids <- sample(id)
	foldall <- cbind(foldall,ids)
}

write.table(foldall,"foldid.csv",sep=",",row.names=F,col.names=F)
