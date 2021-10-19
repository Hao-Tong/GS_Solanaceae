
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
dir <- "D:/04_tomatoperpper/GS_MT/pepper/results/"
setwd(dir)

##########################################################################

resultsall <- NULL
for (i in 1:7){
  
  data <- read.table(paste0("summary_tagsnp.10_",i,".csv"),head=F,sep=",")
  resultsall <- rbind(resultsall,data[-1,])
  
}

write.table(resultsall,"summary_tagsnp.10_all.csv",sep=",",row.names=F,col.names=F)

##########################################################################
