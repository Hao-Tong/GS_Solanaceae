
# Contact: tong@mpimp-golm.mpg.de
 
##########################################################################
## add path and packages
dir <- "D:/04_tomatoperpper/GS/pepper"
setwd(dir)

options(digits = 15)

##########################################################################
## load data
# genotypic data
data <- read.table("rrBLUP/00_summary_all_sets.csv",head=TRUE,sep=",")[1:47,]
category <- read.table("data/trait_group.csv",head=TRUE,sep=",")

data[33,-1] <- NA

group <- category[,2]
group_u <-  unique(group)

outputall <- NULL
for (p in 1:4){
  
cols <- c((8*(p-1)+1):(p*8))
output <- matrix(NA,length(group_u)+1,8*2)
for (i in 1:8){

  datai <- data[,cols[i]+1]
  
  for (j in 1:length(group_u)){
    ids <- which(group==group_u[j])
    ave <- mean(datai[ids],na.rm=T)
    sd <- sd(datai[ids],na.rm=T)
    output[j,i] <- ave
    output[j,i+8] <- sd
  }
  output[length(group_u)+1,i] <- mean(datai,na.rm=T)
  output[length(group_u)+1,i+8] <- sd(datai,na.rm=T)
}

rownames(output) <- c(group_u,"All")

outputall <- rbind(outputall,output)

}

write.table(outputall,"rrBLUP/00_summary_all_sets_group.csv",sep=",",row.names=T,col.names=F)
##########################################################################
