
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
dir <- "D:/04_tomatoperpper/GS/tomato_wild_new"
setwd(dir)

options(digits = 15)

##########################################################################
## load data
# genotypic data
data <- read.table("rrBLUP/00_summary_tagsnp.10_modern_all.csv",head=TRUE,sep=",")[1:47,]
category <- read.table("data/trait_group.csv",head=TRUE,sep=",")

data[28,-1] <- NA

group <- category[,2]
group_u <-  unique(group)

output <- matrix(NA,length(group_u)+1,7*2)
for (i in 1:7){
  datai <- data[,i+1]
  
  for (j in 1:length(group_u)){
    ids <- which(group==group_u[j])
    ave <- mean(datai[ids],na.rm=T)
    sd <- sd(datai[ids],na.rm=T)
    output[j,i] <- ave
    output[j,i+7] <- sd
  }
  output[length(group_u)+1,i] <- mean(datai,na.rm=T)
  output[length(group_u)+1,i+7] <- sd(datai,na.rm=T)
}

rownames(output) <- c(group_u,"All")
write.table(output,"rrBLUP/00_summary_tagsnp.10_modern_all_group.csv",sep=",",row.names=T,col.names=F)

##########################################################################
