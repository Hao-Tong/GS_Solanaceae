##########################################################################
## add path and packages
dir <- "D:/04_tomatoperpper/00_TA_trait_raw/pepper"
setwd(dir)

options(digits = 15)

##########################################################################
## load data

data <- read.table("TA-h2.csv",head=TRUE,sep=",")[,5]
category <- read.table("trait_group.csv",head=TRUE,sep=",")

data[33] <- NA

group <- category[,2]
group_u <-  unique(group)

outputall <- NULL
for (i in 1:length(group_u)){
  ids <- which(group==group_u[i])
  ave <- mean(data[ids],na.rm=T)
  outputall <- rbind(outputall,ave)
}

rownames(outputall) <- group_u

write.table(outputall,"TA-h2-summary.csv",sep=",",row.names=T,col.names=F)

##########################################################################
