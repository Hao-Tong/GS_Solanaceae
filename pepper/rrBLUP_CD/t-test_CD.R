# t-test for significant difference
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
dir <- "D:/04_tomatoperpper/GS/pepper/"
setwd(dir)

##########################################################################
# CD traits
rfdata <- read.table("RF_CD/accuracy_tagsnp.10.csv",sep=",",header=F)
svcdata <- read.table("SVC_CD/accuracy_tagsnp.10.csv",sep=",",header=F)
cddata <- cbind(rfdata,svcdata)

# TA trait 
tadata1 <- read.table("rrBLUP_CD/accuracy_tagsnp.10_CD1.csv",head=F,sep=",")
tadata2 <- read.table("rrBLUP_CD/accuracy_tagsnp.10_CD2.csv",head=F,sep=",")
tadata <- cbind(tadata1,tadata2)

# t test
# trait 1 -- fruit size
cddata.test <- cddata[,c(1,3)]
tadata.test <- tadata1

pvalallall <- NULL
for (i in 1:2){
  cdgs <- cddata.test[c(1:150),i]
  
  pvalall <- NULL
  for (j in 1:9){
    tags <- tadata.test[c(1:150),j]
    test <- t.test(cdgs,tags,alternative="two.sided",na.action=na.omit,paired=TRUE)
    pval <- test$p.value
    pvalall <- rbind(pvalall,pval)
  }
  pvalallall <- cbind(pvalallall,pvalall)
}
colnames(pvalallall) <- c("RF","SVC")
write.table(pvalallall,"rrBLUP_CD/ttest_CD1.csv",sep=",",row.names=F,col.names=T)

# t test
# trait 2 -- fruit yield
cddata.test <- cddata[,c(2,4)]
tadata.test <- tadata2

pvalallall <- NULL
for (i in 1:2){
  cdgs <- cddata.test[c(1:150),i]
  
  pvalall <- NULL
  for (j in 1:9){
    tags <- tadata.test[c(1:150),j]
    test <- t.test(cdgs,tags,alternative="two.sided",na.action=na.omit,paired=TRUE)
    pval <- test$p.value
    pvalall <- rbind(pvalall,pval)
  }
  pvalallall <- cbind(pvalallall,pvalall)
}
colnames(pvalallall) <- c("RF","SVC")
write.table(pvalallall,"rrBLUP_CD/ttest_CD2.csv",sep=",",row.names=F,col.names=T)
 
##########################################################################
