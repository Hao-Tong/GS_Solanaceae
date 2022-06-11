# t-test for significant difference
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
dir <- "D:/04_tomatoperpper/"
setwd(dir)

##########################################################################
# multiple traits
mt <- read.table("GS_MT/pepper/multi_trait.csv",sep=",",header=F)

# single trait predictability 
stdata <- read.table("GS/pepper/rrBLUP/pearson_tagsnp.10.csv",head=F,sep=",")

# multiple trait predictability 
# t test
pvalall <- NULL
for (i in 1:nrow(mt)){
  
  mtdata <- read.table(paste0("GS_MT/pepper/results/pearson_tagsnp.10_",i,".csv"),head=F,sep=",")
  mtn <- as.numeric(mt[i,])
  mtn <- mtn[!is.na(mtn)]
  stdatai <- as.matrix(stdata[,mtn])
  
  for (j in 1:ncol(mtdata)){
    stgs <- stdatai[1:150,j]
    mtgs <- mtdata[1:150,j]
    test <- t.test(stgs,mtgs,alternative="two.sided",na.action=na.omit,paired=TRUE)
    pval <- test$p.value
    pvalall <- rbind(pvalall,pval)
  }
} 

write.table(pvalall,"GS_MT/pepper/ttest.csv",sep=",",row.names=F,col.names=F)

##########################################################################
