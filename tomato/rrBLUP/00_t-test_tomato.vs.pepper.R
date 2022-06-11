# t-test for significant difference
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
dir <- "D:/04_tomatoperpper/GS/"
setwd(dir)

##########################################################################
# tomato traits
todata <- read.table("tomato/rrBLUP/pearson_allsnp.csv",sep=",",header=F)

# pepper trait 
pedata <- read.table("pepper/rrBLUP/pearson_allsnp.csv",sep=",",header=F)

# t test
# each traits
pvalall <- NULL
for (i in 1:47){
  togs <- todata[c(1:150),i]
  pegs <- pedata[c(1:150),i]
  test <- t.test(togs,pegs,alternative="two.sided",na.action=na.omit,paired=TRUE)
  pval <- test$p.value
  pvalall <- rbind(pvalall,pval)
}
pvalall <- cbind(c(1:47),pvalall)
write.table(pvalall,"tomato/rrBLUP/ttest_trait.csv",sep=",",row.names=F,col.names=F)

##########################################################################
# t test
# each category
category <- read.table("tomato/data/trait_group.csv",head=TRUE,sep=",")

group <- category[,2]
group_u <-  unique(group)

pvalall <- NULL
for (j in 1:length(group_u)){
  ids <- which(group==group_u[j])
  togs <- todata[1:150,ids]
  pegs <- pedata[1:150,ids]
  togsc <- c(as.matrix(togs))
  pegsc <- c(as.matrix(pegs))
  test <- t.test(togsc,pegsc,alternative="two.sided",na.action=na.omit,paired=TRUE)
  pval <- test$p.value
  pvalall <- rbind(pvalall,pval)
}
 
togsall <- todata[1:150,1:47]
pegsall <- pedata[1:150,1:47]
togsallc <- c(as.matrix(togsall))
pegsallc <- c(as.matrix(pegsall))
test <- t.test(togsallc,pegsallc,alternative="two.sided",na.action=na.omit,paired=TRUE)
pval <- test$p.value
pvalall <- rbind(pvalall,pval)

pvalall <- cbind(c(1:11),pvalall)
write.table(pvalall,"tomato/rrBLUP/ttest_group.csv",sep=",",row.names=F,col.names=F)

##########################################################################