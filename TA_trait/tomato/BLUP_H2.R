### BLUP estimation for phenotypic traits across multiple environments
### Variance components of G, E, and GxE
### Broad-sense heritability estimation across multiple environments
### Contact: tonghao0605@gmail.com

dir <- ("D:/04_tomatoperpper/00_TA_trait_raw/tomato_all/")
setwd(dir)
rm(list=ls())

library(lme4)

#######################################################################
## combine and check
#######################################################################
# data set 
#dataid <- read.table("tomato_pheno_all_id.csv",sep=",",header=T)[,2]
data1 <- read.table("tomato_rep1.csv",sep=",",header=T)
data2 <- read.table("tomato_rep2.csv",sep=",",header=T)
data3 <- read.table("tomato_rep3.csv",sep=",",header=T)
data4 <- read.table("tomato_rep4.csv",sep=",",header=T)

data1a <- read.table("../tomato_wild/tomato_wild_rep1.csv",sep=",",header=F)
data2a <- read.table("../tomato_wild/tomato_wild_rep2.csv",sep=",",header=F)
data3a <- read.table("../tomato_wild/tomato_wild_rep3.csv",sep=",",header=F)
data4a <- read.table("../tomato_wild/tomato_wild_rep4.csv",sep=",",header=F)
 
colnames(data1a) <- colnames(data1)
write.table(rbind(data1,data1a),"tomato_all_rep1.csv",row.names=F,col.names=T,sep=",",quote = F)
colnames(data2a) <- colnames(data2)
write.table(rbind(data2,data2a),"tomato_all_rep2.csv",row.names=F,col.names=T,sep=",",quote = F)
colnames(data3a) <- colnames(data3)
write.table(rbind(data3,data3a),"tomato_all_rep3.csv",row.names=F,col.names=T,sep=",",quote = F)
colnames(data4a) <- colnames(data4)
write.table(rbind(data4,data4a),"tomato_all_rep4.csv",row.names=F,col.names=T,sep=",",quote = F)

#######################################################################
# data set 
#dataid <- read.table("tomato_pheno_all_id.csv",sep=",",header=T)[,2]
data1 <- read.table("tomato_all_rep1_final.csv",sep=",",header=T)
data2 <- read.table("tomato_all_rep2_final.csv",sep=",",header=T)
data3 <- read.table("tomato_all_rep3_final.csv",sep=",",header=T)
data4 <- read.table("tomato_all_rep4_final.csv",sep=",",header=T)

#######################################################################
### input data format: row for each line, column for each trait

dataall <- rbind(data1,data2,data3,data4)
dataall <- dataall[,-c(1:2)]

n <- nrow(data1) #lines number
p <- ncol(data1) #phenotype number
m <- 4 #environment number

### prepare input data format 
lineid <- rep(1:n,m)
locid <- c(sort(rep(1:m,n)))

xa <- cbind(lineid,locid,dataall)
colnames(xa) <- c("LINE","LOC",1:47)
write.table(xa,"TA-all-rep.csv",sep=',',quote=F,row.names=F)

#######################################################################
#######################################################################

datamean_all <- NULL
for (i in 1:n){
  ids <- which(xa[,1]==i)
  datai <- xa[ids,-c(1:2)]
  datameani <- colMeans(datai,na.rm = T)
  datamean_all <- rbind(datamean_all,datameani)
}
rownames(datamean_all) <- c(1:n)

write.table(datamean_all,"TA-all-rep_mean.csv",sep=',',quote=F,row.names=T,col.names=F)

#######################################################################
#######################################################################

##########################################################
### BLUP with G-by-E
line.blup <- c(1:n)
heritability <- NULL

x <- read.table("TA-all-rep.csv",sep=",",header=T)
colnames(x)[1:2] <- c("LINE","LOC")

for(i in 3:ncol(x)){
  
  print(i-2)
  
  # control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4),
  #                     check.nobs.vs.nlev = "ignore",
  #                     check.nobs.vs.rankZ = "ignore",
  #                     check.nobs.vs.nRE="ignore")
  # varcomp <- lmer(x[,i]~(1|LINE)+(1|LOC)+(1|LINE:LOC),data=x,control=control) 
  # #isSingular(varcomp, tol = 1e-4)
  
  varcomp <- lmer(x[,i]~(1|LINE)+(1|LOC),data=x) 
                  # control=lmerControl(check.nobs.vs.nlev = "ignore",
                  # check.nobs.vs.rankZ = "ignore",
                  # check.nobs.vs.nRE="ignore"))
  
  var.trans <- lme4::VarCorr(varcomp)
  var <- data.frame(Groups=c('LINE','LOC','Residual'),
                    Variance=c(as.numeric(var.trans$LINE),as.numeric(var.trans$LOC),attr(var.trans,'sc')^2),check.names=F)
  #residual standard deviation is stored as attribute "sc"
  Gvar<-as.numeric(as.character(var$Variance))[var$Groups%in%'LINE']
  Evar<-as.numeric(as.character(var$Variance))[var$Groups%in%'LOC']
  evar<-as.numeric(as.character(var$Variance))[var$Groups%in%'Residual']
  
  h2 <- c(Gvar,Evar,evar,Gvar/(Gvar+Evar+evar))
  heritability <- rbind(heritability,h2)
  
  f <- fixef(varcomp)
  r <- ranef(varcomp)$LINE
  blup <- f+r
  line.blup <- cbind(line.blup,blup)
}
colnames(line.blup) <- c('line',names(x)[-c(1:2)])
heritability <- cbind(c(colnames(x)[-c(1:2)]),heritability)
colnames(heritability) <- c("trait","G variance","E variance","residual","h2")

write.table(line.blup,"TA-BLUP.csv",row.names=F,sep=",")
write.table(heritability,"TA-h2.csv",row.names=F,col.names=T,sep=",")

##########################################################
  