# Genome Selection using rrBLUP 
# Prediction accuracy and coincidence from cross-validations
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
dir <- "D:/04_tomatoperpper/GS/pepper/"
setwd(dir)

library(caret)
library(pROC)
options(digits = 15)

##########################################################################
## load data
# phenotypic data --TA
myY <- read.table("data/pepper_pheno_all_blup.csv",head=TRUE,sep=",")
rownames(myY) <- myY[,1]
myY <- myY[,-1]

# selected TA traits
#match_CD <- read.table("data/tomato_match_TA_CP_final.csv",head=FALSE,sep=",")
match <- c(1:7,46:47)
pheno <- myY[,match]

# phenotypic data --CD
myY_CD <- read.table("data/Chilipepper.AG.csv",head=TRUE,sep=",")
rownames(myY_CD) <- myY_CD[,1]
myY_CD <- myY_CD[,-1]

# selected CD traits
#match_CD <- read.table("data/tomato_match_TA_CP_final.csv",head=FALSE,sep=",")
#match_CD <- c(6,7,9,13,26)
match_CD <- c(18,25)
classpheno <- myY_CD[,match_CD]
classpheno <- classpheno[rownames(myY),]

# if more than one classes, then assign as NA
for (i in 1:ncol(classpheno)){
  naid <- which(classpheno[,i]%%1!=0)
  classpheno[naid,i] <- NA
}

# genotypic data
#myG <- read.table("data/tomato_snp_maf5_filter_num.txt",head=TRUE,sep="\t")
#geno <- myG[,-1]
#rownames(geno) <- myG[,1]
#geno <- as.matrix(geno[modern,])

# cross-validation fold id
idsall <- read.table("fold/foldid.csv",sep=",",header=F)
rr <- 30 #replicate number
f <- 5 #fold number

##########################################################################
### function predictability of regression models
ability_regression_classification <- function(name){
  
  for (m in 1:ncol(classpheno)){
    
    classy <- as.numeric(classpheno[,m])
    qq <- as.matrix(table(classy))
    qqid <- cumsum(qq) #Cumulative Sums
    #my <- quantile(y,c(0.3,0.7))
    nle <- length(qqid)
    
    accall <- matrix(NA,rr*f,ncol(pheno))
    aucall <- matrix(NA,rr*f,ncol(pheno))
    match1all <- matrix(NA,rr*f,ncol(pheno))
    match2all <- matrix(NA,rr*f,ncol(pheno))
    
    for (n in 1:ncol(pheno)){
      
      Yqq <- rep(0,length(qqid))
      print(paste("Phenotype ",n," in CD ",m," Start!",sep=""))
      Ypred <- read.table(paste("rrBLUP/predict_",name,"_trait_",match[n],".csv",sep=""),head=F,sep=",")
      Y <- as.numeric(pheno[,n])
      Yrank <- sort(Y)
      
      for (g in 1:(length(qqid)-1)){
        Yqq[g+1] <- mean(Yrank[c(qqid[g]:(qqid[g]+1))])
      }  
      Yqq <- c(Yqq,max(Y))
      
      for (r in 1:rr){
        
        ids <- idsall[,r]
        accp <- NULL
        aucp <- NULL
        match1p <- NULL
        match2p <- NULL
        
        for (p in 1:f){
          teset <- which(ids==p)
          ypred <- Ypred[teset,r]
          yte <- Y[teset]
          
          ypred_c <- rep(NA,length(ypred))
          yte_c <- rep(NA,length(yte))
          
          for (g in 1:(length(Yqq)-1)){
            yte_id <- intersect(which(yte>Yqq[g]),which(yte<=Yqq[g+1]))
            yte_c[yte_id] <- paste0("s",g)
            ypred_id <- intersect(which(ypred>Yqq[g]),which(ypred<=Yqq[g+1]))
            ypred_c[ypred_id] <- paste0("s",g)
          }
          yte <- as.factor(yte_c)
          ypred <- as.factor(ypred_c)
          
          naidte <- which(as.numeric(is.na(yte))==1)
          naidpred <- which(as.numeric(is.na(ypred))==1)
          naid <- c(naidte,naidpred)
          if (length(naid)!=0){
            ypred <- ypred[-naid]
            yte <- yte[-naid]
          }
          
          if (nlevels(yte)>=2){
            ## auc
            #predauc <- multiclass.roc(yte,ypredp)
            aucp <- rbind(aucp,NA)
            
            ## accuracy
            predacc <- defaultSummary(data.frame("obs"=yte,"pred"=ypred))
            accp <- rbind(accp,predacc[1])
            #kappap <- rbind(kappap,predacc[2])
            
            # coincidence in top one class
            ytetop <- which(yte==levels(as.factor(yte))[nle])
            ypredtop <- which(ypred==levels(as.factor(yte))[nle])
            if (length(ytetop) != 0){
              match1 <- sum(ytetop %in% ypredtop)/length(ytetop)
            } else {
              match1 <- NA
            }
            match1p <- rbind(match1p,match1)
            
            # coincidence in top two classes
            ytetop <- which(yte==levels(as.factor(yte))[nle-1])
            ypredtop <- which(ypred==levels(as.factor(yte))[nle-1])
            if (length(ytetop) != 0){
              match2 <- sum(ytetop %in% ypredtop)/length(ytetop)
            } else {
              match2 <- NA
            }
            match12 <- mean(c(match1,match2),na.rm=T)
            match2p <- rbind(match2p,match12)
          } else {
            aucp <- rbind(aucp,NA)
            accp <- rbind(accp,NA)
            match1p <- rbind(match1p,NA)
            match2p <- rbind(match2p,NA)
          }
          
        }
        accall[(f*r-(f-1)):(f*r),n] <- accp
        aucall[(f*r-(f-1)):(f*r),n] <- aucp
        match1all[(f*r-(f-1)):(f*r),n] <- match1p
        match2all[(f*r-(f-1)):(f*r),n] <- match2p
      }
      
    }
    
    aucallm <- rbind(aucall,colMeans(aucall,na.rm=T),apply(aucall,2,sd,na.rm=T))
    accallm <- rbind(accall,colMeans(accall,na.rm=T),apply(accall,2,sd,na.rm=T))
    match1allm <- rbind(match1all,colMeans(match1all,na.rm=T),apply(match1all,2,sd,na.rm=T))
    match2allm <- rbind(match2all,colMeans(match2all,na.rm=T),apply(match2all,2,sd,na.rm=T))
    
    # Output: prediction accuracy and coincidence
    write.table(aucallm,paste("rrBLUP_CD/auc_",name,"_CD",m,".csv",sep=""),sep=",",row.names=F,col.names=F)
    write.table(accallm,paste("rrBLUP_CD/accuracy_",name,"_CD",m,".csv",sep=""),sep=",",row.names=F,col.names=F)
    write.table(match1allm,paste("rrBLUP_CD/coincidence1_",name,"_CD",m,".csv",sep=""),sep=",",row.names=F,col.names=F)
    write.table(match2allm,paste("rrBLUP_CD/coincidence2_",name,"_CD",m,".csv",sep=""),sep=",",row.names=F,col.names=F)
    
    results <- rbind(colMeans(aucall,na.rm=T),colMeans(accall,na.rm=T),colMeans(match1all,na.rm=T),colMeans(match2all,na.rm=T),
                     apply(aucall,2,sd,na.rm=T),apply(accall,2,sd,na.rm=T),apply(match1all,2,sd,na.rm=T),apply(match2all,2,sd,na.rm=T))
    results <- t(results)
    results <- cbind(colnames(pheno),results)
    colnames(results) <- c("Trait","Mean-AUC","Accurancy","Top1","Top2","SD-AUC","Accurancy","Top1","Top2")
    write.table(results,paste("rrBLUP_CD/summary_",name,"_CD",m,".csv",sep=""),sep=",",row.names=F,col.names=T)
  }
}
 
##########################################################################
### test for different tag SNP datasets 

p <- "10"
#rrBLUP_pop(as.matrix(genotag),paste0("tagsnp.",p))
name <- paste0("tagsnp.",p)
ability_regression_classification(paste0("tagsnp.",p))

##########################################################################
