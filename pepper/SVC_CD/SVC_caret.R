#!/usr/bin/Rscript
  
# Genome Selection model SVC
# All models are implemented with R package caret
# Prediction accuracy and coincidence from cross-validations
# Contact: tong@mpimp-golm.mpg.de
 
##########################################################################
## add path and packages
dir <- "D:/04_tomatoperpper/GS/pepper/"
setwd(dir)

library(caret)
#library(pls)
#library(randomForest)
library(kernlab)
#library(RSNNS)
#library(fastDummies)
library(pROC)
library(UBL)
options(digits = 15)

##########################################################################
## load data
myY_CD <- read.table("data/Chilipepper.AG.csv",head=TRUE,sep=",")
rownames(myY_CD) <- myY_CD[,1]
myY_CD <- myY_CD[,-1]
match_CD <- c(18,25)
classpheno <- myY_CD[,match_CD]

## select the same accession as TA set
myY_id <- read.table("data/pepper_pheno_all.csv",head=TRUE,sep=",")[,1]
classpheno <- classpheno[myY_id,]

# if more than one classes in one cell, then assign as NA 
# (i.e. the number must be integer)
for (i in 1:ncol(classpheno)){
  naid <- which(classpheno[,i]%%1!=0)
  classpheno[naid,i] <- NA
}

# genotypic data
myG <- read.table("data/pepper_snp_maf5_filter_num.txt",head=TRUE,sep="\t")
geno <- myG[,-1]
rownames(geno) <- myG[,1]
geno <- as.matrix(geno)

# cross-validation fold id
idsall <- read.table("fold/foldid.csv",sep=",",header=F)
rr <- 30 #replicate number
f <- 5 #fold number

##########################################################################
##########################################################################
### function support vector machine model
SVC <- function(xx,name){
  
  for (n in 1:ncol(classpheno)){
    
    print(paste("Phenotype ",n," Start!",sep=""))
    Y <- classpheno[,n]
    Y[which(Y != "NA")] <- paste0("s",Y[which(Y != "NA")])
    nle <- nlevels(as.factor(Y))
    
    ypred_all <- matrix(NA,nrow(classpheno),rr)
    ypredp_all <- matrix(NA,nrow(classpheno),rr*nle)
    
    for (r in 1:rr){
      
      ids <- idsall[,r]
      
      for (p in 1:f){
        
        print(paste("Phenotype ",n," Replicate ",r,".",p," Start!",sep=""))
        
        trset <- which(ids!=p)
        teset <- which(ids==p)
        xxtr <- xx[trset,]
        xxte <- xx[teset,]
        ytr <- as.factor(Y[trset])
        
        # remove the NA data accessions
        naid <- which(as.numeric(is.na(ytr))==1)
        if (length(naid)!=0){
          xxtr <- xxtr[-naid,]
          ytr <- ytr[-naid]
        }
        
        trdata <- data.frame(ytr,xxtr)
        #table(trdata$ytr)
        #newtrdata <- smote(ytr~., trdata, perc.over=20, perc.under=20)
        newtrdata <- SmoteClassif(ytr~., trdata, k=1, C.perc="balance")
        #table(newtrdata$ytr)
        
        #folds <- createMultiFolds(ytr,k=5,times=1)
        tcontr <- trainControl(method="repeatedcv",number=5,repeats=1,classProbs=TRUE)
        svcfit <- caret::train(ytr~., data=newtrdata, method="svmRadial", metric="Accuracy",
          na.action=na.pass, trControl=tcontr, preProc = c("center","scale"), tuneLength=15)	
        #The tuneLength parameter is set to pick 15 arbitrary values for the C.
        #use its default method of calculating an analytically derived estimate for sigma.
      
        #ytr <- newtrdata[,1]
        #xxtr <- newtrdata[,-1]
        #rffit1 <- caret::train(x=xxtr, y=ytr, method="rf", metric="Accuracy",
        #		trControl=tcontr, tuneGrid=expand.grid(.mtry=sqrt(ncol(xxtr))),ntree=1000)
        
        ypredp <- predict(svcfit,newdata=xxte,type="prob")
        ypred <- predict(svcfit,newdata=xxte)
        
        ypred_all[teset,r] <- as.character(ypred)
        
        if (ncol(ypredp)==nle){
          ypredp_all[teset,c((nle*(r-1)+1):(nle*r))] <- as.matrix(ypredp)
          #colnames(ypredp_all)[c((nle*(r-1)+1):(nle*r))] <- colnames(ypredp)
        } else {
          # in the case less predicted class than the raw data 
          mid <- levels(as.factor(Y)) %in% colnames(ypredp)
          mid <- levels(as.factor(Y))[which(as.numeric(mid)==0)]
          ypredpp <- cbind(as.matrix(ypredp),rep(0,length(teset)))
          colnames(ypredpp) <- c(colnames(ypredp),mid) 
          ypredp_all[teset,c((nle*(r-1)+1):(nle*r))] <- ypredpp[,levels(as.factor(Y))]
          #colnames(ypredp_all) <- c(colnames(ypredp),mid)       
        }
        
      }
    }
    
    write.table(ypred_all,paste("SVC_CD/predict_",name,"_trait_",n,".csv",sep=""),sep=",",row.names=F,col.names=F)
    write.table(ypredp_all,paste("SVC_CD/predict_probability_",name,"_trait_",n,".csv",sep=""),sep=",",row.names=F,col.names=F)
    
  }
}

### function predictability of classification models
ability_classification <- function(name){
  
  accall <- matrix(NA,rr*f,ncol(classpheno))
  aucall <- matrix(NA,rr*f,ncol(classpheno))
  match1all <- matrix(NA,rr*f,ncol(classpheno))
  match2all <- matrix(NA,rr*f,ncol(classpheno))
  
  for (n in 1:ncol(classpheno)){
    
    print(paste("Phenotype ",n," Start!",sep=""))
    Ypred <- read.table(paste("SVC_CD/predict_",name,"_trait_",n,".csv",sep=""),head=F,sep=",")
    Ypredp <- read.table(paste("SVC_CD/predict_probability_",name,"_trait_",n,".csv",sep=""),head=F,sep=",")
    Y <- classpheno[,n]
    Y[which(Y != "NA")] <- paste0("s",Y[which(Y != "NA")])
    nle <- nlevels(as.factor(Y))
    
    for (r in 1:rr){
      
      ids <- idsall[,r]
      accp <- NULL
      aucp <- NULL
      match1p <- NULL
      match2p <- NULL
      
      for (p in 1:f){
        teset <- which(ids==p)
        yte <- as.factor(Y[teset])
        ypred <- as.factor(Ypred[teset,r])
        ypredp <- Ypredp[teset,c((nle*(r-1)+1):(nle*r))]
        colnames(ypredp) <- levels(as.factor(Y))
        
        naid <- which(as.numeric(is.na(yte))==1)
        if (length(naid)!=0){
          ypred <- ypred[-naid]
          ypredp <- ypredp[-naid,]
          yte <- yte[-naid]
        }
        
        if (nlevels(yte)>=2){
          ## auc
          predauc <- multiclass.roc(yte,ypredp)
          aucp <- rbind(aucp,as.numeric(predauc$auc))
          
          ## accuracy
          predacc <- defaultSummary(data.frame("obs"=yte,"pred"=ypred))
          accp <- rbind(accp,predacc[1])
          #kappap <- rbind(kappap,predacc[2])
          
          # coincidence in top one class
          ytetop <- which(yte==levels(as.factor(Y))[nle])
          ypredtop <- which(ypred==levels(as.factor(Y))[nle])
          if (length(ytetop) != 0){
            match1 <- sum(ytetop %in% ypredtop)/length(ytetop)
          } else {
            match1 <- NA
          }
          match1p <- rbind(match1p,match1)
          
          # coincidence in top two classes
          ytetop <- which(yte==levels(as.factor(Y))[nle-1])
          ypredtop <- which(ypred==levels(as.factor(Y))[nle-1])
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
  
  accallm <- rbind(accall,colMeans(accall,na.rm=T))
  aucallm <- rbind(aucall,colMeans(aucall,na.rm=T))
  match1allm <- rbind(match1all,colMeans(match1all,na.rm=T))
  match2allm <- rbind(match2all,colMeans(match2all,na.rm=T))
  
  # Output: prediction accuracy and coincidence
  write.table(aucallm,paste("SVC_CD/auc_",name,".csv",sep=""),sep=",",row.names=F,col.names=F)
  write.table(accallm,paste("SVC_CD/accuracy_",name,".csv",sep=""),sep=",",row.names=F,col.names=F)
  write.table(match1allm,paste("SVC_CD/coincidence1_",name,".csv",sep=""),sep=",",row.names=F,col.names=F)
  write.table(match2allm,paste("SVC_CD/coincidence2_",name,".csv",sep=""),sep=",",row.names=F,col.names=F)
  
  results <- rbind(colMeans(aucall,na.rm=T),colMeans(accall,na.rm=T),colMeans(match1all,na.rm=T),colMeans(match2all,na.rm=T),
                   apply(aucall,2,sd,na.rm=T),apply(accall,2,sd,na.rm=T),apply(match1all,2,sd,na.rm=T),apply(match2all,2,sd,na.rm=T))
  results <- t(results)
  results <- cbind(colnames(classpheno),results)
  colnames(results) <- c("Trait","Mean-AUC","Accurancy","Top1","Top2","SD-AUC","Accurancy","Top1","Top2")
  write.table(results,paste("SVC_CD/summary_",name,".csv",sep=""),sep=",",row.names=F,col.names=T)
  
}

##########################################################################
### test for different tag SNP datasets 

p <- "10"
tagname_all <- NULL
for (i in 1:12){
  allname <- read.table(paste0("tag/chr",i,".100.SNPid.txt"),head=F,sep="\t")[,1]
  tagid <- read.table(paste0("tag/chr",i,".",p,".tagSNP.txt"),head=F,sep="\t")[,1]
  #ntag <- length(tagid)
  #ntag_all <- c(ntag_all,ntag)
  tagid <- tagid+1
  tagname <- allname[tagid]
  tagname_all <- c(tagname_all,tagname)
}

# filter tag snps
tagname_all <- gsub("-",".",tagname_all)
tagname_all <- paste0("X",tagname_all)
genotag <- geno[,tagname_all]

# run the GS model
SVC(as.matrix(genotag),paste0("tagsnp.",p))
ability_classification(paste0("tagsnp.",p))

##########################################################################
# recalculate predictability
p <- "10"
ability_classification(paste0("tagsnp.",p))

##########################################################################


