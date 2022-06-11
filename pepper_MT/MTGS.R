# Multi-trait Genome Selection using sommer
# Prediction accuracy and coincidence from cross-validations
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
#dir <- "D:/04_tomatoperpper/GS_MT/pepper/"
dir <- "/winmounts/tong/winhome/mpidir/11.TAGWAS/GS_MT/pepper/"
setwd(dir)

library(sommer)
library(MatrixCorrelation)
options(digits = 15)

##########################################################################
## load data
# genotypic data
myY <- read.table("data/pepper_pheno_all_blup.csv",head=TRUE,sep=",")
pheno <- myY[,-1]

# phenotypic data
myG <- read.table("data/pepper_snp_maf5_filter_num.txt",head=TRUE,sep="\t")
geno <- myG[,-1]

for (p in 1:ncol(pheno)){
  x <- pheno[,p]
  x[is.na(x)] <- mean(x,na.rm=T)
  pheno[,p] <- x
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
for (p in 1:ncol(pheno)){
  x <- pheno[,p]
  pheno[,p] <- range01(x)
}

# cross-validation fold id
idsall <- read.table("fold/foldid.csv",sep=",",header=F)
rr <- 30 #replicate number
f <- 5 #fold number

mt <- read.table("multi_trait.csv",sep=",",header=F)

##########################################################################

#xx <- as.matrix(genotag)
#name <- paste0("tagsnp.",p)
#pheno <- phenomt

MTGS <- function(xx,name,pheno){
  
    Y <- as.matrix(pheno)
    ypred_all1 <- matrix(NA,nrow(pheno),rr)
    ypred_all2 <- matrix(NA,nrow(pheno),rr)
    ypred_all3 <- matrix(NA,nrow(pheno),rr)
    
    for (r in 1:rr){
      
      print(paste("Replicate ",r," Phenotype Group ",t," Done!",sep=""))
      ids <- idsall[,r]
      
      for (p in 1:f){
        trset <- which(ids!=p)
        teset <- which(ids==p)
        xxtr <- xx[trset,]
        xxte <- xx[teset,]
        ytr <- Y[trset,]
        yte <- Y[teset,]
        
        ### multiple trait GS based on sommer package
        xxtr <- as.matrix(xxtr)
        ytr <- data.frame(id=c(1:length(trset)),ytr)
        npheno <- ncol(ytr)-1
        
        # two trait GS
        if (npheno==2){
        colnames(ytr)[-1] <- paste("y",c(1:2),sep="")
        mix.MTGS <- mmer(fixed=cbind(y1,y2)~1,
                           random=~vs(xxtr,Gtc=unsm(npheno)), 
                           rcov=~vs(units,Gtc=unsm(npheno)), getPEV = FALSE, ##the value change to 2 corresponding to 2 traits
                           data=ytr, verbose = FALSE,tolparinv=1e-1)
        #summary(mix.MTGS)
        
        if (length(mix.MTGS)!=0){
        ypred_y1 <- mix.MTGS[["Beta"]][["Estimate"]][1]+xxte%*%mix.MTGS[["U"]][["u:xxtr"]][["y1"]]
        ypred_y2 <- mix.MTGS[["Beta"]][["Estimate"]][2]+xxte%*%mix.MTGS[["U"]][["u:xxtr"]][["y2"]]
        ypred_all1[teset,r] <- ypred_y1
        ypred_all2[teset,r] <- ypred_y2
        }
        }
        
        # three trait GS
        if (npheno==3){
          colnames(ytr)[-1] <- paste("y",c(1:3),sep="")
          mix.MTGS <- mmer(fixed=cbind(y1,y2,y3)~1,
                           random=~vs(xxtr,Gtc=unsm(npheno)), 
                           rcov=~vs(units,Gtc=unsm(npheno)), getPEV = FALSE, ##the value change to 3 corresponding to 3 traits
                           data=ytr, verbose = FALSE,tolparinv=1e-1)
          #summary(mix.MTGS)
          

          if (length(mix.MTGS)!=0){
          ypred_y1 <- mix.MTGS[["Beta"]][["Estimate"]][1]+xxte%*%mix.MTGS[["U"]][["u:xxtr"]][["y1"]]
          ypred_y2 <- mix.MTGS[["Beta"]][["Estimate"]][2]+xxte%*%mix.MTGS[["U"]][["u:xxtr"]][["y2"]]
          ypred_y3 <- mix.MTGS[["Beta"]][["Estimate"]][3]+xxte%*%mix.MTGS[["U"]][["u:xxtr"]][["y3"]]
          ypred_all1[teset,r] <- ypred_y1
          ypred_all2[teset,r] <- ypred_y2
          ypred_all3[teset,r] <- ypred_y3
          }
          }
        
      }
    }
    
    if (npheno==2){
    write.table(ypred_all1,paste("results/predict_",name,"_trait_",t,"_1.csv",sep=""),sep=",",row.names=F,col.names=F)
    write.table(ypred_all2,paste("results/predict_",name,"_trait_",t,"_2.csv",sep=""),sep=",",row.names=F,col.names=F)
    }
    if (npheno==3){
      write.table(ypred_all1,paste("results/predict_",name,"_trait_",t,"_1.csv",sep=""),sep=",",row.names=F,col.names=F)
      write.table(ypred_all2,paste("results/predict_",name,"_trait_",t,"_2.csv",sep=""),sep=",",row.names=F,col.names=F)
      write.table(ypred_all3,paste("results/predict_",name,"_trait_",t,"_3.csv",sep=""),sep=",",row.names=F,col.names=F)
    }
    
}

### function predictability of regression models
ability_regression_MT <- function(name,pheno){
  
  cor1all <- matrix(NA,rr*f,ncol(pheno))
  cor2all <- matrix(NA,rr*f,ncol(pheno))
  match1all <- matrix(NA,rr*f,ncol(pheno))
  match2all <- matrix(NA,rr*f,ncol(pheno))
  
  for (n in 1:ncol(pheno)){
    
    print(paste("Phenotype ",n," Start!",sep=""))
    Ypred <- read.table(paste("results/predict_",name,"_trait_",t,"_",n,".csv",sep=""),head=F,sep=",")
    Y <- as.numeric(pheno[,n])
    
    for (r in 1:rr){
      
      ids <- idsall[,r]
      cor1p <- NULL
      cor2p <- NULL
      match1p <- NULL
      match2p <- NULL
      
      for (p in 1:f){
        teset <- which(ids==p)
        ypred <- Ypred[teset,r]
        yte <- Y[teset]
        
        if (!all(is.na(ypred))) {
        # accuracy
        cor1 <- cor(yte,ypred,method="pearson",use="complete.obs")
        cor2 <- cor(yte,ypred,method="spearman",use="complete.obs")
        cor1p <- c(cor1p,cor1)
        cor2p <- c(cor2p,cor2)
        
        # coincidence
        ytetop <- order(yte,decreasing=T)[1:round(length(yte)*0.3)] #top 30% genotypes
        ypredtop <- order(ypred,decreasing=T)[1:round(length(yte)*0.3)]
        match1 <- sum(ytetop %in% ypredtop)/length(ytetop)
        match1p <- c(match1p,match1)
        
        # coincidence
        ytetop <- order(yte,decreasing=T)[1:round(length(yte)*0.15)] #top 15% genotypes
        ypredtop <- order(ypred,decreasing=T)[1:round(length(yte)*0.15)]
        match2 <- sum(ytetop %in% ypredtop)/length(ytetop)
        match2p <- c(match2p,match2)
        } else {
          cor1p <- c(cor1p,NA)
          cor2p <- c(cor2p,NA)
          match1p <- c(match1p,NA)
          match2p <- c(match2p,NA)
        }
      }
      
      cor1all[(f*r-(f-1)):(f*r),n] <- cor1p
      cor2all[(f*r-(f-1)):(f*r),n] <- cor2p
      match1all[(f*r-(f-1)):(f*r),n] <- match1p
      match2all[(f*r-(f-1)):(f*r),n] <- match2p
    }
    
  }
  
  cor1allm <- rbind(cor1all,colMeans(cor1all,na.rm=T))
  cor2allm <- rbind(cor2all,colMeans(cor2all,na.rm=T))
  match1allm <- rbind(match1all,colMeans(match1all,na.rm=T))
  match2allm <- rbind(match2all,colMeans(match2all,na.rm=T))
  
  # Output: prediction accuracy and coincidence
  write.table(cor1allm,paste("results/pearson_",name,"_",t,".csv",sep=""),sep=",",row.names=F,col.names=F)
  write.table(cor2allm,paste("results/spearman_",name,"_",t,".csv",sep=""),sep=",",row.names=F,col.names=F)
  write.table(match1allm,paste("results/coincidence30_",name,"_",t,".csv",sep=""),sep=",",row.names=F,col.names=F)
  write.table(match2allm,paste("results/coincidence15_",name,"_",t,".csv",sep=""),sep=",",row.names=F,col.names=F)
  
  results <- rbind(colMeans(cor1all,na.rm=T),colMeans(cor2all,na.rm=T),colMeans(match1all,na.rm=T),colMeans(match2all,na.rm=T),
                   apply(cor1all,2,sd,na.rm=T),apply(cor2all,2,sd,na.rm=T),apply(match1all,2,sd,na.rm=T),apply(match2all,2,sd,na.rm=T))
  results <- t(results)
  results <- cbind(colnames(pheno),results)
  colnames(results) <- c("Trait","Mean-Pearson","Spearman","Top30","Top15","SD-Pearson","Spearman","Top30","Top15")
  write.table(results,paste("results/summary_",name,"_",t,".csv",sep=""),sep=",",row.names=F,col.names=T)
  
}

        
##########################################################################
##########################################################################

Rv_MT <- function(name,pheno){
  
    rvall <- matrix(NA,rr*f,1)
    Y <- as.matrix(pheno)
    
    for (r in 1:rr){
      
      ids <- idsall[,r]
      
      Ypred <- matrix(NA,nrow(pheno),ncol(pheno))
      for (n in 1:ncol(pheno)){
        ypred <- read.table(paste("results/predict_",name,"_trait_",t,"_",n,".csv",sep=""),head=F,sep=",")
        Ypred[,n] <- ypred[,r]
      }
      
      rvp <- NULL
      for (p in 1:f){
        teset <- which(ids==p)
        ypred <- Ypred[teset,]
        yte <- Y[teset,]
        
        if (!all(is.na(ypred))) {
          rv <- RV(ypred,yte,center=TRUE)
          rvp <- c(rvp,rv)
        } else {
          rvp <- c(rvp,NA)
        }
      }
      
      rvall[(f*r-(f-1)):(f*r),1] <- rvp

    }

  rvallm <- rbind(rvall,colMeans(rvall,na.rm=T))

  # Output: prediction accuracy and coincidence
  write.table(rvallm,paste("results/Rv_",name,"_",t,".csv",sep=""),sep=",",row.names=F,col.names=F)

  results <- colMeans(rvall,na.rm=T)
  #results <- t(results)
  results <- cbind(paste0("group",t),results)
  colnames(results) <- c("Trait","Rv")
  write.table(results,paste("results/summary_Rv",name,"_",t,".csv",sep=""),sep=",",row.names=F,col.names=T)
  
}


##########################################################################
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

for (t in 1:nrow(mt)){
  
  print(paste("Phenotype ",t," Start!",sep=""))
  
  mtn <- as.numeric(mt[t,])
  mtn <- mtn[!is.na(mtn)]
  phenomt <- as.matrix(pheno[,mtn])
  MTGS(as.matrix(genotag),paste0("tagsnp.",p),phenomt)
  ability_regression_MT(paste0("tagsnp.",p),phenomt)
  Rv_MT(paste0("tagsnp.",p),phenomt)
  
}

##########################################################################
# recalculate predictability
# p <- "10"
# for (t in 1:nrow(mt)){
#   
#   print(paste("Phenotype ",t," Start!",sep=""))
#   
#   mtn <- as.numeric(mt[t,])
#   mtn <- mtn[!is.na(mtn)]
#   phenomt <- as.matrix(pheno[,mtn])
#   #MTGS(as.matrix(genotag),paste0("tagsnp.",p),phenomt)
#   ability_regression_MT(paste0("tagsnp.",p),phenomt)
#   #Rv_MT(paste0("tagsnp.",p),phenomt)
#   
# }
