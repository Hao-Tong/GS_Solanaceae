# Genome Selection using rrBLUP 
# Prediction accuracy and coincidence from cross-validations
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
#dir <- "D:/04_tomatoperpper/GS/tomato_wild/"
dir <- "/winmounts/tong/winhome/mpidir/11.TAGWAS/GS_all/tomato_wild/"
setwd(dir)

library(rrBLUP)
options(digits = 15)

##########################################################################
## load data
# genotypic data
myY <- read.table("data/tomato_pheno_all_blup.csv",head=TRUE,sep=",")
pheno <- myY[,-1]

# phenotypic data
myG <- read.table("data/tomato_snp_maf5_filter_num.txt",head=TRUE,sep="\t")
geno <- myG[,-1]

# cross-validation fold id
#idsall <- read.table("fold/foldid.csv",sep=",",header=F)
#rr <- 30 #replicate number
#f <- 5 #fold number

# modern and wild accession ids
modern <- read.table("data/modern.txt",header=F,sep="\t")[,1]
wild <- read.table("data/wild.txt",header=F,sep="\t")[,1]

nrepf <- 150

##########################################################################
### function rrBLUP model
nwild <- c(0,7,5*c(2:6))
rrBLUP_pop <- function(xx,name){
  
  for (s in c(0:6)){
    
    print(paste("Wild ",s," Start!",sep=""))
    
    if (s == 0){
      nrep <- 1
      trset_name <- matrix(modern,ncol=1)
      teset_name <- matrix(wild,ncol=1)
    } else {
      nrep <- nrepf
      set_all <- read.table(paste("fold/testid_wild_",nwild[s+1],".csv",sep=""),head=F,sep=",")
      trset_name <- set_all[-c(1:(length(wild)-nwild[s+1])),]
      teset_name <- set_all[c(1:(length(wild)-nwild[s+1])),]
    }
    
    ypred_all <- matrix(NA,(length(wild)-nwild[s+1]),nrep)
    for (n in 1:(ncol(pheno))){
        
      print(paste("Wild ",s," Phenotype ",n," Start!",sep=""))
      Y <- as.numeric(pheno[,n])
        
      for (r in 1:nrep){
          
        trset <- match(myY[,1],trset_name[,r])
        teset <- match(myY[,1],teset_name[,r])
        trset <- which(trset != "NA")
        teset <- which(teset != "NA")
          
        xxtr <- xx[trset,]
        xxte <- xx[teset,]
        ytr <- Y[trset]
        sol <- mixed.solve(ytr,Z=xxtr,K=NULL,SE=F)
        ucoef <- as.matrix(sol$u)
        ypred <- rep(sol$beta,nrow(xxte))+as.numeric(xxte %*% ucoef)
        ypred_all[,r] <- ypred
      }
      
      write.table(ypred_all,paste("rrBLUP/predict_",name,"_trait_",n,"_wild_",nwild[s+1],".csv",sep=""),sep=",",row.names=F,col.names=F)
    }

  }
}

### function predictability of regression models
ability_regression_pop <- function(name){
  
  cor1all <- matrix(NA,nrepf,ncol(pheno))
  cor2all <- matrix(NA,nrepf,ncol(pheno))
  match1all <- matrix(NA,nrepf,ncol(pheno))
  match2all <- matrix(NA,nrepf,ncol(pheno))
  resultsall <- NULL
  for (s in c(0:6)){
    
    print(paste("Wild ",s," Start!",sep=""))

    if (s == 0){
      nrep <- 1
      #trset_name <- matrix(modern,ncol=1)
      teset_name <- matrix(wild,ncol=1)
    } else {
      nrep <- nrepf
      set_all <- read.table(paste("fold/testid_wild_",nwild[s+1],".csv",sep=""),head=F,sep=",")
      #trset_name <- set_all[-c(1:(length(wild)-nwild[s+1])),]
      teset_name <- set_all[c(1:(length(wild)-nwild[s+1])),]
    }
    
    for (n in 1:(ncol(pheno))){
      
      Ypred <- read.table(paste("rrBLUP/predict_",name,"_trait_",n,"_wild_",nwild[s+1],".csv",sep=""),head=F,sep=",")
      Y <- as.numeric(pheno[,n])
      cor1p <- NULL
      cor2p <- NULL
      match1p <- NULL
      match2p <- NULL
      
      for (r in 1:nrep){
        
        #trset <- match(myY[,1],trset_name[,r])
        teset <- match(myY[,1],teset_name[,r])
        #trset <- which(trset != "NA")
        teset <- which(teset != "NA")
        yte <- Y[teset]
        ypred <- Ypred[,r]

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
      }
      
      cor1all[,n] <- cor1p
      cor2all[,n] <- cor2p
      match1all[,n] <- match1p
      match2all[,n] <- match2p
    }
    
    cor1allm <- rbind(cor1all,colMeans(cor1all,na.rm=T))
    cor2allm <- rbind(cor2all,colMeans(cor2all,na.rm=T))
    match1allm <- rbind(match1all,colMeans(match1all,na.rm=T))
    match2allm <- rbind(match2all,colMeans(match2all,na.rm=T))
    
    # Output: prediction accuracy and coincidence
    write.table(cor1allm,paste("rrBLUP/pearson_",name,"_wild_",nwild[s+1],".csv",sep=""),sep=",",row.names=F,col.names=F)
    write.table(cor2allm,paste("rrBLUP/spearman_",name,"_wild_",nwild[s+1],".csv",sep=""),sep=",",row.names=F,col.names=F)
    write.table(match1allm,paste("rrBLUP/coincidence30_",name,"_wild_",nwild[s+1],".csv",sep=""),sep=",",row.names=F,col.names=F)
    write.table(match2allm,paste("rrBLUP/coincidence15_",name,"_wild_",nwild[s+1],".csv",sep=""),sep=",",row.names=F,col.names=F)
    
    results <- rbind(colMeans(cor1all,na.rm=T),colMeans(cor2all,na.rm=T),colMeans(match1all,na.rm=T),colMeans(match2all,na.rm=T))
    results <- t(results)
    results <- cbind(colnames(pheno),results)
    colnames(results) <- c("Trait","Pearson","Spearman","Top30%","Top15%")
    write.table(results,paste("rrBLUP/summary_",name,"_wild_",nwild[s+1],".csv",sep=""),sep=",",row.names=F,col.names=T)
    
    resultsall <- cbind(resultsall,results[,-1])
    }
  
    resultsall_final <- matrix(NA,nrow(resultsall),ncol(resultsall))
    for (i in 1:4){
      resultsall_final[,c((7*(i-1)+1):(7*i))] <- resultsall[,seq(i,4*7,by=4)]
    }
    resultsall_final <- cbind(colnames(pheno),resultsall_final)
    colnames(resultsall_final) <- c("Trait",rep("Pearson",7),rep("Spearman",7),rep("Top30%",7),rep("Top15%",7))
    write.table(resultsall_final,paste("rrBLUP/00_summary_",name,"_all.csv",sep=""),sep=",",row.names=F,col.names=T)

}

##########################################################################
### test for different tag SNP datasets 
### test the different number of wild aceessions

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
rrBLUP_pop(as.matrix(genotag),paste0("tagsnp.",p))
ability_regression_pop(paste0("tagsnp.",p))

##########################################################################
