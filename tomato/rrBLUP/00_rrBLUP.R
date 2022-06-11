# Genome Selection using rrBLUP 
# Prediction accuracy and coincidence from cross-validations
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
dir <- "D:/04_tomatoperpper/GS/tomato/"
#dir <- "/winmounts/tong/winhome/mpidir/11.TAGWAS/GS_all/tomato/"
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

# heterozygous rate
heteroall <- NULL
for (i in 1:nrow(geno)){
  datai <- geno[i,]
  heteroi <- length(which(datai==0.5))
  heteroall <- c(heteroall,heteroi)
}
write.table(heteroall/ncol(geno),"data/heterozygous_all.csv",sep=",",row.names=F,col.names=F)

# only modern accessions
modern <- read.table("data/modern.txt",header=F,sep="\t")[,1]
mset <- match(myY[,1],modern)
mset <- which(mset != "NA")

pheno <- pheno[mset,]
geno <- as.matrix(geno[mset,])

# cross-validation fold id
idsall <- read.table("fold/foldid.csv",sep=",",header=F)
rr <- 30 #replicate number
f <- 5 #fold number

##########################################################################
### function rrBLUP model
rrBLUP <- function(xx,name){
  
  for (n in 1:ncol(pheno)){
      
    ypred_all <- matrix(NA,nrow(pheno),rr)
    print(paste("Phenotype ",n," Start!",sep=""))
    Y <- as.numeric(pheno[,n])
    
    for (r in 1:rr){
      
      print(paste("Replicate ",r," Phenotype ",n," Done!",sep=""))
      ids <- idsall[,r]
      
      for (p in 1:f){
        trset <- which(ids!=p)
        teset <- which(ids==p)
        xxtr <- xx[trset,]
        xxte <- xx[teset,]
        ytr <- Y[trset]
        sol <- mixed.solve(ytr,Z=xxtr,K=NULL,SE=F)
        ucoef <- as.matrix(sol$u)
        ypred <- rep(sol$beta,nrow(xxte))+as.numeric(xxte %*% ucoef)
        ypred_all[teset,r] <- ypred
      }
    }
    
    write.table(ypred_all,paste("rrBLUP/predict_",name,"_trait_",n,".csv",sep=""),sep=",",row.names=F,col.names=F)
    
  }
}

### function predictability of regression models
ability_regression <- function(name){
  
  cor1all <- matrix(NA,rr*f,ncol(pheno))
  cor2all <- matrix(NA,rr*f,ncol(pheno))
  match1all <- matrix(NA,rr*f,ncol(pheno))
  match2all <- matrix(NA,rr*f,ncol(pheno))
  
  for (n in 1:ncol(pheno)){
    
    print(paste("Phenotype ",n," Start!",sep=""))
    Ypred <- read.table(paste("rrBLUP/predict_",name,"_trait_",n,".csv",sep=""),head=F,sep=",")
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
  write.table(cor1allm,paste("rrBLUP/pearson_",name,".csv",sep=""),sep=",",row.names=F,col.names=F)
  write.table(cor2allm,paste("rrBLUP/spearman_",name,".csv",sep=""),sep=",",row.names=F,col.names=F)
  write.table(match1allm,paste("rrBLUP/coincidence30_",name,".csv",sep=""),sep=",",row.names=F,col.names=F)
  write.table(match2allm,paste("rrBLUP/coincidence15_",name,".csv",sep=""),sep=",",row.names=F,col.names=F)
  
  results <- rbind(colMeans(cor1all,na.rm=T),colMeans(cor2all,na.rm=T),colMeans(match1all,na.rm=T),colMeans(match2all,na.rm=T),
                   apply(cor1all,2,sd,na.rm=T),apply(cor2all,2,sd,na.rm=T),apply(match1all,2,sd,na.rm=T),apply(match2all,2,sd,na.rm=T))
  results <- t(results)
  results <- cbind(colnames(pheno),results)
  colnames(results) <- c("Trait","Mean-Pearson","Spearman","Top30","Top15","SD-Pearson","Spearman","Top30","Top15")
  write.table(results,paste("rrBLUP/summary_",name,".csv",sep=""),sep=",",row.names=F,col.names=T)

}

##########################################################################
## predict agronomic traits from snp 
# filter with non-variance SNPs
nonvar <- apply(geno,2,var)
nonid <- which(nonvar==0)
geno_non <- geno[,-nonid]

rrBLUP(as.matrix(geno_non),"allsnp")
ability_regression("allsnp")

##########################################################################
### test for different tag SNP datasets 

percentage <- c("05","1","10","25","50","75","100")

ntag_allall <- NULL
for (p in percentage){

ntag_all <- NULL
tagname_all <- NULL
for (i in 1:12){
  allname <- read.table(paste0("tag/chr",i,".100.SNPid.txt"),head=F,sep="\t")[,1]
  tagid <- read.table(paste0("tag/chr",i,".",p,".tagSNP.txt"),head=F,sep="\t")[,1]
  ntag <- length(tagid)
  ntag_all <- c(ntag_all,ntag)
  tagid <- tagid+1
  tagname <- allname[tagid]
  tagname_all <- c(tagname_all,tagname)
}

write.table(ntag_all,paste0("tag/tagsnp_number.",p,".csv"),sep=",",row.names=F,col.names=F)
ntag_allall <- rbind(ntag_allall,ntag_all)

# filter tag snps
tagname_all <- gsub("-",".",tagname_all)
tagname_all <- paste0("X",tagname_all)
genotag <- geno[,tagname_all]

# run the GS model
rrBLUP(as.matrix(genotag),paste0("tagsnp.",p))
ability_regression(paste0("tagsnp.",p))

}

write.table(t(ntag_allall),"tag/tagsnp_number_all.csv",sep=",",row.names=F,col.names=F)

##########################################################################
### summary all results

percentage <- c("05","1","10","25","50","75","100")
tagsnps <- paste0("tagsnp.",percentage)
allname <- c(tagsnps,"allsnp")

resultsall_final <- matrix(NA,ncol(pheno),4*8)
for (m in 1:8){
  name <- allname[m]
  data <- read.table(paste0("rrBLUP/summary_",name,".csv"),head=T,sep=",")[,2:5]
  resultsall_final[,seq(m,4*8,by=8)] <- as.matrix(data)
}

resultsall_final <- cbind(colnames(pheno),resultsall_final)
colnames(resultsall_final) <- c("Trait",rep("Pearson",8),rep("Spearman",8),rep("Top30%",8),rep("Top15%",8))
write.table(resultsall_final,"rrBLUP/00_summary_all_sets.csv",sep=",",row.names=F,col.names=T)
## without sd
##########################################################################
