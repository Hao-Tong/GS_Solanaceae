# Genome Selection using rrBLUP 
# Prediction accuracy and coincidence from cross-validations
# Contact: tong@mpimp-golm.mpg.de

##########################################################################
## add path and packages
dir <- "D:/04_tomatoperpper/GS/tomato_wild/"
#dir <- "/winmounts/tong/winhome/mpidir/11.TAGWAS/GS_all/tomato_wild/"
setwd(dir)

##########################################################################
## load data
# modern and wild accession ids
modern <- read.table("data/modern.txt",header=F,sep="\t")[,1]
#wild <- read.table("data/wild.txt",header=F,sep="\t")[,1]
wild_all <- read.table("data/wild_species_final.csv",header=F,sep=",")
wild <- wild_all[,1]
wild_species <- wild_all[,2]

nrepf <- 150

# sampling data with wild accessions
nwild <- c(7,5*c(2:6))

for (i in c(1:6)){
  
    trset_all <- matrix(NA,(length(modern)+nwild[i]),nrepf)
    teset_all <- matrix(NA,(length(wild)-nwild[i]),nrepf)
    for (r in 1:nrepf){
      
      # sample from each wild species
      # ensure one accession of one wild species
      trset.wild_s_all <- NULL
      for (s in 1:7){
        wild_s <- wild[which(wild_species==s)]
        trset.wild_s <- sample(wild_s,1,replace=FALSE)
        trset.wild_s_all <- c(trset.wild_s_all,trset.wild_s)
      }
      
      # sample additional wild accessions from remains
      wild_remain <- setdiff(wild,trset.wild_s_all)
      trset.wild_add <- sample(wild_remain,nwild[i]-7,replace=FALSE)
      
      # combine with modern accessions
      trset_name <- c(modern,trset.wild_s_all,trset.wild_add)
      teset_name <- setdiff(wild,c(trset.wild_s_all,trset.wild_add))
      trset_all[,r] <- trset_name
      teset_all[,r] <- teset_name
    }
    
    # first testing set, then training set!! 
    set_all <- rbind(teset_all,trset_all)
    write.table(set_all,paste("fold/testid_wild_",nwild[i],".csv",sep=""),sep=",",row.names=F,col.names=F)
    
  }
   
##########################################################################
