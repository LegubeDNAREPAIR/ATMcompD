require(glue)
require(tidyverse)
#bash code
# expe="CTRLDIvA CTRLOHT siSCC1DIvA siSCC1OHT"
# obserOE="observed"
# res=25000
# reskb="25kb"
# for exp in $expe; do for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do  java -jar programs/juicer/juicer_tools.jar dump $obserOE KR /home/rochevin/Documents/PROJET_INGE/HiC_Coline/REVISIONS/hic_files/${exp}.hic $chr $chr BP $res data/KR/dump_${obserOE}_KR_${exp}_chr${chr}_${reskb}.txt; gzip data/KR/dump_${obserOE}_KR_${exp}_chr${chr}_${reskb}.txt; done; done

options(scipen=999)
dirHiC <- "~/ownCloud/Documents/HiC_D"
expe <- list.files(dirHiC,recursive = T,pattern=".hic")
my_res <- c(1000000,500000,250000,100000,50000,25000,10000,5000)
names(my_res) <- str_c((my_res/1000),"kb")
obserOE <- c("observed","OE")

juicer <- "/media/HDD_ROCHER/PROJET_INGE/HiC_Coline/programs/juicer/juicer_tools.jar"

for(exp in expe){
  expe_name <- dirname(exp)
  message(exp)
  for(myObs in obserOE){
    message(myObs)
    for(res in names(my_res)){
      message(res)
      mclapply(c(1:22,"X"),function(chr){
        my_out <- glue("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/dump_{myObs}_KR_{expe_name}_chr{chr}_{res}")
        cmd <- glue("java -jar {juicer} dump {myObs} KR {dirHiC}/{exp} {chr} {chr} BP {my_res[[res]]} {my_out}.txt; gzip {my_out}.txt")
        message(cmd)
        system(cmd)
      },mc.cores=8)
    }
  }
}

# WORKING DIRECTORY ----------------------------------------------------

setwd("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/")
library(Matrix)
library(MASS)
library(HiTC)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)


# Chromosomes
SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
Chr.V=paste0("chr",c(1:22,"X"))
seqlen=seqlengths(SeqInfo[Chr.V])


# Build HiTC per chr
for(obserOE in c("OE","observed")){
  for(binSizeTxt in names(my_res)){
    binSize=my_res[[binSizeTxt]]
    for(k in 1:length(expe)){
      HTCl=list()
      for(i in 1:length(Chr.V)){
        
        # Load data
        chrendi=seqlengths(SeqInfo[Chr.V[i]])
        binendi=ceiling(chrendi/binSize)
        
        fileCounti=paste0("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/dump_",obserOE,"_KR_",dirname(expe[k]),"_",Chr.V[i],"_",binSizeTxt,".txt.gz")
        datai=as.matrix(fread(paste0("zcat ",fileCounti),sep='\t',header=F))
        if(nrow(datai) == 0)
          next
        datai=datai[!is.na(datai[,3]),]
        datai=rbind(datai,c((binendi-1)*binSize,(binendi-1)*binSize,0))
        if(obserOE=="observed"){
          datai[,3]=round(datai[,3])
        }else{
          datai[,3]=round(datai[,3],2)
        }
        
        data.Mati=sparseMatrix(i=(datai[,1]/binSize)+1,j=(datai[,2]/binSize)+1,x=datai[,3])
        data.MatSymi=data.Mati+t(data.Mati)
        diag(data.MatSymi)=diag(data.MatSymi)/2
        rm(data.Mati,datai)
        
        starti=seq(1,chrendi,by=binSize)
        endi=c(seq(binSize,chrendi,by=binSize),chrendi)
        chri.GR=GRanges(Chr.V[i],IRanges(starti,endi))
        names(chri.GR)=paste0("bin",1:length(chri.GR))
        
        HTC=HTCexp(data.MatSymi,chri.GR,chri.GR)
        HTCl[[i]]=HTC
        fileouti=paste0("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/HTC/",obserOE,"/HTC_",obserOE,"_KR_",dirname(expe[k]),"_",Chr.V[i],"_",binSizeTxt,".RData")
        save(HTC,file=fileouti)
        
        rm(data.MatSymi,chri.GR)
      }
      #HTCL=HTClist(HTCl)
      # fileouti=paste0("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/HiC/HTC/HTC_",obserOE,"_KR_",expe[k],"_intrachr_",binSizeTxt,".RData")
      # save(HTCL,file=fileouti)
    }
    
  }
}


#plot HiC

















