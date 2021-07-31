require(HiTC)
require(MASS)
require(Matrix)
require(tidyverse)
require(plyranges)
require(reshape2)
require(rtracklayer)
require(cowplot)
require(ggforce)
require(keras)
require(BSgenome.Hsapiens.UCSC.hg19.masked)
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
get_data <- function(file,AsiSI.GR,mode="AsiSI",dist = c(5000,100e3),normFact = NULL){
  HTC <- loadRData(file)
  mat_DIvAi=intdata(HTC)
  if(is.null(normFact)){
    mat_DIvAi=mat_DIvAi/sum(triu(mat_DIvAi))
  }else{
    facOHTi=sum(mat_DIvAi)/normFact
    mat_DIvAi=mat_DIvAi/facOHTi
  }
  bin.GRi=x_intervals(HTC)
  idxAsiSI=which(countOverlaps(bin.GRi,AsiSI.GR)>0)
  
  tab_DIvAi=summary(mat_DIvAi)
  if(mode=="_noAsiSI"){
    tab_DIvAi=tab_DIvAi[!tab_DIvAi$i%in%idxAsiSI & !tab_DIvAi$j%in%idxAsiSI,]
  }
  tab_DIvAi$dist=(tab_DIvAi$j-tab_DIvAi$i)*binSize
  tab_DIvAi=tab_DIvAi[tab_DIvAi$dist>=dist[1] & tab_DIvAi$dist<dist[2],]
  by_DIvAi=by(tab_DIvAi$x,tab_DIvAi$dist,mean)
  meanCount_DIvAi=as.vector(by_DIvAi)
  meanDist_DIvAi=as.numeric(names(by_DIvAi))
  meanDistKb_DIvAi=meanDist_DIvAi/1e3
  
  res_tibble <- tibble(meanCount = meanCount_DIvAi,meanDistKb=meanDistKb_DIvAi)
  return(list(res_tibble,sum(mat_DIvAi)))
}


seqlens <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 %>% seqlengths()

mes_res <- c(
  # "5kb_100kb"=5e3,
  # "20kb_500kb"=10e3,
  "50kb_1000kb"=10e3,
  # "50kb_5000kb"=100e3,
  # "100kb_5000kb"=100e3,
  "500kb_50000kb"=100e3
)
mes_zoom <- list(
  # "5kb_100kb"=c(1,15),
  "50kb_1000kb"=c(200,250),
  "500kb_50000kb"=c(10000,15000)
)


mes_breaks <- list(
  # "5kb_100kb"=list("Kb",1,c(10,20,50),20),
  "50kb_1000kb"=list("Kb",1,c(100,200,500,1000),100),
  "500kb_50000kb"=list("Mb",1e-3,c(1,2,5,10,20,50)*1e3,1000)
)


Chr.V=glue::glue("chr{c(1,13)}")
ctrlCOND <- "DIvA"
obserOE <- "observed"

DSB174 <- read_bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/ASIsites_hg19_174_clived_IQ1.5.bed")
for(manips in glue::glue("manip{c('A','D')}")){
  pdf(glue::glue("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/heatmap_HiC_D/plot_Distance_chr1_chr13_{manips}.pdf"),height=4,width=19)
  for(one_res in names(mes_res)){
    binSize=mes_res[[one_res]]
    
    sub_zoom <- mes_zoom[[one_res]]
    
    unit_breaks <- mes_breaks[[one_res]][[1]]
    unit_scale <- mes_breaks[[one_res]][[2]]
    unit_breaks_points <- mes_breaks[[one_res]][[3]]
    unit_breaks_linetype <- mes_breaks[[one_res]][[4]]
    
    binSizeTxt=paste0(binSize/1e3,"kb")
    dist=str_extract_all(one_res,"([0-9]+)kb") %>% unlist() %>% str_remove("kb") %>% as.numeric()
    dist <- dist * 1e3
    message(dist)
    if(obserOE == "observed"){
      if(manips == "manipA"){
        HiC.files <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/HTC",pattern=binSizeTxt,full.names=T) %>% str_subset("manipA")
      }else{
        HiC.files <- list.files("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/HTC/observed",pattern=str_c("_",binSizeTxt),full.names=T)
      }
      
      
    }else{
      HiC.files <- list.files("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/HTC/OE",pattern=str_c("_",binSizeTxt),full.names=T)
    }
    
    
    res.DSB <- lapply(Chr.V,function(chrom){
      message(chrom)
      f <- HiC.files[str_detect(HiC.files,str_c("_",chrom,"_"))]
      Replicate <- f %>% map(.%>% basename() %>% str_remove_all("HTC_observed_KR_|_manip[AB]|HTC_(observed|OE)_KR_HiC_D_|_chr[0-9A-Z]+_[0-9]+kb.RData")) %>% unlist()
      names(f) <- Replicate
      # Replicate <- f %>% map(.%>% basename() %>% str_extract("manip[A-Z]+|rep[0-9]+"))
      diva.chr <- DSB174 %>%filter(seqnames ==chrom) %>% 
        anchor_center() %>% mutate(width = 1e6)
      if(length(diva.chr)<2)
        return(NULL)
      
      f_ctrl <- Replicate %>% str_subset(glue::glue("{ctrlCOND}$"))
      f_exp <- Replicate %>% str_subset(glue::glue("{ctrlCOND}$"),negate = T)
      
      
      c(res_CTRL,normFactDIVA) %<-% get_data(file = f[f_ctrl],AsiSI.GR = diva.chr,mode="_noAsiSI",dist = dist)
      res <- lapply(f[f_exp],function(s){
        c(sres_exp,normFactexpe) %<-% get_data(file = s,AsiSI.GR = diva.chr,mode="_noAsiSI",dist = dist,normFact = normFactDIVA)
        
        sres_exp
      })
      
      
      
      res[["CTRL"]] <- res_CTRL
      
      res <- res %>% bind_rows(.id ="Condition")
      
      
      
      res %>% ggplot(aes(x=meanDistKb,y=meanCount,col=Condition)) +
        geom_line() +
        facet_zoom(xy = meanDistKb > sub_zoom[1] & meanDistKb< sub_zoom[2],zoom.size = 1) +
        ylab("HiC contact probability") + xlab(glue::glue("Distance ({unit_breaks})")) + ggtitle(chrom) +
        scale_y_log10() + scale_x_log10(breaks = unit_breaks_points,labels = scales::unit_format(unit = unit_breaks, scale = unit_scale,accuracy=1)) 
      
    })
    cowplot::plot_grid(plotlist = res.DSB) %>% print()
    
  }
  dev.off()
}



setwd("/media/mourad/diskSave/MCF_Toulouse/recherche/LegubeTeam/")
library(Matrix)
library(MASS)
library(HiTC)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(HiCblock)


# LOAD DATA AND RESULTS ----------------------------------------------------------

# Exp
expe="manipA"
obserOE="observed" # "observed" "OE"
mode="_noAsiSI" # "" "_AsiSI" "_noAsiSI"

# Resolution
# binSize=100e3
binSize=10e3
binSizeTxt=paste0(binSize/1e3,"kb")
# dist=c(500e3,50000e3)
dist=c(50e3,1000e3)

# Chromosomes
SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
Chr.V=paste0("chr",c(1:22,"X"))
seqlen=seqlengths(SeqInfo[Chr.V])

# AsiSI
AsiSI.GR=read_bed("/mnt/NAS/DATA/AsiSI/ASIsites_hg19.bed")
AsiSI.GR=resize(AsiSI.GR,fix="center",width=1e6)

# Folder
# dir.create(paste0("results/HiC/distance/res",binSizeTxt,"_",dist[1]/1e3,"kb_",dist[2]/1e3,"kb"))

# Build HiTC per chr

HiC.files <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/HTC",pattern=binSizeTxt,full.names=T) %>% str_subset(expe)



# fileDIvA=paste0("data/HiC_Legube/HTC/HTC_",obserOE,"_KR_DIvA_",expe,"_intrachr_",binSizeTxt,".RData")
# load(fileDIvA)
# HTCL_DIvA=HTCL
# fileOHT=paste0("data/HiC_Legube/HTC/HTC_",obserOE,"_KR_OHT_",expe,"_intrachr_",binSizeTxt,".RData")
# load(fileOHT)
# HTCL_OHT=HTCL
matCoef=NULL
for(i in 1:length(Chr.V)){
  f <- HiC.files[str_detect(HiC.files,str_c("_",i,"_"))]
  # Load data
  HTC_DIvAi=loadRData(f[1])
  mat_DIvAi=intdata(HTC_DIvAi)
  mat_DIvAi=mat_DIvAi/sum(triu(mat_DIvAi))
  bin.GRi=x_intervals(HTC_DIvAi)
  idxAsiSI=which(countOverlaps(bin.GRi,AsiSI.GR)>0)
  
  # Normalize across conditions
  HTC_OHTi=loadRData(f[2])
  mat_OHTi=intdata(HTC_OHTi)
  facOHTi=sum(mat_OHTi)/sum(mat_DIvAi)
  mat_OHTi=mat_OHTi/facOHTi
  
  # Compute slopes
  tab_DIvAi=summary(mat_DIvAi)
  if(mode=="_noAsiSI"){
    tab_DIvAi=tab_DIvAi[!tab_DIvAi$i%in%idxAsiSI & !tab_DIvAi$j%in%idxAsiSI,]
  }
  tab_DIvAi$dist=(tab_DIvAi$j-tab_DIvAi$i)*binSize
  tab_DIvAi=tab_DIvAi[tab_DIvAi$dist>dist[1] & tab_DIvAi$dist<dist[2],]
  by_DIvAi=by(tab_DIvAi$x,tab_DIvAi$dist,mean)
  meanCount_DIvAi=as.vector(by_DIvAi)
  meanDist_DIvAi=as.numeric(names(by_DIvAi))
  meanDistKb_DIvAi=meanDist_DIvAi/1e3
  lm_DIvAi=lm(log10(meanCount_DIvAi)~log10(meanDistKb_DIvAi))
  coef_DIvAi=summary(lm_DIvAi)$coefficients[2,1]
  
  tab_OHTi=summary(mat_OHTi)
  if(mode=="_noAsiSI"){
    tab_OHTi=tab_OHTi[!tab_OHTi$i%in%idxAsiSI & !tab_OHTi$j%in%idxAsiSI,]
  }
  tab_OHTi$dist=(tab_OHTi$j-tab_OHTi$i)*binSize
  tab_OHTi=tab_OHTi[tab_OHTi$dist>dist[1] & tab_OHTi$dist<dist[2],]
  by_OHTi=by(tab_OHTi$x,tab_OHTi$dist,mean)
  meanCount_OHTi=as.vector(by_OHTi)
  meanDist_OHTi=as.numeric(names(by_OHTi))
  meanDistKb_OHTi=meanDist_OHTi/1e3
  lm_OHTi=lm(log10(meanCount_OHTi)~log10(meanDistKb_OHTi))
  coef_OHTi=summary(lm_OHTi)$coefficients[2,1]
  
  # Test difference of slopes
  dataDiffDist=data.frame(meanDistKb=c(meanDistKb_DIvAi,meanDistKb_OHTi),meanCount=c(meanCount_DIvAi,meanCount_OHTi),
                          expe=c(rep("DIvA",length(meanDistKb_DIvAi)),rep("OHT",length(meanDistKb_OHTi))))
  lmi=lm(log10(meanCount)~log10(meanDistKb)*expe,data=dataDiffDist)
  coef_diffi=summary(lmi)$coefficients[4,]
  
  matCoefi=c(coefDIvA=coef_DIvAi,coefOHT=coef_OHTi,coef_diff=coef_diffi)
  matCoef=rbind(matCoef,matCoefi) 
  
  # Plot
  file_ploti=paste0("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/plot_distance_",expe,"_",Chr.V[i],"_",dist[1]/1e3,"kb_",dist[2]/1e3,"kb",mode,".pdf")
  pdf(file_ploti,6,6)
  plot(meanDistKb_DIvAi,meanCount_DIvAi,log="xy",col="blue",type="l",xlab="Distance (kb)",ylab="Hi-C contacts",
       main=paste0("Slope DIvA=",round(coef_DIvAi,3),", Slope OHT=",round(coef_OHTi,3),", p=",round(coef_diffi[4],5)))
  lines(meanDistKb_OHTi,meanCount_OHTi,col="red")
  legend("topright",legend=c("DIvA","OHT"),col=c("blue","red"),lty=1)
  dev.off()
  
  rm(data.MatSymi,chri.GR)
  print(i)
}

matCoef=data.frame(chr=Chr.V,matCoef)

filemat=paste0("results/HiC/distance/res",binSizeTxt,"_",dist[1]/1e3,"kb_",dist[2]/1e3,"kb","/matCoef_distance_",expe,"_",dist[1]/1e3,"kb_",dist[2]/1e3,"kb",mode,".csv")
write.table(matCoef,file=filemat,sep='\t',row.names=F,quote=F)





