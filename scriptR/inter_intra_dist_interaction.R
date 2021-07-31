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


TAD_borders <- "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC/TAD/TopDom/TAD_topdom_DIvA_manipA_hg19_50kb_ws50.bed" %>% read_bed()
# TAD_borders <- "/media/HDD_ROCHER/PROJET_INGE/HiC_Coline/TADborder_topdom_DIvA_manipA_hg19_50kb_ws10_strong_t-0.05.bed" %>% read_bed()



mes_res <- c(
  "5kb_2000kb"=5e3
)




Chr.V=glue::glue("chr{c(1,16)}")
ctrlCOND <- "DIvA"
obserOE <- "observed"

DSB174 <- read_bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/ASIsites_hg19_174_clived_IQ1.5.bed")
for(manips in glue::glue("manip{c('A','D')}")){
  pdf(glue::glue("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/heatmap_HiC_D/plot_Distance_chr1_chr13_{manips}.pdf"),height=4,width=19)
  for(one_res in names(mes_res)){
    binSize=mes_res[[one_res]]
    
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