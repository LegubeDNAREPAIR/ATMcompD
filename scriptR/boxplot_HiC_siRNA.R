#FROM Heatma_HiC_data_DSB.R
require(HiTC)
require(tidyverse)
require(plyranges)
require(reshape2)
require(keras)
require(rtracklayer)
require(cowplot)
require(BSgenome.Hsapiens.UCSC.hg19.masked)
require(regioneR)
seqlens <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 %>% seqlengths()
#Functions

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


get_sub_matrix_Hic <- function(file,bed,binsize){
  HTC <- loadRData(file)
  rows.chr <- x_intervals(HTC)
  val1 <- HTC %>% intdata %>% sum
  # HTC <- normPerReads(HTC)
  
  overlaps <- as(findOverlaps(bed,rows.chr), "List")
  overlaps.bin <- lapply(overlaps,length) %>% unlist()
  overlaps <- overlaps[overlaps.bin == binsize]
  
  res <- lapply(1:length(overlaps),function(i){
    x <- overlaps[[i]]
    if(as.character(strand(bed[i])) == "-"){
      x <- rev(x)
    }
    intdata(HTC)[x,x]
  })
  return(list(res,val1))
}

get_nb_read_by_chr <- function(file){
  HTC <- loadRData(file)
  HTC %>% intdata %>% sum
}

#ADD binned 5kb genomic region with no damage and lot of damage (excluding damages)
##USING DSB
bless80 = read_bed("/mnt/NAS/DATA/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
HR =  read_bed("/mnt/NAS/DATA/AsiSI/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed")
NHEJ =  read_bed("/mnt/NAS/DATA/AsiSI/BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed")
Random80 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/80random.bed")
uncut <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/ASIsites_hg19_without_bless80.bed")


obserOE = "observed" # observed, OE


#boxplot

MES_DSB <- list(
  "bless80"=bless80,
  "HR"=HR,
  "NHEJ"=NHEJ,
  "uncut"=uncut,
  "Random80"=Random80
)



mes_windows <- c(
  "100kb"=100000
  # "50kb"=50000
  # "10kb"=10000
  # "5kb"=5000
)

mes_limites <- list(
  "100kb"=c(1500000)
  # "50kb"=c(500000)
  # "10kb"=c(1000000)
  # "5kb"=c(500000,250000)
)
res.boxplot <- lapply(names(MES_DSB),function(DSB.n){
  message(DSB.n)
  DSB174 <- MES_DSB[[DSB.n]]
  
  
  
  lapply(names(mes_windows),function(m.w){
    message(m.w)
    sw <- mes_windows[[m.w]]
    
    
    if(obserOE == "observed"){
      HiC.files <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/REVISIONS/matrix_generation/data/HTC/observed",pattern=m.w,full.names=T)
      # HiC.files <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/REVISIONS/matrix_generation/test/observed",pattern=m.w,full.names=T)
    }else{
      HiC.files <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/REVISIONS/matrix_generation/data/HTC/OE",pattern=m.w,full.names=T)
    }
    lapply(mes_limites[[m.w]],function(window){
      message(window)
      bin <- ((window*2)/sw)+1
      
      res.DSB <- lapply(unique(seqnames(DSB174)),function(chrom){
        message(chrom)
        f <- HiC.files[str_detect(HiC.files,str_c("_",chrom,"_"))]
        Replicate <- f %>% map(.%>% basename() %>% str_extract("CTRL|siSCC1"))
        
        
        diva.chr <- DSB174 %>% filter(seqnames ==chrom) %>%
          anchor_center() %>% mutate(width = mes_limites[[m.w]]*2)
        
        if(length(diva.chr)<2)
          return(NULL)
        
        res <- lapply(unique(Replicate),function(j){
          s <- f[grep(j,f)]
          if(length(s)<2){
            NULL
          }else{
            message(s)
            Cond <- s %>% map(.%>% basename() %>% str_extract("DIvA|OHT"))
            c(res1,val1) %<-% get_sub_matrix_Hic(file = s[[1]],bed = diva.chr,binsize = bin)
            c(res2,val2) %<-% get_sub_matrix_Hic(file = s[[2]],bed = diva.chr,binsize = bin)
            
            res1 <- lapply(res1,function(x){
              colnames(x) <- rownames(x) <- 1:bin
              x[lower.tri(x)] <- NA
              x %>% as.matrix() %>% melt() %>% drop_na()
            })  %>% bind_rows(.id="DSB")
            
            res2 <- lapply(res2,function(x){
              x <- x * (val1/val2)
              colnames(x) <- rownames(x) <- 1:bin
              x[lower.tri(x)] <- NA
              x %>% as.matrix() %>% melt() %>% drop_na()
            })  %>% bind_rows(.id="DSB")
            
            
            full <- list(res1,res2)
            names(full) <- Cond
            full %>% bind_rows(.id="Condition") %>% mutate(DSB = str_c(chrom,DSB,sep="_")) %>% mutate(Type = DSB.n) %>% mutate(mylim = window) %>% mutate(win = m.w)
          }
          
          
        })
        names(res) <- unique(Replicate)
        res%>% bind_rows(.id="Replicate")
        
        
        
      })
      
      res.DSB <- res.DSB %>% bind_rows()
      res.DSB %>% 
        spread(key = Condition,value = value) %>%
        mutate(ratio = OHT/DIvA) %>%
        mutate(log2ratio = log2(ratio)) %>% 
        mutate(log10ratio = log10(ratio)) %>% 
        drop_na()
      
      
    }) %>% bind_rows() 
  }) %>% bind_rows() 
})%>% bind_rows() 

# p1 <- res.boxplot %>% ggplot(aes(x=Type,y=ratio,fill=Type)) + geom_violin() + stat_summary(
#     geom = "point",
#     fun.y = "mean",
#     col = "black"
# ) + theme_classic()
bin <- ((1500000*2)/100000)+1
TADlim <- (bin-1)/3
TADs_middle <- 1:bin %>% median

res.boxplot <- res.boxplot %>% 
  mutate(TAD1= case_when(
    Var1 < (TADs_middle-(TADlim/2)) ~ "leftTAD",
    Var1 > (TADs_middle+(TADlim/2)) ~ "rightTAD",
    TRUE ~ "DamagedTAD",
    
  )) %>% 
  mutate(TAD2= case_when(
    Var2 < (TADs_middle-(TADlim/2)) ~ "leftTAD",
    Var2 > (TADs_middle+(TADlim/2)) ~ "rightTAD",
    TRUE ~ "DamagedTAD",
  )) %>% 
  unite("TAD",c("TAD1","TAD2")) %>% 
  filter(!TAD %in% c("leftTAD_rightTAD")) %>%  #On a pas toutes les possibilitées car j'ai tronqué la partie supérieur gauche de la matrice de comptage
  mutate(Comparison = case_when(
    TAD %in% c("leftTAD_leftTAD","rightTAD_rightTAD") ~ "AdjacentTAD_cis",
    TAD %in% c("DamagedTAD_DamagedTAD") ~ "DamagedTAD_cis",
    TAD %in% c("DamagedTAD_rightTAD","leftTAD_DamagedTAD") ~ "AdjacentTAD_DamagedTAD",
  ))



cc <- res.boxplot %>% group_by(Replicate,Var1,Var2,Type,mylim,win,Comparison) %>% summarise(DIvA=mean(DIvA,na.rm=T),OHT=mean(OHT,na.rm=T)) %>% mutate(ratio = OHT/DIvA) %>%
  mutate(log2ratio = log2(ratio))

p1 <- cc %>% filter(Comparison != "AdjacentTAD_cis",Type == "bless80",is.finite(log2ratio)) %>% 
  ggplot(aes(x=Type,y=log2ratio,fill=Replicate)) +facet_wrap(~Comparison,nrow=1) + geom_boxplot(width=0.5) + theme_classic() +theme(legend.position = "none",axis.text.x = element_text(angle = 90)) + coord_cartesian(ylim=c(-0.5,0.5))+ coord_cartesian(ylim=c(-0.3,0.3))
print(p1)

#Les pvals pour fig 4C papier coline
hihi1 <- cc %>% filter(Comparison != "AdjacentTAD_cis") %>% group_by(Comparison,Type) %>% nest()%>% 
  mutate(pval = map_dbl(data,function(x){
    wilcox.test(log2ratio~Replicate,data=x)$p.value
  })) %>% dplyr::select(-data) 
hihi1 %>% filter(Type =="bless80")
hihi2 <- cc %>% filter(Comparison != "AdjacentTAD_cis") %>% group_by(Type,Comparison,Replicate) %>% nest() %>% 
  mutate(pval = map_dbl(data,function(x){
    wilcox.test(x$DIvA,x$OHT,paired=T)$p.value
  })) %>% dplyr::select(-data)
