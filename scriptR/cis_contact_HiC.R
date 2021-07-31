#FROM /home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/scripts/HeatMap_and_boxplot_HiC.R
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

get_trans_matrix_Hic <- function(file,bed,binsize){
  c(my.ranges,intdata) %<-% file
  
  val1 <- intdata %>% sum
  
  overlaps <- as(findOverlaps(bed,my.ranges), "List")
  overlaps.bin <- lapply(overlaps,length) %>% unlist()
  overlaps <- overlaps[overlaps.bin == binsize]
  
  res <- lapply(1:length(overlaps),function(i){
    x <- overlaps[[i]]
    if(as.character(strand(bed[i])) == "-"){
      x <- rev(x)
    }
    intdata[x,x] %>% round()
  })
  return(list(res,val1))
}

get_nb_read_by_chr <- function(file){
  HTC <- loadRData(file)
  HTC %>% intdata %>% sum
}

plot_ratio <-function(data.ratio,my.quantile = 0.95,facet=TRUE,window = 1000000,fixed = F,fixed.limite = c(0.1,0.1),DSB = "DSB"){
  if(fixed ==T){
    limite <- fixed.limite
    if(length(limite)<2){
      limite <- rep(limite,2)
    }
  }else{
    limite <- data.ratio %>% pull(ratio) %>% quantile(my.quantile) %>% as.numeric()
    limite <- rep(limite,2)
  }
  
  debut <- data.ratio$Var1 %>% levels %>% as.numeric() %>% min()
  milieu <- data.ratio$Var1 %>% levels %>% as.numeric() %>% median()
  fin <- data.ratio$Var1 %>% levels %>% as.numeric() %>% max()
  p.ratio <- data.ratio %>%
    mutate(ratio = ifelse(ratio > limite[2],limite[2],ratio)) %>%
    mutate(ratio = ifelse(ratio < -limite[1],-limite[1],ratio)) %>%
    ggplot(aes(x=Var1,y=Var2,fill=ratio)) + geom_tile() +
    scale_fill_gradient2(low = "#2980b9",high = "#f1c40f",midpoint = 0,mid = "black",limits = c(-limite[1],limite[2])) +
    scale_x_discrete(name = 'Position',
                     breaks = c(debut,milieu,fin),
                     labels = c(-window,DSB,window)
    ) +
    scale_y_discrete(name = 'Position',
                     breaks = c(debut,milieu,fin),
                     labels = c(-window,DSB,window)
    ) +
    geom_vline(xintercept = milieu,linetype ="dashed",col="red")+
    theme_classic() +
    theme(legend.position="bottom")
  if(facet == TRUE){
    p.ratio + facet_wrap(~Replicate,ncol=1)
  }else{
    p.ratio
  }
}


plot_diff <-function(data.diff,my.quantile = 0.90,facet=TRUE,window = 1000000,fixed = F,fixed.limite = 10,DSB = "DSB"){
  if(fixed ==T){
    limite <- fixed.limite
  }else{
    limite <- data.diff %>% pull(diff) %>% quantile(my.quantile) %>% as.numeric()
  }
  
  debut <- data.diff$Var1 %>% levels %>% as.numeric() %>% min()
  milieu <- data.diff$Var1 %>% levels %>% as.numeric() %>% median()
  fin <- data.diff$Var1 %>% levels %>% as.numeric() %>% max()
  p.diff <- data.diff %>%
    mutate(diff = ifelse(diff > limite,limite,diff)) %>%
    mutate(diff = ifelse(diff < -limite,-limite,diff)) %>%
    ggplot(aes(x=Var1,y=Var2,fill=diff)) + geom_tile() +
    scale_fill_gradient2(low = "#2980b9",high = "#f1c40f",midpoint = 0,mid = "black",limits = c(-limite,limite)) +
    scale_x_discrete(name = 'Position',
                     breaks = c(debut,milieu,fin),
                     labels = c(-window,DSB,window)
    ) +
    scale_y_discrete(name = 'Position',
                     breaks = c(debut,milieu,fin),
                     labels = c(-window,DSB,window)
    ) +
    geom_vline(xintercept = milieu,linetype ="dashed",col="red")+
    theme_classic() +
    theme(legend.position="bottom")
  
  if(facet == TRUE){
    p.diff + facet_wrap(~Replicate,ncol=1)
  }else{
    p.diff
  }
}
#ADD binned 5kb genomic region with no damage and lot of damage (excluding damages)
##USING DSB
bless80 = read_bed("/mnt/NAS/DATA/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
HR =  read_bed("/mnt/NAS/DATA/AsiSI/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed")
NHEJ =  read_bed("/mnt/NAS/DATA/AsiSI/BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed")
Random80 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/80random.bed")
uncut <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/ASIsites_hg19_without_bless80.bed")


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
    
    
    HiC.files <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/HTC",pattern=m.w,full.names=T)
    # HiC.files <- list.files("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/HTC/observed",pattern=str_c("_",m.w),full.names=T) %>% str_subset("_OHT_|_DIvA_")
    lapply(mes_limites[[m.w]],function(window){
      message(window)
      bin <- ((window*2)/sw)+1
      message(bin)
      
      res.DSB <- lapply(unique(seqnames(DSB174)),function(chrom){
        message(chrom)
        f <- HiC.files[str_detect(HiC.files,str_c("_",chrom,"_"))]
        Replicate <- f %>% map(.%>% basename() %>% str_extract("manip[A-Z]+|rep[0-9]+|HiC_D"))
        
        
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



cc <- res.boxplot %>% group_by(Var1,Var2,Type,mylim,win,Comparison)  %>% summarise(DIvA=mean(DIvA,na.rm=T),OHT=mean(OHT,na.rm=T)) %>% mutate(ratio = OHT/DIvA) %>%
  mutate(log2ratio = log2(ratio))

pdf("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/boxplot_Hich_average_ratio_meanReplicate_DSBvsAdj_log2ratio.pdf",width=6,height=6)
for(my.type in unique(cc$Type)){
  p1 <- cc %>% filter(Comparison != "AdjacentTAD_cis",Type == my.type) %>% 
    ggplot(aes(x=Comparison,y=log2ratio,fill=Comparison)) +facet_wrap(~Type,nrow=1) + geom_boxplot(width=0.5) + theme_classic() +theme(legend.position = "none",axis.text.x = element_text(angle = 90)) + coord_cartesian(ylim=c(-0.5,0.5))
  print(p1)
}

dev.off()

pdf("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/boxplot_Hich_average_ratio_meanReplicate_TADComparisonAdj_log2ratio.pdf",width=6,height=6)
for(my.type in unique(cc$Type)){
  p1 <- cc %>% filter(Comparison != "AdjacentTAD_DamagedTAD",Type == my.type) %>% 
    ggplot(aes(x=Comparison,y=log2ratio,fill=Comparison)) +facet_wrap(~Type,nrow=1)  + geom_boxplot(width=0.5) + theme_classic() +theme(legend.position = "none",axis.text.x = element_text(angle = 90)) + coord_cartesian(ylim=c(-0.3,0.3))
  print(p1)
}
dev.off()

#Les pvals pour fig 4C papier coline
hihi1 <- cc %>% filter(Comparison != "AdjacentTAD_cis") %>% group_by(Type) %>% nest()%>% 
  mutate(pval = map_dbl(data,function(x){
    wilcox.test(log2ratio~Comparison,data=x)$p.value
  })) %>% dplyr::select(-data)
hihi15 <- cc %>% filter(Comparison != "AdjacentTAD_DamagedTAD",Type %in% c("bless80","Random80")) %>% group_by(Type) %>% nest()%>% 
  mutate(pval = map_dbl(data,function(x){
    wilcox.test(log2ratio~Comparison,data=x)$p.value
  })) %>% dplyr::select(-data)


hihi2 <- cc %>% filter(Comparison != "AdjacentTAD_cis") %>% group_by(Type,Comparison) %>% nest() %>% 
  mutate(pval = map_dbl(data,function(x){
    wilcox.test(x$DIvA,x$OHT,paired=T)$p.value
  })) %>% dplyr::select(-data)

hihi3 <- cc %>% filter(Comparison != "AdjacentTAD_cis",Type %in% c("bless80","uncut")) %>% group_by(Comparison) %>% nest()%>% 
  mutate(pval = map_dbl(data,function(x){
    wilcox.test(log2ratio~Type,data=x)$p.value
  })) %>% dplyr::select(-data)

hihi4 <- cc  %>% group_by(Type,Comparison)  %>% nest() %>% 
  mutate(pval = map_dbl(data,function(x){
    wilcox.test(x$log2ratio)$p.value
  })) %>% dplyr::select(-data)

#SAME BUT WITH D EXPERIMENTS

ctrlCOND <- "DIvA"
res.boxplot <- lapply(names(MES_DSB),function(DSB.n){
  message(DSB.n)
  DSB174 <- MES_DSB[[DSB.n]]
  
  
  
  lapply(names(mes_windows),function(m.w){
    message(m.w)
    sw <- mes_windows[[m.w]]
    
    
    HiC.files <- list.files("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/HTC/observed",pattern=str_c("_",m.w),full.names=T)
    lapply(mes_limites[[m.w]],function(window){
      message(window)
      bin <- ((window*2)/sw)+1
      
      res.DSB <- lapply(unique(seqnames(DSB174)),function(chrom){
        message(chrom)
        f <- HiC.files[str_detect(HiC.files,str_c("_",chrom,"_"))]
        Replicate <- f %>% map(.%>% basename() %>% str_remove_all("HTC_observed_KR_|_manip[AB]|HTC_observed_KR_HiC_D_|_chr[0-9A-Z]+_[0-9]+kb.RData")) %>% unlist()
        names(f) <- Replicate
        # Replicate <- f %>% map(.%>% basename() %>% str_extract("manip[A-Z]+|rep[0-9]+"))
        
        diva.chr <- DSB174 %>% filter(seqnames ==chrom) %>%
          anchor_center() %>% mutate(width = mes_limites[[m.w]]*2)
        
        if(length(diva.chr)<2)
          return(NULL)
        
        
        f_ctrl <- Replicate %>% str_subset(glue::glue("{ctrlCOND}$"))
        f_exp <- Replicate %>% str_subset(glue::glue("{ctrlCOND}$"),negate = T)
        
        
        c(res_CTRL,val_CTRL) %<-% get_sub_matrix_Hic(file = f[f_ctrl],bed = diva.chr,binsize = bin)
        res_CTRL <- lapply(res_CTRL,function(x){
          colnames(x) <- rownames(x) <- 1:bin
          x[lower.tri(x)] <- NA
          x %>% as.matrix() %>% melt() %>% drop_na()
        })  %>% bind_rows(.id="DSB")
        res <- lapply(f[f_exp],function(s){
          c(sres_exp,sval_exp) %<-% get_sub_matrix_Hic(file = s,bed = diva.chr,binsize = bin)
          
          sres_exp <- lapply(sres_exp,function(x){
            x <- x * (val_CTRL/sval_exp)
            colnames(x) <- rownames(x) <- 1:bin
            x[lower.tri(x)] <- NA
            x %>% as.matrix() %>% melt() %>% drop_na()
          })  %>% bind_rows(.id="DSB")
          
          sres_exp
        })
        
        
        
        
        
        res[["CTRL"]] <- res_CTRL
        
        
        
        res%>% bind_rows(.id="Condition") %>% mutate(DSB = str_c(chrom,DSB,sep="_")) %>% mutate(Type = DSB.n) %>% mutate(mylim = window) %>% mutate(win = m.w)
      })
      res.DSB <- res.DSB %>% bind_rows()
      # res.DSB %>%
      #   spread(key = Condition,value = value) %>%
      #   mutate(ratio = HiC_D_OHT/CTRL) %>%
      #   mutate(log2ratio = log2(ratio)) %>%
      #   mutate(log10ratio = log10(ratio)) %>%
      #   drop_na()
      res.DSB %>%
        spread(key = Condition,value = value) %>%
        gather(key = "Condition",value = value,-DSB:-CTRL) %>%
        # mutate(OHT = OHT * (487/393)) %>%
        mutate(ratio = (value/(CTRL))) %>%
        mutate(log2ratio = log2(ratio)) %>%
        mutate(log10ratio = log10(ratio)) %>%
        mutate(Condition = glue::glue("{Condition}/{ctrlCOND}")) %>%
        drop_na()
      
      
    }) %>% bind_rows() 
  }) %>% bind_rows() 
})%>% bind_rows() 
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



cc <- res.boxplot %>% group_by(Var1,Var2,Type,mylim,win,Comparison,Condition)  %>% summarise(CTRL=mean(CTRL,na.rm=T),value=mean(value,na.rm=T)) %>% mutate(ratio = value/CTRL) %>%
  mutate(log2ratio = log2(ratio))

pdf("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/boxplot_Hich_average_ratio_meanReplicate_DSBvsAdj_log2ratio_manipD.pdf",width=12,height=6)
for(my.type in unique(cc$Type)){
  p1 <- cc %>% filter(Comparison != "AdjacentTAD_cis",Type == my.type) %>% 
    ggplot(aes(x=Comparison,y=log2ratio,fill=Comparison)) +facet_wrap(Condition~Type,nrow=1) + geom_boxplot(width=0.5) + theme_classic() +theme(legend.position = "none",axis.text.x = element_text(angle = 90)) + coord_cartesian(ylim=c(-0.5,0.5))
  print(p1)
}
# 
# hu <- lapply(unique(cc$Type),function(my.type){
#   cc %>% filter(Comparison != "AdjacentTAD_cis",Type == my.type) %>%
#     ggplot(aes(x=Comparison,y=log2ratio,fill=Comparison)) +facet_wrap(Condition~Type,nrow=1) + geom_boxplot(width=0.5) + theme_classic() +theme(legend.position = "none",axis.text.x = element_text(angle = 90)) + coord_cartesian(ylim=c(-0.5,0.5))
# })
dev.off()


#Les pvals pour fig 4C papier coline
hihi1 <- cc %>% filter(Comparison != "AdjacentTAD_cis") %>% group_by(Type,Condition) %>% nest()%>% 
  mutate(pval = map_dbl(data,function(x){
    wilcox.test(log2ratio~Comparison,data=x)$p.value
  })) %>% dplyr::select(-data)

hihi4 <- cc  %>% group_by(Type,Condition,Comparison)  %>% nest() %>% 
  mutate(pval = map_dbl(data,function(x){
    wilcox.test(x$log2ratio)$p.value
  })) %>% dplyr::select(-data)

hihi2 <- cc %>% filter(Comparison != "AdjacentTAD_cis") %>% group_by(Type,Condition,Comparison) %>% nest() %>% 
  mutate(pval = map_dbl(data,function(x){
    wilcox.test(x$CTRL,x$value,paired=T)$p.value
  })) %>% dplyr::select(-data) %>% filter(Type =="bless80")


hihi3 <- cc %>% filter(Comparison == "DamagedTAD_cis") 
hihi3 %>% filter(Condition %in% c("HiC_D_OHT/DIvA","HiC_D_OHTATMi/DIvA")) %>% wilcox.test(log2ratio ~ Condition,data=.) %>% .$p.value 
hihi3 %>% filter(Condition %in% c("HiC_D_OHT/DIvA","HiC_D_OHTDNAPKi/DIvA")) %>% wilcox.test(log2ratio ~ Condition,data=.) %>% .$p.value%>% format.pval(digits = 2)
