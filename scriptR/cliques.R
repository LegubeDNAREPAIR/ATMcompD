# /home/rochevin/Documents/PROJET_INGE/HiC_Coline/REVISIONS/ORIGINAL_SCRIPTS/igraph_cliques_trans_DSB.R
require(HiTC)
require(tidyverse)
require(plyranges)
require(reshape2)
require(keras)
require(rtracklayer)
require(cowplot)
require(igraph)
require(ggnetwork)
# require(igraph)
# require(ggnetwork)

Get1val <- function(my.wigs,one.w,x){
  lapply(split(x,droplevels(seqnames(x))),function(zz){
    message(unique(as.character(seqnames(zz))))
    cov <- one.w[[unique(as.character(seqnames(zz)))]]
    score <- Views( cov, start = start(zz), end = end(zz) ) %>% sum()
    tibble(wig = my.wigs,value = score,rowname = zz$name)
  }) %>% bind_rows()
}
seqlens <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 %>% seqlengths()
DSB174 <- read_bed("/mnt/NAS1/DATA/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")
bless80 <- read_bed("/mnt/NAS1/DATA/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
HR <-  read_bed("/mnt/NAS1/DATA/AsiSI/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed")
NHEJ <-  read_bed("/mnt/NAS1/DATA/AsiSI/BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed")
random80 <- read_bed("/mnt/NAS1/DATA/AsiSI/80random.bed")
random30 <- read_bed("/mnt/NAS1/DATA/AsiSI/30random.bed")
res <- read_tsv("../DSB_DSB_interaction_1mbEach.tsv.gz")

mes_windows <- c(500000)


# mes_treshold <- c(0.8,0.9,1,1.2,1.4)
mes_treshold <- c(0.9)
#REPLI SEQ
# FILTRER PAR ORI 
bw.files <- c(
  "repli-Seq_1_10_pOHTvsmOHT"="/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Repli-seq/BAMCOMPARE/repli-Seq_1_10_pOHTvsmOHT_frombam.bw",
  "repli-Seq_3_10_pOHTvsmOHT"="/mnt/NAS1/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Repli-seq/BAMCOMPARE/repli-Seq_3_10_pOHTvsmOHT_frombam.bw"
) %>% purrr::map(import.bw,as="RleList")



res.bw <- lapply(names(bw.files),function(one.n){
  message(one.n)
  lapply(mes_windows,function(window){
    diva.chr <- bless80 %>%
      anchor_center() %>% mutate(width = 500000)
    Get1val(one.n,bw.files[[one.n]],diva.chr) %>% mutate(w=window)
  }) %>% bind_rows()
}) %>% bind_rows()


my_list_of_special_bw <- c(
  "repli-Seq_1_10_pOHTvsmOHT",
  "repli-Seq_3_10_pOHTvsmOHT"
)

# mes_treshold <- c(0.8,0.9,1,1.2,1.4)
mes_treshold <- c(0.9)

for(treshold.FC in mes_treshold){
  message(treshold.FC)
  data.ratio.AB <- res %>%
    filter(bin1 %in% bless80$name,bin2 %in% bless80$name) %>%
    mutate(score = score * (val1/val2)) %>%
    dplyr::select(score,bin1,bin2,Condition,Manip) %>%
    group_by(bin1,bin2,Condition) %>%
    summarise(meanValue = mean(score)) %>%
    spread(key = Condition,value = meanValue) %>%
    mutate(OHT = OHT +1 ) %>%
    mutate(DIvA = DIvA + 1) %>%
    mutate(ratio = log2(OHT/DIvA)) %>% ungroup() %>%
    mutate(Diff = ifelse(abs(ratio) < treshold.FC,"NoDiff","Diff"))
  
  
  nodes <- data.ratio.AB %>% filter(Diff == "Diff") %>% dplyr::rename(from = bin1) %>% dplyr::rename(to = bin2) %>%
    mutate(NType = ifelse(ratio < 0,"Downregulated","Upregulated")) %>%
    filter(NType == "Upregulated")
  #cr <- spatstat::colourmap(colorRampPalette(c("yellow","red"))(1000), range=c(0,max(Infos$BLESS_pOHT)))
  Infos <- bless80 %>% as_tibble() %>%
    mutate(Type = case_when(
      name %in% HR$name ~ "HR",
      name %in% NHEJ$name ~ "NHEJ",
      TRUE ~ "None"
    ))  %>%
    dplyr::select(-seqnames:-strand,-score) %>%
    # left_join(res.LADS,by = c("name"="rowname")) %>%
    filter(name %in% unique(c(nodes$from,nodes$to)))
  # mutate(BlessColor = cr(BLESS_pOHT))
  
  
  g <- graph_from_data_frame(nodes, vertices=Infos,directed=F)
  
  my_cliques <- max_cliques(g,min=2)
  
  
  histoCliques <- my_cliques %>%
    map(function(clique){myGroupName = str_c(names(clique),collapse ="_");tibble(Group = myGroupName,name = names(clique),size = length(clique))}) %>%
    bind_rows() %>%
    # right_join(bless80 %>% as_tibble() %>% dplyr::select(name)) %>%
    # mutate(Group = ifelse(is.na(Group),name,Group)) %>%
    # mutate(size = ifelse(is.na(size),0,size)) %>%
    mutate(Type = case_when(
      name %in% HR$name ~ "HR",
      name %in% NHEJ$name ~ "NHEJ",
      TRUE ~ "None"
    ))
  
  for(my_wig in unique(res.bw$wig)){
    message(my_wig)
    sub.res.bw <- res.bw %>% filter(wig == my_wig)
    ccdicc <- ifelse(my_wig %in% my_list_of_special_bw,TRUE,FALSE)
    for(my_w in unique(sub.res.bw$w)){
      message(my_w)
      sub.sub.res.bw <- sub.res.bw %>% filter(w == my_w)
      
      p3 <-  histoCliques %>% dplyr::select(-Group) %>% group_by(size) %>% distinct() %>%
        right_join(bless80 %>% as_tibble() %>% dplyr::select(name)) %>%
        ungroup() %>%
        mutate(size = ifelse(is.na(size),0,size)) %>%
        left_join(sub.sub.res.bw,by = c("name"="rowname")) %>%
        group_by(size)%>%
        summarise(n = dplyr::n(),meanvalue=mean(value)) %>%
        ggplot(aes(x=as.factor(size),y=n))  + coord_flip() + scale_fill_viridis_c()
      if(ccdicc){
        p3 <- p3 + geom_histogram(aes(fill=meanvalue),stat="identity",position = "stack",col="black")
      }else{
        p3 <- p3 +  geom_histogram(aes(fill=log2(meanvalue)),stat="identity",position = "stack",col="black")
      }
      
      histoCounts <- histoCliques %>% dplyr::select(-Group) %>%
        group_by(size) %>% distinct() %>%
        left_join(sub.sub.res.bw,by = c("name"="rowname")) %>%
        group_by(size,Type)%>%
        summarise(n = dplyr::n(),meanvalue=mean(value)) %>%
        group_by(size) %>% mutate( ntot = sum(n)) %>% mutate(percent = (n/ntot)*100)
      
      p4 <- histoCounts%>%
        ggplot(aes(x=as.factor(size),y=percent,col=Type,label=Type))+
        scale_color_manual(values = c("#c0392b","#2980b9","#7f8c8d")) +
        scale_fill_gradient2(low = "#2980b9",high = "#f1c40f",midpoint = 0,mid = "black",limits = c(-1,1))
      
      
      
      
      histoNone <- histoCliques %>% dplyr::select(-Group) %>%
        filter(Type !="None") %>%
        group_by(size) %>% distinct() %>%
        left_join(sub.sub.res.bw,by = c("name"="rowname")) %>%
        group_by(size,Type)%>%
        summarise(n = dplyr::n(),meanvalue=mean(value)) %>%
        group_by(size) %>% mutate( ntot = sum(n)) %>% mutate(percent = (n/ntot)*100)
      
      p5 <- histoNone %>%
        ggplot(aes(x=as.factor(size),y=percent,col=Type,label=Type))+
        scale_color_manual(values = c("#c0392b","#2980b9","#7f8c8d"))
      
      histoMerge <- histoCliques %>% dplyr::select(-Group) %>%
        mutate(Type = ifelse(Type == "HR","HR","NHEJ+None")) %>%
        group_by(size) %>% distinct() %>%
        left_join(sub.sub.res.bw,by = c("name"="rowname")) %>%
        group_by(size,Type)%>%
        summarise(n = dplyr::n(),meanvalue=mean(value)) %>%
        group_by(size) %>% mutate( ntot = sum(n)) %>% mutate(percent = (n/ntot)*100)
      p6 <- histoMerge %>%
        ggplot(aes(x=as.factor(size),y=percent,col=Type,label=Type))+
        scale_color_manual(values = c("#c0392b","#2980b9","#7f8c8d"))
      
      
      
      
      if(ccdicc){
        p4 <- p4 + geom_histogram(aes(fill=meanvalue),stat="identity",position = "stack") +
          geom_label(position = position_stack(vjust = 0.5),fill="white") +coord_flip() +
          scale_fill_gradient2(low = "#2980b9",high = "#f1c40f",midpoint = 0,mid = "black")
        p5 <- p5 + geom_histogram(aes(fill=meanvalue),stat="identity",position = "stack") +
          geom_label(position = position_stack(vjust = 0.5),fill="white") +coord_flip()+
          scale_fill_gradient2(low = "#2980b9",high = "#f1c40f",midpoint = 0,mid = "black")
        p6 <- p6 + geom_histogram(aes(fill=meanvalue),stat="identity",position = "stack") +
          geom_label(position = position_stack(vjust = 0.5),fill="white") +coord_flip()+
          scale_fill_gradient2(low = "#2980b9",high = "#f1c40f",midpoint = 0,mid = "black")
        
      }else{
        p4 <- p4 +  geom_histogram(aes(fill=log2(meanvalue)),stat="identity",position = "stack") +
          geom_label(position = position_stack(vjust = 0.5),fill="white") +coord_flip() +
          scale_fill_viridis_c()
        p5 <- p5 +  geom_histogram(aes(fill=log2(meanvalue)),stat="identity",position = "stack") +
          geom_label(position = position_stack(vjust = 0.5),fill="white") +coord_flip() +
          scale_fill_viridis_c()
        p6 <- p6 +  geom_histogram(aes(fill=log2(meanvalue)),stat="identity",position = "stack") +
          geom_label(position = position_stack(vjust = 0.5),fill="white") +coord_flip() +
          scale_fill_viridis_c()
      }
      
      
      pdf(str_c("PLOTS/cliques/Distrib_by_clique_size_DSB_DSB_interaction_1mbEach",round(10^treshold.FC),my_wig,my_w,"All80bless_by.pdf",sep="_"),height=5,width=10)
      print(p3 + ggtitle(str_c("FC: ",round(10^treshold.FC)," wig: ",my_wig," window: ",my_w)))
      dev.off()
      pdf(str_c("PLOTS/cliques/NO_0_Distrib_HR_NHEJ_DSB_DSB_interaction_1mbEach",round(10^treshold.FC),my_wig,my_w,"All80bless_by.pdf",sep="_"),height=5,width=10)
      print(p4 + ggtitle(str_c("FC: ",round(10^treshold.FC)," wig: ",my_wig," window: ",my_w)))
      print(p5 + ggtitle(str_c("FC: ",round(10^treshold.FC)," wig: ",my_wig," window: ",my_w)))
      print(p6 + ggtitle(str_c("FC: ",round(10^treshold.FC)," wig: ",my_wig," window: ",my_w)))
      dev.off()
    }
    
  }
  
  
  
}