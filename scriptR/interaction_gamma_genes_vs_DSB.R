#From /home/rochevin/Documents/PROJET_INGE/HiC_Coline/REVISIONS/ORIGINAL_SCRIPTS/Heatmap_HiC_data_DSB.R
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


#SITE133,SITE1160,SITE1159,SITE379,SITE495,SITE533,SITE628,SITE713,SITE716
get_sub_matrix_Hic_DSB_genes <- function(file,bed1,bed2){
  HTC <- loadRData(file)
  rows.chr <- x_intervals(HTC)
  val1 <- HTC %>% intdata %>% sum
  # HTC <- normPerReads(HTC)
  
  overlaps1 <- findOverlaps(bed1,rows.chr)
  overlaps2 <- findOverlaps(bed2,rows.chr)
  my_test <- !subjectHits(overlaps1) == subjectHits(overlaps2)
  if(FALSE %in% my_test)
    message(bed1$name)
  overlaps2 <- overlaps2[my_test]
  overlaps2 <- as(overlaps2, "List")
  overlaps1 <- as(overlaps1, "List")
  res <- lapply(1:length(overlaps1),function(i){
    x <- overlaps1[[i]]
    lapply(1:length(overlaps2),function(j){
      y <- overlaps2[[j]]
      intdata(HTC)[x,y] %>% enframe() %>% dplyr::select(-name) %>% mutate(name1 = bed1[i]$name)%>% mutate(name2 = bed2[j]$gene_id) %>% mutate(bin_pos = as.character(rows.chr[y]))
    }) %>% bind_rows()
    
  }) %>% bind_rows()
  return(list(res,val1))
}

get_nb_read_by_chr <- function(file){
  HTC <- loadRData(file)
  HTC %>% intdata %>% sum
}
#Prepare dataset
DE_DIVA <- PhDfunc::GetDE_DIvA()
bless80 <- read_bed("/mnt/NAS/DATA/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed") 
domains <- bless80 %>% anchor_center() %>% mutate(width=2000000)
ens.genes <- read.table("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed",sep="\t",h=T) %>% GRanges()
ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))
ens.genes$name <- ens.genes$gene_id
ens.extended <- ens.genes  %>% anchor_start() %>% stretch(1000) %>% anchor_end() %>% stretch(1000)
direct.damaged <- ens.extended %>% filter_by_overlaps(bless80)
genes.domains <- ens.genes %>% filter(gene_id %in% DE_DIVA$rowname) %>% filter(!gene_id %in% direct.damaged$gene_id) %>% 
  filter_by_overlaps(domains)
TSS.domains <- genes.domains %>% anchor_5p() %>% mutate(width=1) %>% anchor_center() %>% mutate(width=1000)

m.w <- "10kb"
sw <- 10000
HiC.files <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/HTC/OE",pattern=m.w,full.names=T)
res.HiC <- mclapply(split(bless80,bless80$name),function(oneDSB){
  
  chrom <- unique(seqnames(oneDSB))
  message(chrom)
  f <- HiC.files[str_detect(HiC.files,str_c("_",chrom,"_"))]
  Replicate <- f %>% map(.%>% basename() %>% str_extract("manip[A-Z]+|rep[0-9]+"))
  
  my_genes <- genes.domains %>% filter_by_overlaps(domains[domains$name == oneDSB$name])
  if(length(my_genes) < 1)
    return(NULL)
  res <- lapply(unique(Replicate),function(j){
    s <- f[grep(j,f)]
    if(length(s)<2){
      NULL
    }else{
      message(s)
      Cond <- s %>% map(.%>% basename() %>% str_extract("DIvA|OHT"))
      c(res1,val1) %<-% get_sub_matrix_Hic_DSB_genes(file = s[[1]],bed1 = oneDSB,bed2 = my_genes)
      c(res2,val2) %<-% get_sub_matrix_Hic_DSB_genes(file = s[[2]],bed1 = oneDSB,bed2 = my_genes)
      
      res2 <- res2 %>% mutate(value = value * (val1/val2))
      
      

      full <- list(res1,res2)
      names(full) <- Cond
      full %>% bind_rows(.id="Condition")
    }
    
    
  })
  names(res) <- unique(Replicate)
  res%>% bind_rows(.id="Replicate")
},mc.cores=6) %>% bind_rows()

res.HiC <- res.HiC %>% group_by(name1,name2,Replicate,Condition) %>%
  summarise(value = mean(value)) %>%
  ungroup() %>% group_by(name1,name2,Condition) %>%
  summarise(meanVal = mean(value))%>% 
  spread(key = Condition,value=meanVal) %>% 
  filter(DIvA != 0) %>% 
  filter(OHT != 0) %>% 
  mutate(ratio = OHT/DIvA)

#Get category of gene interactions
distinct_genes <- res.HiC %>% distinct(name2)
res.HiC <- res.HiC %>% filter(name2 %in% distinct_genes$name2)
res.HiC <- res.HiC %>% ungroup() %>% 
  mutate(Cat_Interaction_20_DIVA = ntile(DIvA,20)) %>% 
  mutate(Cat_Interaction_100_DIVA = ntile(DIvA,100)) %>% 
  mutate(CatVal_DIVA = case_when(
    DIvA <= 0.3 ~ "Bottom",
    DIvA >= 3  ~ "Top"
  )) %>% 
  mutate(Cat_Interaction_20_OHT = ntile(OHT,20)) %>% 
  mutate(Cat_Interaction_100_OHT = ntile(OHT,100)) %>% 
  mutate(CatVal_OHT = case_when(
    OHT <= 0.3 ~ "Bottom",
    OHT >= 3  ~ "Top"
  )) %>% 
  mutate(Cat_Interaction_20_ratio = ntile(ratio,20)) %>% 
  mutate(Cat_Interaction_100_ratio = ntile(ratio,100)) %>% 
  mutate(CatVal_ratio = case_when(
    ratio <= 0.5 ~ "Bottom",
    ratio >= 2  ~ "Top"
  )) 


#On remove les genes qui pourraient etre dans plusieurs categories car overlap entre domaines

wigs_for_HiC <- c(
  # PhDfunc::GetBWList()[c("Tot_Pol2_mOHT","Tot_Pol2_pOHT","Ser2P_Pol2_mOHT","Ser2P_Pol2_pOHT","T1P_Pol2_mOHT","T1P_Pol2_pOHT")],
  "BruSeq_30vsU"="/mnt/NAS/DATA/DATA_FROM_OTHER_LABS_and_DATABASES/Bru-Seq/Bru-Seq_Ianelli/BAMCOMPARE/BruSeq_30vsU_frombw.bw",
  "BruSeq_60vsU"="/mnt/NAS/DATA/DATA_FROM_OTHER_LABS_and_DATABASES/Bru-Seq/Bru-Seq_Ianelli/BAMCOMPARE/BruSeq_60vsU_frombw.bw",
  "BruSeq_150vsU"="/mnt/NAS/DATA/DATA_FROM_OTHER_LABS_and_DATABASES/Bru-Seq/Bru-Seq_Ianelli/BAMCOMPARE/BruSeq_150vsU_frombw.bw",
  "BruSeq_240vsU"="/mnt/NAS/DATA/DATA_FROM_OTHER_LABS_and_DATABASES/Bru-Seq/Bru-Seq_Ianelli/BAMCOMPARE/BruSeq_240vsU_frombw.bw",
  # "Ser5P_Pol2_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Ser5P-Pol2/PROCESSED/BY_US/WIGGLE/H37HYBGX3_Ser5P_Pol2_DIVA_17s001626-1-1_Clouaire_lane117s001626_sequence_normalized.bw",
  # "Ser5P_Pol2_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Ser5P-Pol2/PROCESSED/BY_US/WIGGLE/H37HYBGX3_Ser5P_Pol2_OHT_17s001627-1-1_Clouaire_lane117s001627_sequence_normalized.bw",
  # "Ser7P_Pol2_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Ser7P-Pol2/PROCESSED/BY_US/WIGGLE/H37HYBGX3_Ser7P_Pol2_DIVA_17s001628-1-1_Clouaire_lane117s001628_sequence_normalized.bw",
  # "Ser7P_Pol2_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Pol2/EMBL_Ser7P-Pol2/PROCESSED/BY_US/WIGGLE/H37HYBGX3_Ser7P_Pol2_OHT_17s001629-1-1_Clouaire_lane117s001629_sequence_normalized.bw",
  "RNA_C1_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYL5LBGX2/HYL5LBGX2_RNA_C1_1_17s002464-1-1_Clouaire_lane117s002464_hg19plusERCC_normalized.bw",
  "RNA_C1_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYL5LBGX2/HYL5LBGX2_RNA_C1_OHT_1_17s002465-1-1_Clouaire_lane117s002465_hg19plusERCC_normalized.bw",
  "RNA_C2_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYMHNBGX2/HYMHNBGX2_RNA_C1_2_17s002468-1-1_Clouaire_lane117s002468_hg19plusERCC_normalized.bw",
  "RNA_C2_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYMHNBGX2/HYMHNBGX2_RNA_C1_OHT_2_17s002469-1-1_Clouaire_lane117s002469_hg19plusERCC_normalized.bw",
  "RNA_C1_pvsmOHT"="/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/BAMCOMPARE/HYL5LBGX2_RNA_C1_bamcompare_logFC_normalized.bw",
  "RNA_C2_pvsmOHT"="/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/BAMCOMPARE/HYMHNBGX2_RNA_C2_bamcompare_logFC_normalized.bw",
  "Bru-seq_Induced_30min"="/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/DDR_Ianelli/PROCESSED_08102018/mapping/BIGWIG/SRR5441994_normalized.bw",
  "Bru-seq_Uninduced"="/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/DDR_Ianelli/PROCESSED_08102018/mapping/BIGWIG/SRR5441998_normalized.bw"
) %>% mclapply(import.bw,as="RleList",mc.cores=4)
wigs_histones <- c(
  "h3k79me2_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K79me2/BGI_RUN1_201203/PROCESSED/BIGWIG/h3k79me2_mOHT_normalized_hg19.bw",
  "h3k79me2_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K79me2/BGI_RUN1_201203/PROCESSED/BIGWIG/h3k79me2_pOHT_normalized_hg19.bw",
  "H3K4me3_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K4me3/BGI_RUN3_201412/PROCESSED/BIGWIG/H3K4me3_normalized_hg19.bw",
  "H3K4me3_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K4me3/BGI_RUN3_201412/PROCESSED/BIGWIG/H3K4me3_OHT_normalized_hg19.bw",
  "h2bub_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H2BK120Ub/BGI_RUN1_201203/PROCESSED/BIGWIG/h2bub_mOHT_normalized_hg19.bw",
  "h2bub_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H2BK120Ub/BGI_RUN1_201203/PROCESSED/BIGWIG/h2bub_pOHT_normalized_hg19.bw"
)%>% map(import.bw,as="RleList")
wigs_histones_bamcompare <- c(
  "h3k79me2_pvsmOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K79me2/BGI_RUN1_201203/PROCESSED/BIGWIG/h3k79me2_BGI_RUN_1_bin50_log2ratio.bw",
  "H3K4me3_pvsmOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K4me3/BGI_RUN3_201412/PROCESSED/BIGWIG/H3K4me3_BGI_RUN_3_bin50_log2ratio.bw",
  "h2bub_pvsmOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H2BK120Ub/BGI_RUN1_201203/PROCESSED/BIGWIG/h2bub_BGI_RUN_1_bin50_log2ratio.bw"
)%>% map(import.bw,as="RleList")

wigs_histones_bwcompare <- c(
  "h3k79me2_pvsmOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K79me2/BGI_RUN1_201203/PROCESSED/BIGWIG/bigwigCompare_h3k79me2_BGI_RUN_1_pvsmOHT.bw",
  "H3K4me3_pvsmOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K4me3/BGI_RUN3_201412/PROCESSED/BIGWIG/bigwigCompare_H3K4me3_BGI_RUN_3_pvsmOHT.bw",
  "h2bub_pvsmOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H2BK120Ub/BGI_RUN1_201203/PROCESSED/BIGWIG/bigwigCompare_h2bub_BGI_RUN_1_pOHTvsmOHT.bw"
)%>% map(import.bw,as="RleList")

#H2BK120ub / H3K4me3 / H3K79me2 -> sur -1TSS/+1 TSS
#RNA seq m et p OHT, Bruseq Untreated et 30m 
#Tester la quantité sur les bin directement
#Ajouter boxplot de l'interaction
p1 <- res.HiC %>% gather(key=Cat,value=value,-name1:-OHT) %>% filter(Cat == "Cat_Interaction_100")%>% ggplot(aes(x=as.factor(value),y=DIvA)) + geom_boxplot() + facet_wrap(~Cat,ncol=1,scales="free_y") + theme_classic()
p2 <- res.HiC %>% gather(key=Cat,value=value,-name1:-OHT) %>% filter(Cat == "Cat_Interaction_20")%>% ggplot(aes(x=as.factor(value),y=DIvA)) + geom_boxplot() + facet_wrap(~Cat,ncol=1,scales="free_y") + theme_classic()

pdf("../results/Interaction_gene_gamma_DSB_by_cat.pdf")
print(p1)
print(p2)
dev.off()
Get1valMean <- function (Name, one.w, x) 
{
  require(magrittr)
  lapply(split(x, droplevels(seqnames(x))), function(zz) {
    message(unique(as.character(seqnames(zz))))
    cov <- one.w[[unique(as.character(seqnames(zz)))]]
    score <- IRanges::Views(cov, start = start(zz), end = end(zz)) %>% 
      mean()
    tibble::tibble(wig = Name, value = score, rowname = zz$name)
  }) %>% bind_rows()
}
res.boxplot.HiC <- lapply(wigs_for_HiC,function(mon_wig){
  Get1valMean(Name = "",x = genes.domains,one.w = mon_wig) %>% dplyr::select(-wig)
}) %>% bind_rows(.id="wig")%>% mutate(Condition = str_extract(wig,"mOHT|pOHT|Induced_30min|Uninduced|pvsmOHT")) %>% mutate(wig = str_remove(wig,"_mOHT|_pOHT|_Induced_30min|_Uninduced|_pvsmOHT"))

res.boxplot.HiC.histone <- lapply(wigs_histones,function(mon_wig){
  Get1valMean(Name = "",x = TSS.domains,one.w = mon_wig) %>% dplyr::select(-wig)
}) %>% bind_rows(.id="wig")%>% mutate(Condition = str_extract(wig,"mOHT|pOHT|Induced_30min|Uninduced|pvsmOHT")) %>% mutate(wig = str_remove(wig,"_mOHT|_pOHT|_Induced_30min|_Uninduced|_pvsmOHT"))
res.boxplot.HiC.histone.bc <- lapply(wigs_histones_bamcompare,function(mon_wig){
  Get1valMean(Name = "",x = TSS.domains,one.w = mon_wig) %>% dplyr::select(-wig)
}) %>% bind_rows(.id="wig")%>% mutate(Condition = str_extract(wig,"mOHT|pOHT|Induced_30min|Uninduced|pvsmOHT")) %>% mutate(wig = str_remove(wig,"_mOHT|_pOHT|_Induced_30min|_Uninduced|_pvsmOHT"))

res.boxplot.HiC.histone.bwc <- lapply(wigs_histones_bwcompare,function(mon_wig){
  Get1valMean(Name = "",x = TSS.domains,one.w = mon_wig) %>% dplyr::select(-wig)
}) %>% bind_rows(.id="wig")%>% mutate(Condition = str_extract(wig,"mOHT|pOHT|Induced_30min|Uninduced|pvsmOHT")) %>% mutate(wig = str_remove(wig,"_mOHT|_pOHT|_Induced_30min|_Uninduced|_pvsmOHT"))


ccplot <- res.boxplot.HiC %>% left_join(dplyr::select(res.HiC,-name1,-DIvA:-ratio),by=c("rowname"="name2"))

p1 <- cowplot::plot_grid(plotlist = list(
  ccplot %>% filter(Cat_Interaction_20_DIVA %in% c(1,20))%>% filter(wig %in% c("Bru-seq")) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_DIVA),y=value,fill=Condition)) + geom_boxplot()  +scale_y_log10()+ facet_wrap(~wig,ncol=2,scales="free") + theme_classic() + theme(legend.position = "bottom"),
  ccplot %>% filter(Cat_Interaction_100_DIVA %in% c(1,100)) %>% filter(wig %in% c("Bru-seq")) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_DIVA),y=value,fill=Condition)) + geom_boxplot()  +scale_y_log10()+ facet_wrap(~wig,ncol=2,scales="free") + theme_classic()+ theme(legend.position = "bottom"),
  ccplot %>% filter(!is.na(CatVal_DIVA)) %>% filter(wig %in% c("Bru-seq")) %>% ggplot(aes(x=as.factor(CatVal_DIVA),y=value,fill=Condition)) + geom_boxplot()  +scale_y_log10()+ facet_wrap(~wig,ncol=2,scales="free") + theme_classic()+ theme(legend.position = "bottom"),
  ccplot %>% filter(Cat_Interaction_20_OHT %in% c(1,20))%>% filter(wig %in% c("Bru-seq")) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_OHT),y=value,fill=Condition)) + geom_boxplot()  +scale_y_log10()+ facet_wrap(~wig,ncol=2,scales="free") + theme_classic()+ theme(legend.position = "bottom"),
  ccplot %>% filter(Cat_Interaction_100_OHT %in% c(1,100)) %>% filter(wig %in% c("Bru-seq")) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_OHT),y=value,fill=Condition)) + geom_boxplot()  +scale_y_log10()+ facet_wrap(~wig,ncol=2,scales="free") + theme_classic()+ theme(legend.position = "bottom"),
  ccplot %>% filter(!is.na(CatVal_OHT)) %>% filter(wig %in% c("Bru-seq")) %>% ggplot(aes(x=as.factor(CatVal_OHT),y=value,fill=Condition)) + geom_boxplot()  +scale_y_log10()+ facet_wrap(~wig,ncol=2,scales="free") + theme_classic()+ theme(legend.position = "bottom")
),ncol=3)

pbruseq_OHT <- cowplot::plot_grid(plotlist = list(
  ccplot %>% filter(Cat_Interaction_20_OHT %in% c(1,20))%>% filter(wig %in% c("Bru-seq")) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_OHT),y=value,fill=Condition)) + geom_boxplot()  +scale_y_log10()+ facet_wrap(~wig,ncol=2,scales="free") + theme_classic()+ theme(legend.position = "bottom"),
  ccplot %>% filter(Cat_Interaction_100_OHT %in% c(1,100)) %>% filter(wig %in% c("Bru-seq")) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_OHT),y=value,fill=Condition)) + geom_boxplot()  +scale_y_log10()+ facet_wrap(~wig,ncol=2,scales="free") + theme_classic()+ theme(legend.position = "bottom")
),ncol=2)

ccplot %>% filter(Cat_Interaction_100_DIVA %in% c(1,100)) %>% filter(wig %in% c("Bru-seq")) %>% group_by(Cat_Interaction_100_DIVA) %>% nest() %>% pull(data) %>% .[[1]] %>% wilcox.test(value ~ Condition,data=.) %>% .$p.value
ccplot %>% filter(Cat_Interaction_100_DIVA %in% c(1,100)) %>% filter(wig %in% c("Bru-seq")) %>% group_by(Cat_Interaction_100_DIVA) %>% nest() %>% pull(data) %>% .[[2]] %>% wilcox.test(value ~ Condition,data=.) %>% .$p.value
ccplot %>% filter(Cat_Interaction_100_OHT %in% c(1,100)) %>% filter(wig %in% c("Bru-seq")) %>% group_by(Cat_Interaction_100_OHT) %>% nest() %>% pull(data) %>% .[[1]] %>% wilcox.test(value ~ Condition,data=.) %>% .$p.value
ccplot %>% filter(Cat_Interaction_100_OHT %in% c(1,100)) %>% filter(wig %in% c("Bru-seq")) %>% group_by(Cat_Interaction_100_OHT) %>% nest() %>% pull(data) %>% .[[2]] %>% wilcox.test(value ~ Condition,data=.) %>% .$p.value
ccplot %>% filter(Cat_Interaction_20_OHT %in% c(1,20)) %>% filter(wig %in% c("Bru-seq")) %>% group_by(Cat_Interaction_20_OHT) %>% nest() %>% pull(data) %>% .[[1]] %>% wilcox.test(value ~ Condition,data=.) %>% .$p.value
ccplot %>% filter(Cat_Interaction_20_OHT %in% c(1,20)) %>% filter(wig %in% c("Bru-seq")) %>% group_by(Cat_Interaction_20_OHT) %>% nest() %>% pull(data) %>% .[[2]] %>% wilcox.test(value ~ Condition,data=.) %>% .$p.value


p2bis <- cowplot::plot_grid(plotlist = list(
  ccplot %>% filter(Cat_Interaction_20_DIVA %in% c(1,20))%>% filter(wig %in% c("BruSeq_150vsU","BruSeq_240vsU","BruSeq_30vsU","BruSeq_60vsU")) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_DIVA),y=value)) + geom_boxplot()  + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic(),
  ccplot %>% filter(Cat_Interaction_100_DIVA %in% c(1,100)) %>% filter(wig %in% c("BruSeq_150vsU","BruSeq_240vsU","BruSeq_30vsU","BruSeq_60vsU")) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_DIVA),y=value)) + geom_boxplot()  + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic(),
  ccplot %>% filter(!is.na(CatVal_DIVA)) %>% filter(wig %in% c("BruSeq_150vsU","BruSeq_240vsU","BruSeq_30vsU","BruSeq_60vsU")) %>% ggplot(aes(x=as.factor(CatVal_DIVA),y=value)) + geom_boxplot()  + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic(),
  ccplot %>% filter(Cat_Interaction_20_OHT %in% c(1,20))%>% filter(wig %in% c("BruSeq_150vsU","BruSeq_240vsU","BruSeq_30vsU","BruSeq_60vsU")) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_OHT),y=value)) + geom_boxplot()  + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic(),
  ccplot %>% filter(Cat_Interaction_100_OHT %in% c(1,100)) %>% filter(wig %in% c("BruSeq_150vsU","BruSeq_240vsU","BruSeq_30vsU","BruSeq_60vsU")) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_OHT),y=value)) + geom_boxplot()  + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic(),
  ccplot %>% filter(!is.na(CatVal_OHT)) %>% filter(wig %in% c("BruSeq_150vsU","BruSeq_240vsU","BruSeq_30vsU","BruSeq_60vsU")) %>% ggplot(aes(x=as.factor(CatVal_OHT),y=value)) + geom_boxplot()  + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()
),nrow=1)


ccplothistone <- res.boxplot.HiC.histone.bc %>% left_join(dplyr::select(res.HiC,-name1,-DIvA:-OHT),by=c("rowname"="name2"))

hihi <- title <- ggdraw() +
  draw_label(
    "BamCompare",
    fontface = 'bold',
    x = 0,
    hjust = 0
  )
p2 <- plot_grid(hihi,
                cowplot::plot_grid(plotlist = list(
                  ccplothistone %>% filter(Cat_Interaction_20_DIVA %in% c(1,20)) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_DIVA),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "bottom"),
                  ccplothistone %>% filter(Cat_Interaction_100_DIVA %in% c(1,100)) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_DIVA),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "bottom"),
                  ccplothistone %>% filter(!is.na(CatVal_DIVA)) %>% ggplot(aes(x=as.factor(CatVal_DIVA),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "bottom"),
                  ccplothistone %>% filter(Cat_Interaction_20_OHT %in% c(1,20)) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_OHT),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "bottom"),
                  ccplothistone %>% filter(Cat_Interaction_100_OHT %in% c(1,100)) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_OHT),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "bottom"),
                  ccplothistone %>% filter(!is.na(CatVal_OHT)) %>% ggplot(aes(x=as.factor(CatVal_OHT),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "bottom")
                ),nrow=2),nrow=2,rel_heights = c(0.05, 1) )

ccplothistone <- res.boxplot.HiC.histone.bwc %>% left_join(dplyr::select(res.HiC,-name1,-DIvA:-OHT),by=c("rowname"="name2"))

hihi <- title <- ggdraw() +
  draw_label(
    "BwCompare",
    fontface = 'bold',
    x = 0,
    hjust = 0
  )
p3 <- plot_grid(hihi,
                cowplot::plot_grid(plotlist = list(
                  ccplothistone %>% filter(Cat_Interaction_20_DIVA %in% c(1,20)) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_DIVA),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "bottom"),
                  ccplothistone %>% filter(Cat_Interaction_100_DIVA %in% c(1,100)) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_DIVA),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "bottom"),
                  ccplothistone %>% filter(!is.na(CatVal_DIVA)) %>% ggplot(aes(x=as.factor(CatVal_DIVA),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "bottom"),
                  ccplothistone %>% filter(Cat_Interaction_20_OHT %in% c(1,20)) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_OHT),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "bottom"),
                  ccplothistone %>% filter(Cat_Interaction_100_OHT %in% c(1,100)) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_OHT),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "bottom"),
                  ccplothistone %>% filter(!is.na(CatVal_OHT)) %>% ggplot(aes(x=as.factor(CatVal_OHT),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "bottom")
                ),nrow=2),nrow=2,rel_heights = c(0.05, 1) )


#Get exons for RNA seq only
# library(EnsDb.Hsapiens.v75)
ens.exons.domains <- ensembldb::exons(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,filter = AnnotationFilter::GeneIdFilter(genes.domains$gene_id))
ens.exons.domains <- regioneR::filterChromosomes(ens.exons.domains,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.exons.domains) <- paste0("chr",seqlevels(ens.exons.domains))
ens.exons.domains$name <- ens.exons.domains$gene_id

res.boxplot.HiC.RNA <- lapply(wigs_for_HiC[c("RNA_C1_mOHT","RNA_C1_pOHT","RNA_C1_pvsmOHT","RNA_C2_mOHT","RNA_C2_pOHT","RNA_C2_pvsmOHT")],function(mon_wig){
  Get1valMean(Name = "",x = ens.exons.domains,one.w = mon_wig) %>% dplyr::select(-wig)
}) %>% bind_rows(.id="wig")%>% mutate(Condition = str_extract(wig,"mOHT|pOHT|Induced_30min|Uninduced")) %>% mutate(wig = str_remove(wig,"_mOHT|_pOHT|_Induced_30min|_Uninduced"))
res.boxplot.HiC.RNA <- res.boxplot.HiC.RNA %>% group_by(wig,rowname,Condition) %>% summarise(value = mean(value))

ccplot <- res.boxplot.HiC.RNA %>% left_join(dplyr::select(res.HiC,-name1,-DIvA:-OHT),by=c("rowname"="name2"))
# ccplot %>% drop_na() %>% dplyr::filter(wig == "RNA_C1") %>% ggplot(aes(x=as.factor(CatVal),y=value,fill=Condition)) + scale_y_log10() + geom_boxplot()+ facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()
# ccplot %>% drop_na() %>% dplyr::filter(wig != "RNA_C1") %>% ggplot(aes(x=as.factor(CatVal),y=value,fill=Condition))  + geom_boxplot()+ facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()
pRNA <- cowplot::plot_grid(plotlist = list(
  ccplot %>% filter(Cat_Interaction_20_DIVA %in% c(1,20)) %>% dplyr::filter(wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_DIVA),y=value,fill=Condition)) + scale_y_log10() + geom_boxplot()+ facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "bottom"),
  ccplot %>% filter(Cat_Interaction_100_DIVA %in% c(1,100)) %>% dplyr::filter(wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_DIVA),y=value,fill=Condition)) + scale_y_log10() +geom_boxplot()+ facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "bottom"),
  ccplot %>% filter(!is.na(CatVal_DIVA)) %>% dplyr::filter(wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(CatVal_DIVA),y=value,fill=Condition)) + scale_y_log10() +geom_boxplot()+ facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "bottom"),
  ccplot %>% filter(Cat_Interaction_20_OHT %in% c(1,20)) %>% dplyr::filter(wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_OHT),y=value,fill=Condition)) + scale_y_log10() + geom_boxplot()+ facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "bottom"),
  ccplot %>% filter(Cat_Interaction_100_OHT %in% c(1,100)) %>% dplyr::filter(wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_OHT),y=value,fill=Condition)) + scale_y_log10() +geom_boxplot()+ facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "bottom"),
  ccplot %>% filter(!is.na(CatVal_OHT)) %>% dplyr::filter(wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(CatVal_OHT),y=value,fill=Condition)) + scale_y_log10() +geom_boxplot()+ facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "bottom")
),nrow=2)


pRNAbis <-  cowplot::plot_grid(plotlist = list(
  ccplot %>% filter(Cat_Interaction_20_DIVA %in% c(1,20)) %>% dplyr::filter(!wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_DIVA),y=value,fill=Condition)) + geom_boxplot()+ facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic() + theme(legend.position = "none"),
  ccplot %>% filter(Cat_Interaction_100_DIVA %in% c(1,100)) %>% dplyr::filter(!wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_DIVA),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "none"),
  ccplot %>% filter(!is.na(CatVal_DIVA)) %>% dplyr::filter(!wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(CatVal_DIVA),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "none"),
  ccplot %>% filter(Cat_Interaction_20_OHT %in% c(1,20)) %>% dplyr::filter(!wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_OHT),y=value,fill=Condition)) + geom_boxplot()+ facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "none"),
  ccplot %>% filter(Cat_Interaction_100_OHT %in% c(1,100)) %>% dplyr::filter(!wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_OHT),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "none"),
  ccplot %>% filter(!is.na(CatVal_OHT)) %>% dplyr::filter(!wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(CatVal_OHT),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "none")
),nrow=1)


pRNAbis_OHT <-  cowplot::plot_grid(plotlist = list(
  ccplot %>% filter(Cat_Interaction_20_OHT %in% c(1,20)) %>% dplyr::filter(!wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(Cat_Interaction_20_OHT),y=value,fill=Condition)) + geom_boxplot()+ facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "none"),
  ccplot %>% filter(Cat_Interaction_100_OHT %in% c(1,100)) %>% dplyr::filter(!wig %in% c("RNA_C1","RNA_C2")) %>% ggplot(aes(x=as.factor(Cat_Interaction_100_OHT),y=value,fill=Condition)) + geom_boxplot() + facet_wrap(~wig,ncol=1,scales="free_y") + theme_classic()+ theme(legend.position = "none")
),nrow=1)

ccplot %>% filter(Cat_Interaction_100_DIVA %in% c(1,100)) %>% dplyr::filter(wig %in% c("RNA_C1_pvsmOHT")) %>% group_by(Cat_Interaction_100_DIVA) %>% nest() %>% mutate(pval = map_dbl(data,function(x){wilcox.test(x$value)$p.value})) %>% 
  dplyr::select(-data)
x <- ccplot %>% filter(Cat_Interaction_100_DIVA %in% c(1,100)) %>% dplyr::filter(wig %in% c("RNA_C1_pvsmOHT"))
wilcox.test(value~Cat_Interaction_100_DIVA,data=x)$p.value


ccplot %>% filter(Cat_Interaction_100_OHT %in% c(1,100)) %>% dplyr::filter(wig %in% c("RNA_C1_pvsmOHT")) %>% group_by(Cat_Interaction_100_OHT) %>% nest() %>% mutate(pval = map_dbl(data,function(x){wilcox.test(x$value)$p.value})) %>% 
  dplyr::select(-data)

x <- ccplot %>% filter(Cat_Interaction_100_OHT %in% c(1,100)) %>% dplyr::filter(wig %in% c("RNA_C1_pvsmOHT"))
wilcox.test(value~Cat_Interaction_100_OHT,data=x)$p.value


ccplot %>% filter(Cat_Interaction_20_OHT %in% c(1,20)) %>% dplyr::filter(wig %in% c("RNA_C1_pvsmOHT")) %>% group_by(Cat_Interaction_20_OHT) %>% nest() %>% mutate(pval = map_dbl(data,function(x){wilcox.test(x$value)$p.value})) %>% 
  dplyr::select(-data)

x <- ccplot %>% filter(Cat_Interaction_20_OHT %in% c(1,20)) %>% dplyr::filter(wig %in% c("RNA_C1_pvsmOHT"))
wilcox.test(value~Cat_Interaction_20_OHT,data=x)$p.value

pdf("../results/Pol2Exp_gene_gamma_DSB_by_cat_20_100.pdf",height=10,width=12)
print(p1)
print(p2)
print(p3)
print(p2bis)
print(pRNA)
print(pRNAbis)
dev.off()



pdf("../../results/Pol2Exp_gene_gamma_DSB_by_cat_20_100_OHT_only.pdf",height=5,width=6)
print(pbruseq_OHT)
print(pRNAbis_OHT)
dev.off()
