require(tidyverse)
require(plyranges)
require(rtracklayer)
require(cowplot)
require(BSgenome.Hsapiens.UCSC.hg19.masked)
require(regioneR)
seqlens <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 %>% seqlengths()
#Functions
bless80 <- read_bed("/mnt/NAS/DATA/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed") 
DE_DIVA <- PhDfunc::GetDE_DIvA()
domains <- bless80 %>% anchor_center() %>% mutate(width=2000000)
ens.genes <- read.table("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed",sep="\t",h=T) %>% GRanges()
ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))
ens.genes$name <- ens.genes$gene_id

ens.filtered <- ens.genes%>% filter(gene_id %in% DE_DIVA$rowname)

direct.damaged <- ens.filtered %>% filter_by_overlaps(bless80)
gamma.genes <- ens.filtered %>% filter_by_overlaps(domains) %>% filter(!gene_id %in% direct.damaged$gene_id) 

ctrl.genes <- regioneR::resampleRegions(c(direct.damaged,gamma.genes),universe = ens.filtered%>% filter_by_non_overlaps(domains))



# domains <- bless80 %>% anchor_center() %>% mutate(width=1000000)


cc.genes <- ens.genes %>% filter(gene_id %in% DE_DIVA$rowname) %>% filter_by_overlaps(domains) %>% 
  join_overlap_left(domains)

cc.genes <- split(cc.genes,cc.genes$name.y)
cc.genes <- lapply(cc.genes,function(x){
  dsb <- bless80 %>% filter(name == unique(x$name.y))
  y <- promoters(x,1,1)
  x %>% mutate(distanceToDSB = distance(y,dsb))
})

cc.genes <- cc.genes %>% map(as_tibble) %>% bind_rows() %>% 
  mutate(GroupDistance = case_when(
    distanceToDSB <= 100000 ~ 1,
    distanceToDSB > 100000 & distanceToDSB <= 200000~ 2,
    distanceToDSB > 200000 & distanceToDSB <= 300000~ 3,
    distanceToDSB > 300000 & distanceToDSB <= 400000~ 4,
    distanceToDSB > 400000 & distanceToDSB <= 500000~ 5,
    distanceToDSB > 500000 & distanceToDSB <= 600000~ 6,
    distanceToDSB > 600000 & distanceToDSB <= 700000~ 7,
    distanceToDSB > 700000 & distanceToDSB <= 800000~ 8,
    distanceToDSB > 800000 & distanceToDSB <= 900000~ 9,
    distanceToDSB > 900000 ~ 10,
    
  ))  

cc.genes <- cc.genes %>% left_join(DE_DIVA,by = c("gene_id"="rowname"))


#Default plot (all genes)

cc.summary <- cc.genes %>% group_by(GroupDistance) %>% summarise(gene_n = dplyr::n(),meanFC = mean(logFC))
p2 <- cc.summary %>% filter(GroupDistance !=1)%>% ggplot(aes(x=GroupDistance,y=meanFC)) +
  geom_point(aes(size=gene_n)) +geom_smooth(aes(x=GroupDistance,y=meanFC),method = "lm") +
  scale_x_continuous(breaks = c(2,3,4,5,6,7,8,9,10),labels = c("100-200kb","200-300kb","300-400kb","400-500kb","500-600kb","600-700kb","700-800kb","800-900kb","900-1Mb"))


pdf("../results/testplot_ianelli_gene_distance_log2FC.pdf",width=10,height=6)
print(p2)
dev.off()

#Ianelli plot (rerun until l55 before running theses lines)
cc.genes <- cc.genes %>% map(as_tibble) %>% bind_rows() %>%
  mutate(GroupDistance = case_when(
    distanceToDSB <= 100 ~ "0.1kb",
    distanceToDSB > 100 & distanceToDSB <= 1000 ~ "1kb",
    distanceToDSB > 1000 & distanceToDSB <= 10000 ~ "10kb",
    distanceToDSB > 10000 & distanceToDSB <= 100000 ~ "100kb",
    distanceToDSB > 100000 ~ "1Mb"
    
  ))
cc.genes <- cc.genes %>%
  mutate(GroupDistance = case_when(
    gene_id %in% direct.damaged$gene_id ~ "0.1kb",
    TRUE ~ GroupDistance
  )) %>%
  mutate(GroupDistance = factor(GroupDistance,levels = c("0.1kb","1kb","10kb","100kb","1Mb")))

cc.genes <- cc.genes %>% left_join(DE_DIVA,by = c("gene_id"="rowname"))
cc.summary <- cc.genes %>% group_by(GroupDistance) %>% summarise(gene_n = dplyr::n(),meanFC = mean(logFC))
p <- cc.summary %>% ggplot(aes(x=GroupDistance,y=meanFC)) +
  geom_point(aes(size=gene_n)) +geom_smooth(aes(x=as.numeric(GroupDistance),y=meanFC),method = "lm")
library(gridExtra)
library(grid)
pdf("../results/testplot_ianelli_gene_distance_log2FC.pdf")

grid.table(cc.summary)
print(p)
dev.off()