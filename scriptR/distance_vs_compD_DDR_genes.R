require(tidyverse)
require(TFBSTools)
require(motifmatchr)
require(plyranges)
require(BSgenome.Hsapiens.UCSC.hg19)
# PCs.ratio.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio.bw"
# PCs.ratio.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHT.bw"
PCs.ratio.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHTDNAPKi.bw"
PCs.ratio <- import.bw(PCs.ratio.path)
out_name_file <- basename(PCs.ratio.path) %>% str_remove(".bw")
PCs.ratio$name <- str_c("bin",1:length(PCs.ratio))

chr.to.study <- keep.chr <- glue::glue("chr{c(2,6,9,13,18,1,17,20,'X')}")
# chr.to.study <- keep.chr <- glue::glue("chr{c(1,17,'X')}")
# chr.to.study <- keep.chr <- glue::glue("chr{c(1)}")

PCs.ratio <- PCs.ratio %>% filter(seqnames %in% chr.to.study)
mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
  arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 1000000)

PCs.ratio <- PCs.ratio %>% filter_by_non_overlaps(gamma_region)

#Gene compD vs gene noCompD

cc.cov <- PCs.ratio.path %>% import.bw(as="RleList")
DE_DIVA <- PhDfunc::GetDE_DIvA()

DE_DIVA <- DE_DIVA %>% 
  mutate(Type = case_when(
    # logFC < 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Downregulated",
    # logFC < 0 & FILTER.FC == 1~ "Downregulated",
    logFC < 0 ~ "Downregulated",
    # logFC > 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Upregulated",
    # logFC > 0 & FILTER.FC == 1~ "Upregulated",
    logFC > 0~ "Upregulated",
    TRUE ~ "None"
  )) 

ens.genes <- read.table("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed",sep="\t",h=T) %>% GRanges()
ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))
mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
  arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
# gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 1000000)
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 2000000)
genes_of_interest <- ens.genes %>% filter(gene_id %in% DE_DIVA$rowname)%>%
  filter_by_non_overlaps(gamma_region) %>% 
  regioneR::filterChromosomes(keep.chr = keep.chr)


Get1val <- function (Name, one.w, x) 
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


mes_genes <- list("Genes"=genes_of_interest
)

cc.dsb <- mes_genes %>% 
  map(function(one_dsb){
    one_dsb$name <- one_dsb$gene_id
    Get1val(Name = "cc",one.w = cc.cov,x = one_dsb)
  }) %>% bind_rows(.id = "Name")

cc.dsb <- cc.dsb %>%
  mutate(Group.PC1 = case_when(value < 0 ~ "noCompD",value > 0 ~"compD")) 

##Distance to nearest DSB

res_distance <- mes_genes %>% map(function(x){
  sres <- x %>% distanceToNearest(mes_DSB)
  x[queryHits(sres)] %>% mutate(distance = mcols(sres)$distance) %>% as_tibble()
})%>% bind_rows(.id = "Name")


res_f <- res_distance %>% left_join(dplyr::select(DE_DIVA,rowname,Type),by = c("gene_id"="rowname") ) %>% left_join(cc.dsb,by = c("gene_id"="rowname") )  %>% drop_na()


p1 <- res_f%>% ggplot(aes(x=log10(distance),y=value,col=Type)) + geom_point() + facet_grid(Type~Group.PC1)
p1.bis <-  res_f%>% ggplot(aes(x=log10(distance),y=value,col=Group.PC1)) + geom_point() + facet_grid(Type~.)
pf <- p1/p1.bis + plot_layout(guide="collect")

pb1 <- res_f %>% filter(Type == "Upregulated") %>% ggplot(aes(x=Group.PC1,y=(distance),fill=Group.PC1)) + geom_boxplot(outlier.shape=NA) + geom_jitter(alpha=0.2) +ggtitle("Upregulated Genes")
pb2 <- res_f %>% filter(Type == "Upregulated") %>% ggplot(aes(x=Group.PC1,y=log10(distance),fill=Group.PC1)) + geom_boxplot(outlier.shape=NA) + geom_jitter(alpha=0.2)
pf2 <- pb1/pb2 + plot_layout(guide="collect")
pdf("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/compD_value_vs_DSB_distance.pdf",height=12,width=16)
print(pf)
print(pf2)
dev.off()


#Nombre de gènes à chaque distance 
cc <- res_f %>% filter(Type == "Upregulated") %>% mutate(distance_tile = ntile(distance,20)) %>% ggplot(aes(x=distance_tile,y=value)) +geom_bar(stat="summary",fun="mean")

res_f %>% dplyr::select(-Name.x,-Name.y,-seq_coord_system)%>% write_tsv("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/genes_compD_value_vs_DSB_distance.tsv")

# corriger A ET B pour saddle + tsv + françois plot
# voir si on voit toujours l'effet up regulated avec fort niveau de D par rapport aux autres + motifs en p53
# + effet distance seulement ceux au dessus de 0 dans le treshold