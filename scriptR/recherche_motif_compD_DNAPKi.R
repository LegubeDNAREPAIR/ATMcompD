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

PCs.ratio <- PCs.ratio %>% filter(seqnames %in% chr.to.study)
mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
  arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 1000000)

PCs.ratio <- PCs.ratio %>% filter_by_non_overlaps(gamma_region)
PCs.ratio <- PCs.ratio %>% mutate(Group.PC1 = ntile(score,100))
PCs.ratio.seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19,PCs.ratio)

MotifList <- readJASPARMatrix("../data/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt")

#Gene compD vs gene noCompD

cc.cov <- PCs.ratio.path %>% import.bw(as="RleList")
DE_DIVA <- PhDfunc::GetDE_DIvA()

DE_DIVA <- DE_DIVA %>% 
  mutate(Type = case_when(
    logFC < 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Downregulated",
    # logFC < 0 & FILTER.FC == 1~ "Downregulated",
    logFC > 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Upregulated",
    # logFC > 0 & FILTER.FC == 1~ "Upregulated",
    TRUE ~ "None"
  )) 

ens.genes <- read.table("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed",sep="\t",h=T) %>% GRanges()
ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))
mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
  arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 1000000)
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



#matchMotifs function
# up.distance <- 1000
up.distance <- 500
down.distance <- 0
up.compD <- cc.dsb %>% filter(Group.PC1=="compD")
down.compD <- cc.dsb %>% filter(Group.PC1=="noCompD")
# down.compD <- cc.dsb %>% filter(Group.PC1=="Everything_else")
genes.up <- ens.genes[ens.genes$gene_id %in% c(as.vector(up.compD$rowname),as.vector(down.compD$rowname))] %>% 
  promoters(up.distance,down.distance)
# genes.up$Group.PC1 <- c(rep("TopCompD",nrow(up.compD)),rep("Everything_else",nrow(down.compD)))
genes.up$Group.PC1 <- c(rep("compD",nrow(up.compD)),rep("noCompD",nrow(down.compD)))

genes.up.seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19,genes.up)
names(genes.up.seq) <- genes.up$gene_name
# genes.up.seq[genes.up$Group.PC1 =="CompD"] %>% writeXStringSet(glue::glue("../results/DNA_sequences_in_compD_up_genes_{out_name_file}_{up.distance}_{down.distance}.fa"))
# genes.up.seq[!genes.up$Group.PC1 =="CompD"] %>% writeXStringSet(glue::glue("../results/DNA_sequences_in_NOcompD_up_genes_{out_name_file}_{up.distance}_{down.distance}.fa"))
ix2 = motifmatchr::matchMotifs(MotifList,genes.up.seq,out = "scores")

#count number of motifs in each sequence
mat2 = motifmatchr::motifCounts(ix2)
hihi = motifmatchr::motifScores(ix2)
colnames(mat2) <- names(MotifList)
rownames(mat2) <- paste(genes.up)

countMat <- mat2 %>% as.matrix() %>% reshape2::melt()
genes.up.count <- genes.up %>% as_tibble() %>% mutate(pos = paste(genes.up)) %>% right_join(countMat, by = c("pos"="Var1")) %>% dplyr::rename(Motif="Var2",count=value)

genes.up.count.group <- genes.up.count %>% group_by(Group.PC1,Motif) %>% summarise(count = sum(count))

genes.up.count.group12 <- genes.up.count.group %>% spread(key = Group.PC1,value = count)
total <- colSums(genes.up.count.group12[,2:3])
genes.up.count.group12$stat.fisher.universe <- apply(genes.up.count.group12[,2:3],1,function(x){fisher.test(matrix(c(x[[1]],x[[2]],total-x),byrow=T,ncol=2),alternative = "greater")$p.value})
genes.up.count.group12$adj.fisher.universe <- p.adjust(genes.up.count.group12$stat.fisher.universe,method = "BH")

#Final table
genes.up.count.final <- genes.up.count.group12 %>%
  # mutate(total_Low = sum(`Everything_else`),total_High = sum(`TopCompD`)) %>% mutate(Enrichment = (`TopCompD`/total_High)/(`Everything_else`/total_Low)) %>%
  mutate(total_Low = sum(`noCompD`),total_High = sum(`compD`)) %>% mutate(Enrichment = (`compD`/total_High)/(`noCompD`/total_Low)) %>%
  filter(adj.fisher.universe < 0.4) %>%
  dplyr::select(-total_Low,-total_High) %>%
  dplyr::select(Motif,Enrichment,everything())

genes.up.count.final %>%
  filter(is.finite(Enrichment)) %>% 
  ggplot(aes(x=fct_reorder(Motif,Enrichment),y=Enrichment,col=-log10(adj.fisher.universe),size=Enrichment)) + geom_point() +coord_flip() +
  xlab("Motifs")+scale_color_viridis_c() +
  # ggtitle(glue::glue("Motif enrichment of top100 compD vs top1 compD genes\nupstream: {up.distance}bp and downstream: {down.distance}bp"))-> pmotifs
  ggtitle(glue::glue("Motif enrichment of compD vs nocompD genes\nupstream: {up.distance}bp and downstream: {down.distance}bp"))-> pmotifs
ggsave(glue::glue("../results/Motifs_enrichment_in_compDvsnocompD_genes_{out_name_file}_nb_genes_{nrow(up.compD)}_{nrow(down.compD)}.pdf"),pmotifs,height=10,width=10)
