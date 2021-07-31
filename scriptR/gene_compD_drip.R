require(tidyverse)
require(plyranges)
require(BSgenome.Hsapiens.UCSC.hg19)
require(TFBSTools)
require(motifmatchr)
require(regioneR)
DE_DIVA <- PhDfunc::GetDE_DIvA()

DE_DIVA <- DE_DIVA %>% 
  mutate(Type = case_when(
    # logFC < 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Downregulated",
    logFC < 0 & FILTER.FC == 1 ~ "Downregulated",
    # logFC < -0.15 ~ "Downregulated",
    # logFC > 0 & FILTER.FC == 1 & FILTER.P == 1~ "Upregulated",
    logFC > 0 & FILTER.FC == 1~ "Upregulated",
    # logFC > 0.15 ~ "Upregulated",
    TRUE ~ "None"
  )) 

ens.genes <- read.table("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed",sep="\t",h=T) %>% GRanges()
ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))
mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
  arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 1000000)
# chr.to.study <- glue::glue("chr{c(2,6,9,13,18,1,17,20,'X')}")
chr.to.study <- glue::glue("chr{c(1,17,'X')}")
# chr.to.study <- glue::glue("chr{c(1,17)}")
genes_DIVA <- ens.genes %>% plyranges::filter(gene_id %in% DE_DIVA$rowname)%>%
  filter_by_non_overlaps(gamma_region)


DE_DIVA_INFO <- DE_DIVA %>% dplyr::select(rowname,logFC,Type)

compD.bw <- c(
  "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHT.bw",
  "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipA.bw",
  "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipB.bw"
) %>%
  setNames(str_remove_all(basename(.),"PC1_all_chr_log2ratio_100kb_|.bw")) %>% 
  map(import.bw) %>%
  map(filter,seqnames %in% chr.to.study) %>% 
  map(as_tibble) %>% 
  map(dplyr::select,-strand,-width) %>% 
  purrr::reduce(left_join,by=c("seqnames","start","end")) %>% 
  mutate_at(vars(starts_with("score")),function(x){
    as.numeric(x>0)
  }) %>% 
  mutate(score = rowSums(across(starts_with("score")))) %>% filter(score == 3)


mes_genes.PC1 <- genes_DIVA %>%
  # filter(seqnames %in% chr.to.study) %>% 
  filter_by_overlaps(as_granges(compD.bw)) %>%
  as_tibble() %>% 
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(TypeD = glue::glue("CompD"))

mes_genes.PC1_nocompD <- genes_DIVA %>%
  filter(seqnames %in% chr.to.study) %>%
  filter_by_non_overlaps(as_granges(compD.bw)) %>%
  as_tibble() %>% 
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(TypeD = glue::glue("NoCompD"))

mes_genes <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD) %>% as_granges()
mes_genes$name <- mes_genes$gene_id
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

mes_bw <- list(
  "DRIP2_C1_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/HYG77BGX2_DRIPseq_SETX/PROCESSED/WIGGLE/hg19/HYG77BGX2_DRIP2_C1_DIVA_17s002460-1-1_Clouaire_lane117s002460_normalized_hg19_nodups.bw",
  "DRIP2_C1_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/HYG77BGX2_DRIPseq_SETX/PROCESSED/WIGGLE/hg19/HYG77BGX2_DRIP2_C1_OHT_17s002461-1-1_Clouaire_lane117s002461_normalized_hg19_nodups.bw"
) %>% map(import.bw,as="RleList")
res.plot <- lapply(names(mes_bw),function(one_wig){
  
  PhDfunc::Get1val(Name = one_wig,one.w = mes_bw[[one_wig]],x = mes_genes %>% promoters(1000,1000))
}) %>% bind_rows()

res <- as_tibble(mes_genes) %>% right_join(res.plot,by = c("name"="rowname"))
p1 <- res %>% ggplot(aes(x=TypeD,y=value,fill=wig)) + geom_boxplot() + scale_y_log10()
p2 <- res %>% ggplot(aes(x=Type,y=value,fill=wig)) + geom_boxplot() + scale_y_log10() + facet_wrap(~TypeD)
pdf("../results/Drip_seq_compD_ABD_genes.pdf",height=6,width=8)
print(p1)
print(p2)
dev.off()
