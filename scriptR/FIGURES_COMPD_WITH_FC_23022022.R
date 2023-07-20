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

require(parallel)
require(tidyverse)
require(plyranges)
require(BSgenome.Hsapiens.UCSC.hg19)
require(regioneR)
DE_DIVA <- PhDfunc::GetDE_DIvA()

DE_DIVA <- DE_DIVA %>% 
  mutate(Type = case_when(
    # logFC < 0 & FILTER.P == 1 ~ "Downregulated",
    # logFC < 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Downregulated",
    # logFC < 0 & FILTER.FC == 1 ~ "Downregulated",
    logFC < -0.3 ~ "Downregulated",
    # logFC < -0.15 ~ "Downregulated",
    # logFC > 0 & FILTER.FC == 1 & FILTER.P == 1~ "Upregulated",
    # logFC > 0 & FILTER.P == 1~ "Upregulated",
    # logFC > 0 & FILTER.FC == 1~ "Upregulated",
    # logFC > 0.15 ~ "Upregulated",
    logFC > 0.3 ~ "Upregulated",
    TRUE ~ "None"
  )) 

ens.genes <- read.table("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed",sep="\t",h=T) %>% GRanges()
ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))
mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
  arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 2000000)
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
  map(plyranges::filter,seqnames %in% chr.to.study) %>% 
  map(as_tibble) %>% 
  map(dplyr::select,-strand,-width) %>% 
  bind_rows(.id="Name") %>% 
  mutate(binscore = as.numeric(score>0)) %>% 
  pivot_wider(names_from = c(Name),values_from = c(score,binscore)) %>% 
  mutate(score = rowSums(across(starts_with("binscore")))) 
# %>%  filter(score == 3)
# %>%  filter(score == 3)



DNAPKichr <- glue::glue("chr{c(2,6,9,13,18,1,17,20,'X')}")
compD.bw.DNAPKi <- c(
  "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHTDNAPKi.bw"
) %>%
  setNames(str_remove_all(basename(.),"PC1_all_chr_log2ratio_100kb_|.bw")) %>% 
  map(import.bw) %>%
  map(plyranges::filter,seqnames %in% DNAPKichr) %>% 
  map(as_tibble) %>% 
  map(dplyr::select,-strand,-width) %>% 
  bind_rows(.id="Name") %>% 
  mutate(binscore = as.numeric(score>0)) %>% 
  pivot_wider(names_from = c(Name),values_from = c(score,binscore)) %>% 
  mutate(score = rowSums(across(starts_with("binscore"))))

# mes_genes <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD) %>% as_granges()
mes_genes <- genes_DIVA
mes_genes$name <- mes_genes$gene_id



#### TRANSLOC

mes_genes.PC1 <- genes_DIVA %>%
  filter(seqnames %in% chr.to.study) %>%
  filter_by_overlaps(as_granges(compD.bw) %>% filter(score==3)) %>%
  as_tibble() %>%
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("CompD_{Type}"))


mes_genes.PC1_nocompD <- genes_DIVA %>%
  filter(seqnames %in% chr.to.study) %>%
  filter_by_overlaps(as_granges(compD.bw) %>% filter(score==0)) %>%
  as_tibble() %>%
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("NoCompD_{Type}")) 


mes_genes.PC1_nocompD.duplicate <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD) %>% dplyr::count(gene_id) %>% filter(n>1) %>% pull(gene_id)

mes_genes.PC1 <- mes_genes.PC1 %>% filter(!gene_id %in% mes_genes.PC1_nocompD.duplicate)
mes_genes.PC1_nocompD <- mes_genes.PC1_nocompD %>% filter(!gene_id %in% mes_genes.PC1_nocompD.duplicate)




A.list <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD) 

require(regioneR)


A.list.compD.nocompD <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD) %>% separate(Type,into=c("compD","Type")) %>% 
  mutate(Type=case_when( 
    Type == "Upregulated" ~ "Upregulated",
    TRUE ~ "NONE")) %>% 
  group_by(compD,Type) %>% nest()

my_transloc.GR <- data_files_breakpoints.GR[[1]]

A.list.compD.nocompD.res <- A.list.compD.nocompD %>% mutate(permtest = map(data,
                                                                           function(x){
                                                                             A <-as_granges(x)
                                                                             set.seed(1234)
                                                                             pt <- permTest(A=A, ntimes=500, randomize.function=resampleRegions, universe=filter(genes_DIVA,seqnames %in% glue::glue("chr{c(1,17,'X')}")),
                                                                                            evaluate.function=numOverlaps, B=my_transloc.GR, verbose=F)
                                                                             
                                                                             tibble(
                                                                               n = length(A),
                                                                               permuted = mean(pt$numOverlaps$permuted),
                                                                               observed = pt$numOverlaps$observed,
                                                                               pval = format.pval(pt$numOverlaps$pval,2)
                                                                             )
                                                                             
                                                                           })) %>% dplyr::select(-data) %>% unnest()
my.p.bar <- A.list.compD.nocompD.res %>% 
  mutate(observed = observed/permuted) %>% 
  mutate(permuted = permuted/permuted) %>% 
  gather(key=permtype,value = value,-compD,-Type,-pval,-n) %>% 
  mutate(mycolor = case_when(
    permtype == "permuted" ~ "#7f8c8d",
    permtype == "observed" & value > 1 ~ "#27ae60",
    permtype == "observed" & value < 1 ~ "#27ae60"
  )) %>% 
  mutate(Type = glue::glue("{Type} ({n})")) %>% 
  ggplot(aes(x=Type,y=value,fill=mycolor)) + geom_bar(stat="identity",position = "dodge",col="black") + theme_classic() + theme(axis.text.x = element_markdown()) +
  scale_fill_identity() + scale_y_continuous(labels = scales::percent,limits = c(0,2.6)) + facet_wrap(~compD,scales="free_x")
ggsave("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/after_FC03/barplot_permutatio_observed_permuted_percentage_manipABD_compD_nocompD_bygenetype.pdf",my.p.bar,width=8,height=4)


#### TRANSLOC DNAPKi

mes_genes.PC1.dnapki <- genes_DIVA %>%
  filter(seqnames %in% DNAPKichr) %>%
  filter_by_overlaps(as_granges(compD.bw.DNAPKi) %>% filter(score==1)) %>%
  as_tibble() %>%
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("CompD_{Type}"))


mes_genes.PC1_nocompD.dnapki <- genes_DIVA %>%
  filter(seqnames %in% DNAPKichr) %>%
  filter_by_overlaps(as_granges(compD.bw.DNAPKi) %>% filter(score==0)) %>%
  as_tibble() %>%
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("NoCompD_{Type}")) 


mes_genes.PC1_nocompD.duplicate <- rbind(mes_genes.PC1.dnapki,mes_genes.PC1_nocompD.dnapki) %>% dplyr::count(gene_id) %>% filter(n>1) %>% pull(gene_id)

mes_genes.PC1.dnapki <- mes_genes.PC1.dnapki %>% filter(!gene_id %in% mes_genes.PC1_nocompD.duplicate)
mes_genes.PC1_nocompD.dnapki <- mes_genes.PC1_nocompD.dnapki %>% filter(!gene_id %in% mes_genes.PC1_nocompD.duplicate)



A.list.compD.nocompD.dnapki <- rbind(mes_genes.PC1.dnapki,mes_genes.PC1_nocompD.dnapki) %>% separate(Type,into=c("compD","Type")) %>% 
  mutate(Type=case_when( 
    Type == "Upregulated" ~ "Upregulated",
    TRUE ~ "NONE")) %>% 
  group_by(compD,Type) %>% nest()

my_transloc.GR <- data_files_breakpoints.GR[[1]]

A.list.compD.nocompD.res.dnapki <- A.list.compD.nocompD.dnapki %>% mutate(permtest = map(data,
                                                                           function(x){
                                                                             A <-as_granges(x)
                                                                             set.seed(1234)
                                                                             pt <- permTest(A=A, ntimes=1000, randomize.function=resampleRegions, universe=filter(genes_DIVA,seqnames %in% glue::glue("chr{c(1,17,'X')}")),
                                                                                            evaluate.function=numOverlaps, B=my_transloc.GR, verbose=F)
                                                                             
                                                                             tibble(
                                                                               n = length(A),
                                                                               permuted = mean(pt$numOverlaps$permuted),
                                                                               observed = pt$numOverlaps$observed,
                                                                               pval = format.pval(pt$numOverlaps$pval,2)
                                                                             )
                                                                             
                                                                           })) %>% dplyr::select(-data) %>% unnest()


my.p.bar_Dnapki <- A.list.compD.nocompD.res.dnapki %>% 
  mutate(observed = observed/permuted) %>% 
  mutate(permuted = permuted/permuted) %>% 
  gather(key=permtype,value = value,-compD,-Type,-pval,-n) %>% 
  mutate(mycolor = case_when(
    permtype == "permuted" ~ "#7f8c8d",
    permtype == "observed" & value > 1 ~ "#27ae60",
    permtype == "observed" & value < 1 ~ "#27ae60"
  )) %>% 
  mutate(Type = glue::glue("{Type} ({n})")) %>% 
  ggplot(aes(x=Type,y=value,fill=mycolor)) + geom_bar(stat="identity",position = "dodge",col="black") + theme_classic() + theme(axis.text.x = element_markdown()) +
  scale_fill_identity() + scale_y_continuous(labels = scales::percent,limits = c(0,2.6)) + facet_wrap(~compD,scales="free_x")
ggsave("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/after_FC03/barplot_permutatio_observed_permuted_percentage_manipABD_compD_nocompD_bygenetype.pdf",my.p.bar,width=8,height=4)


# # A tibble: 6 × 6
# # Groups:   compD, Type [6]
# compD   Type              n permuted observed pval 
# <chr>   <chr>         <int>    <dbl>    <int> <chr>
#   1 CompD   None            449    312.       323 0.35 
# 2 CompD   Downregulated    32     22.3       18 0.36 
# 3 CompD   Upregulated      45     31.0       66 0.001
# 4 NoCompD None           1157    804.       795 0.41 
# 5 NoCompD Upregulated     117     81.0       66 0.17 
# 6 NoCompD Downregulated    74     51.6       34 0.072


##FIG 4B

mes_bw.list <- list(
  "DRIP2_C1_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/HYG77BGX2_DRIPseq_SETX/PROCESSED/WIGGLE/hg19/HYG77BGX2_DRIP2_C1_DIVA_17s002460-1-1_Clouaire_lane117s002460_normalized_hg19_nodups.bw"
  ,"DRIP2_C1_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/HYG77BGX2_DRIPseq_SETX/PROCESSED/WIGGLE/hg19/HYG77BGX2_DRIP2_C1_OHT_17s002461-1-1_Clouaire_lane117s002461_normalized_hg19_nodups.bw"
  # ,"qDRIP_mOHT_HWGF7BGXJ"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/qDRIP/Clouaire_HWGF7BGXJ_4h_24h_OHT_ALINE/PROCESSED/mapping/BIGWIG/BIGWIG/HWGF7BGXJ_qDRIP_4_24_21s004024-1-1_Clouaire_lane1qDRIPDIVA_total_normalizedReadCount.bw"
  # ,"qDRIP_pOHT_HWGF7BGXJ"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/qDRIP/Clouaire_HWGF7BGXJ_4h_24h_OHT_ALINE/PROCESSED/mapping/BIGWIG/BIGWIG/HWGF7BGXJ_qDRIP_4_24_21s004024-1-1_Clouaire_lane1qDRIPOHT_total_normalizedReadCount.bw"
  # ,"qDRIP_mOHT_HH7TKBGXK"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/qDRIP/Clouaire_HH7TKBGXK_SETX_FTO_FLORIAN/Clouaire_HH7TKBGXK_SETX_FTO/PROCESSED/mapping/BIGWIG/BIGWIG/HH7TKBGXK_qDRIP_siC_siS_siFTO_21s004481-1-1_Clouaire_lane1qDRIPsiCtrlDIVA_total_normalizedReadCount.bw"
  # ,"qDRIP_pOHT_HH7TKBGXK"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/qDRIP/Clouaire_HH7TKBGXK_SETX_FTO_FLORIAN/Clouaire_HH7TKBGXK_SETX_FTO/PROCESSED/mapping/BIGWIG/BIGWIG/HH7TKBGXK_qDRIP_siC_siS_siFTO_21s004481-1-1_Clouaire_lane1qDRIPsiCtrlOHT_total_normalizedReadCount.bw"
  # ,"DRIP2_SETX_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/HYG77BGX2_DRIPseq_SETX/PROCESSED/WIGGLE/hg19/HYG77BGX2_DRIP2_STX_DIVA_17s002462-1-1_Clouaire_lane117s002462_normalized_hg19_nodups.bw",
  # ,"DRIP2_SETX_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/HYG77BGX2_DRIPseq_SETX/PROCESSED/WIGGLE/hg19/HYG77BGX2_DRIP2_STX_OHT_17s002463-1-1_Clouaire_lane117s002463_normalized_hg19_nodups.bw"
) 
mes_bw <- mes_bw.list %>% map(import.bw,as="RleList")
res.plot <- mclapply(names(mes_bw),function(one_wig){
  # my_x <- mes_genes %>% promoters(500,500)
  my_x <- mes_genes %>% anchor_start() %>% plyranges::stretch(1000)%>% anchor_end() %>% plyranges::stretch(1000)
  dat1 <- Get1valMean(Name = one_wig,one.w = mes_bw[[one_wig]],x = my_x)
  my_x <- mes_genes %>% promoters(1000,1000)
  dat2 <- Get1valMean(Name = one_wig,one.w = mes_bw[[one_wig]],x = my_x)
  rbind(
    dat1 %>% mutate(GeneType = "GeneBody"),
    dat2 %>% mutate(GeneType = "Promoters")
  )
},mc.cores=length(mes_bw)) %>% bind_rows()


res <- as_tibble(mes_genes) %>% left_join(DE_DIVA_INFO,by = c("gene_id"="rowname"))%>% right_join(res.plot,by = c("name"="rowname"))



theme_set(theme_classic(base_size=12))
theme_update(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
##First plot -> 3 ABD with 1,17,X

res4 <- 
  rbind(
    res %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)) %>% as_tibble() %>% mutate(TypeD = glue::glue("CompD")),
    res %>% as_granges() %>% filter_by_non_overlaps(as_granges(compD.bw)) %>% as_tibble() %>% mutate(TypeD = glue::glue("NoCompD")) %>% filter(seqnames %in% chr.to.study)
  ) %>% filter(GeneType == "GeneBody")


res4.duplicate <- res4 %>% dplyr::count(gene_id,wig,GeneType) %>% filter(n>1) %>% pull(gene_id) %>% unique()

res4 <- res4 %>% filter(!gene_id %in% res4.duplicate)

p7 <- res4 %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() + scale_y_log10()  + facet_wrap(~GeneType)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))

pdf("../../results/after_FC03/Drip_seq_compD_ABD_genes_14012022.pdf",height=6,width=5)
print(p7)
dev.off()

##PVALUES
numbers <- res4 %>% dplyr::count(wig,TypeD) %>% spread(key = TypeD,value=n)
pvalues.dat <- res4 %>% group_by(wig) %>% nest() %>% 
  mutate(data = map(data,function(x){
    wilcox.test(value~TypeD,data=x) %>% broom::tidy()
  })) %>% unnest(data) %>% mutate(p.value = format.pval(p.value,3))
pvalues.dat %>% left_join(numbers,by="wig")
# # Groups:   wig [2]
# wig           statistic p.value method                                            alternative CompD NoCompD
# <chr>             <dbl> <chr>   <chr>                                             <chr>       <int>   <int>
#   1 DRIP2_C1_mOHT   2601099 0.00276 Wilcoxon rank sum test with continuity correction two.sided     526    9179
# 2 DRIP2_C1_pOHT   2577032 0.00912 Wilcoxon rank sum test with continuity correction two.sided     526    9179
## FIG4C


###sort by qDRIP signal and then plot compD
##NE BOUGE PAS

qDRIPgenes <- res %>% filter(seqnames %in% chr.to.study)

gene_list <- qDRIPgenes %>% pull(gene_name) %>% unique


compD.signal <- c(
  "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHT.bw",
  "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipA.bw",
  "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipB.bw",
  # "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHTDNAPKi.bw"
  "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/compD/PC1_all_chr_log2ratio_100kb_HIC_siSETX_OHT_HIC_siSETX_DIVA.bw"
) %>%
  setNames(str_remove_all(basename(.),"PC1_all_chr_log2ratio_100kb_|.bw")) %>% 
  map(import.bw,as="RleList")

res.compD <- mclapply(names(compD.signal),function(one_wig){
  # my_x <- mes_genes %>% promoters(500,500)
  my_x <- mes_genes %>% filter(gene_name %in% gene_list) %>% anchor_start() %>% plyranges::stretch(1000)%>% anchor_end() %>% plyranges::stretch(1000)
  dat1 <- Get1valMean(Name = one_wig,one.w = compD.signal[[one_wig]],x = my_x)
  dat1
},mc.cores=length(compD.signal)) %>% bind_rows() %>% spread(key=wig,value=value)


qDRIPgenes.compD <- qDRIPgenes %>% left_join(res.compD,by=c("gene_id"="rowname"))

qDRIPgenes.compD <- qDRIPgenes.compD %>%  group_by(GeneType,wig) %>% arrange(desc(value)) %>% mutate(Group = ntile(value,100))






pDRIP.compd <- qDRIPgenes.compD %>%
  filter(GeneType == "GeneBody") %>% 
  mutate(Group = case_when(
    Group %in% 1:10 ~ "Bottom10",
    Group %in% 91:100 ~ "Top10",
    TRUE ~ "Other"
  )) %>% filter(Group != "Other") %>% dplyr::select(seqnames,gene_name,wig:Group) %>% gather(key=CompD_Source,value=compD,-seqnames:-GeneType,-Group)

plotqDRIP14 <- pDRIP.compd %>% ggplot(aes(x=GeneType,y=compD,fill=Group)) + geom_boxplot() + facet_grid(CompD_Source~wig) + scale_fill_manual(values=MetBrewer::met.brewer("Tsimshian", 2))


pdf("../../results/after_FC03/compD_by_Drip_Cat_14012022.pdf",height=12,width=8)
print(plotqDRIP14 + ggtitle("10"))
dev.off()


numbers <- pDRIP.compd %>% group_by(CompD_Source,wig,Group) %>% dplyr::count() %>% spread(key=Group,value=n)

pvalues.dat <- pDRIP.compd %>% group_by(CompD_Source,wig) %>% nest() %>% 
  mutate(data = map(data,function(x){
    wilcox.test(compD~Group,data=x) %>% broom::tidy()
  })) %>% unnest(data) %>% mutate(p.value = format.pval(p.value,2))
pvalues.dat %>% left_join(numbers,by=c("wig","CompD_Source"))
# 
# # A tibble: 8 × 8
# # Groups:   wig, CompD_Source [8]
# wig           CompD_Source                   statistic p.value method              alternative Bottom10 Top10
# <chr>         <chr>                              <dbl> <chr>   <chr>               <chr>          <int> <int>
#   1 DRIP2_C1_mOHT HiC_D_OHT                          15932 0.26    Wilcoxon rank sum … two.sided        190   180
# 2 DRIP2_C1_pOHT HiC_D_OHT                          16293 0.43    Wilcoxon rank sum … two.sided        190   180
# 3 DRIP2_C1_mOHT HIC_siSETX_OHT_HIC_siSETX_DIVA     12715 2e-05   Wilcoxon rank sum … two.sided        190   180
# 4 DRIP2_C1_pOHT HIC_siSETX_OHT_HIC_siSETX_DIVA     13232 0.00017 Wilcoxon rank sum … two.sided        190   180
# 5 DRIP2_C1_mOHT OHT_manipA                         12273 2.7e-06 Wilcoxon rank sum … two.sided        190   180
# 6 DRIP2_C1_pOHT OHT_manipA                         14234 0.0053  Wilcoxon rank sum … two.sided        190   180
# 7 DRIP2_C1_mOHT OHT_manipB                         14141 0.004   Wilcoxon rank sum … two.sided        190   180
# 8 DRIP2_C1_pOHT OHT_manipB                         14431 0.0095  Wilcoxon rank sum … two.sided        190   180

## FIG 4e



mes_bw.rna <- list(
  # "RNA_C1_1_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYL5LBGX2/HYL5LBGX2_RNA_C1_OHT_1_17s002465-1-1_Clouaire_lane117s002465_hg19plusERCC_normalized.bw"
  # ,"RNA_STX_1_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYL5LBGX2/HYL5LBGX2_RNA_STX_OHT_1_17s002467-1-1_Clouaire_lane117s002467_hg19plusERCC_normalized.bw"
  # ,"RNA_C1_2_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYMHNBGX2/HYMHNBGX2_RNA_C1_OHT_2_17s002469-1-1_Clouaire_lane117s002469_hg19plusERCC_normalized.bw"
  # ,"RNA_STX_2_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYMHNBGX2/HYMHNBGX2_RNA_STX_OHT_2_17s002471-1-1_Clouaire_lane117s002471_hg19plusERCC_normalized.bw"
  "HYL5LBGX2_RNA_STX_vs_C1_DIVA.log2.bamCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_STX_vs_C1_DIVA.log2.bamCompare.bw"
  ,"HYL5LBGX2_RNA_STX_vs_C1_DIVA.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_STX_vs_C1_DIVA.log2.bigwigCompare.bw"
  ,"HYL5LBGX2_RNA_STX_vs_C1_OHT.log2.bamCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_STX_vs_C1_OHT.log2.bamCompare.bw"
  ,"HYL5LBGX2_RNA_STX_vs_C1_OHT.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_STX_vs_C1_OHT.log2.bigwigCompare.bw"
  # ,"HYMHNBGX2_RNA_STX_vs_C1_DIVA.log2.bamCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYMHNBGX2_RNA_STX_vs_C1_DIVA.log2.bamCompare.bw"
  # ,"HYMHNBGX2_RNA_STX_vs_C1_DIVA.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYMHNBGX2_RNA_STX_vs_C1_DIVA.log2.bigwigCompare.bw"
  # ,"HYMHNBGX2_RNA_STX_vs_C1_OHT.log2.bamCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYMHNBGX2_RNA_STX_vs_C1_OHT.log2.bamCompare.bw"
  # ,"HYMHNBGX2_RNA_STX_vs_C1_OHT.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYMHNBGX2_RNA_STX_vs_C1_OHT.log2.bigwigCompare.bw"
  # ,"HYL5LBGX2_RNA_OHT_vs_DIVA_STX.log2.bamCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_OHT_vs_DIVA_STX.log2.bamCompare.bw"
  # ,"HYL5LBGX2_RNA_OHT_vs_DIVA_STX.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_OHT_vs_DIVA_STX.log2.bigwigCompare.bw"
  # ,"HYL5LBGX2_RNA_OHT_vs_DIVA_CTRL.log2.bamCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_OHT_vs_DIVA_CTRL.log2.bamCompare.bw"
  # ,"HYL5LBGX2_RNA_OHT_vs_DIVA_CTRL.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_OHT_vs_DIVA_CTRL.log2.bigwigCompare.bw"
) 
mes_bw.rna <- mes_bw.rna %>% mclapply(import.bw,as="RleList",mc.cores=length(mes_bw.rna))


ens.exons.domains <- ensembldb::exons(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,filter = AnnotationFilter::GeneIdFilter(mes_genes$gene_id))
ens.exons.domains <- regioneR::filterChromosomes(ens.exons.domains,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.exons.domains) <- paste0("chr",seqlevels(ens.exons.domains))
ens.exons.domains$name <- ens.exons.domains$gene_id


res.plot.rna <- mclapply(names(mes_bw.rna),function(one_wig){
  # my_x <- mes_genes %>% promoters(500,500)
  my_x <- mes_genes
  dat1 <- Get1valMean(Name = one_wig,one.w = mes_bw.rna[[one_wig]],x = my_x)
  my_x <- ens.exons.domains 
  dat2 <- Get1valMean(Name = one_wig,one.w = mes_bw.rna[[one_wig]],x = my_x)%>% group_by(wig,rowname) %>% summarise(value = mean(value))
  rbind(
    dat1 %>% mutate(GeneType = "GeneBody")
    # ,dat2 %>% mutate(GeneType = "Exons")
  )
},mc.cores=length(mes_bw.rna)) %>% bind_rows()

res.rna <- as_tibble(mes_genes) %>% left_join(DE_DIVA_INFO,by = c("gene_id"="rowname"))%>% right_join(res.plot.rna,by = c("name"="rowname"))

#plot 1 compD vs nocompD
# DNAPKI=TRUE
DNAPKI=FALSE
if(DNAPKI==TRUE){

  res_comp <-
    rbind(
      res.rna %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw.DNAPKi) %>% filter(score==1)) %>% as_tibble() %>% mutate(TypeD = glue::glue("CompD")),
      res.rna %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw.DNAPKi) %>% filter(score==0)) %>% as_tibble() %>% mutate(TypeD = glue::glue("NoCompD"))
    )
}else{
  res_comp <-
    rbind(
      res.rna %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw) %>% filter(score==3)) %>% as_tibble() %>% mutate(TypeD = glue::glue("CompD")),
      res.rna %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw) %>% filter(score==0)) %>% as_tibble() %>% mutate(TypeD = glue::glue("NoCompD"))
    )
}




res_comp.duplicate <- res_comp %>% dplyr::count(gene_id,wig,GeneType) %>% filter(n>1) %>% pull(gene_id)

res_comp <- res_comp %>% filter(!gene_id %in% res_comp.duplicate)

regroup_data_per_plot <- res_comp %>% mutate(GroupPlot = case_when(str_detect(wig,"STX_vs_C1") ~ "STX_vs_C1",TRUE ~ "OHT_vs_DIVA")) %>%
  mutate(bwType = str_extract(wig,"bigwigCompare|bamCompare"))%>% mutate(wig = str_remove_all(wig,"HYL5LBGX2_|log2.bamCompare|log2.bigwigCompare")) %>%
  group_by(GroupPlot,bwType) %>% nest()

theme_set(theme_classic(base_size=12))
theme_update(axis.text.x = element_text(angle = 90))




plot_name <- glue::glue("{regroup_data_per_plot$GroupPlot[1]}_{regroup_data_per_plot$bwType[1]}")
p2 <- regroup_data_per_plot$data[[1]] %>% filter(str_detect(wig,"OHT")) %>%
  mutate(Type = case_when(Type == "Upregulated" ~ "Upregulated",TRUE ~ "NONE/DOWN")) %>% 
  ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() +facet_grid(~Type)  + scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2)) +ggtitle(plot_name)
##PVALUES
datapvals <- regroup_data_per_plot$data[[1]] %>% filter(str_detect(wig,"OHT")) %>%
  mutate(Type = case_when(Type == "Upregulated" ~ "Upregulated",TRUE ~ "NONE/DOWN"))


numbers <- datapvals %>% group_by(TypeD,Type) %>% dplyr::count() %>% spread(key=TypeD,value=n)

pvalues.dat <- datapvals %>% group_by(Type) %>% nest() %>% 
  mutate(data = map(data,function(x){
    wilcox.test(value~TypeD,data=x) %>% broom::tidy()
  })) %>% unnest(data) %>% mutate(p.value = format.pval(p.value,2))
pvalues.dat %>% left_join(numbers,by=c("Type"))

# # A tibble: 2 × 7
# # Groups:   Type [2]
# Type        statistic p.value method                                            alternative CompD NoCompD
# <chr>           <dbl> <chr>   <chr>                                             <chr>       <int>   <int>
#   1 NONE/DOWN     2052241 0.18    Wilcoxon rank sum test with continuity correction two.sided     481    8233
# 2 Upregulated     23809 0.18    Wilcoxon rank sum test with continuity correction two.sided      45     946


regroup_data_per_plot_3 <- res.rna %>% filter(str_detect(wig,"OHT")) %>%  mutate(GroupPlot = case_when(str_detect(wig,"STX_vs_C1") ~ "STX_vs_C1",TRUE ~ "OHT_vs_DIVA")) %>%
  mutate(bwType = str_extract(wig,"bigwigCompare|bamCompare"))%>% mutate(wig = str_remove_all(wig,"HYL5LBGX2_|log2.bamCompare|log2.bigwigCompare")) %>%
  mutate(Type = case_when(Type == "Upregulated" ~ "Upregulated",TRUE ~ "NONE/DOWN")) %>% 
  group_by(GroupPlot,bwType) %>% nest()

plot_name <- glue::glue("{regroup_data_per_plot_3$GroupPlot[1]}_{regroup_data_per_plot_3$bwType[1]}")
p1 <- regroup_data_per_plot_3$data[[1]] %>% ggplot(aes(x=wig,y=value,fill=Type)) + geom_boxplot()  + scale_fill_manual(values=MetBrewer::met.brewer("Egypt", 2)) +ggtitle(plot_name)

##PVALUES

datapvals <- regroup_data_per_plot_3$data[[1]] %>% filter(str_detect(wig,"OHT")) %>%
  mutate(Type = case_when(Type == "Upregulated" ~ "Upregulated",TRUE ~ "NONE/DOWN"))


numbers <- datapvals %>% group_by(Type) %>% dplyr::count() %>% spread(key=Type,value=n)
pvalues.dat <- wilcox.test(value~Type,data=datapvals) %>% broom::tidy()
pvalues.dat %>% cbind(numbers)
# statistic      p.value                                            method alternative NONE/DOWN Upregulated
# 1   2831163 8.762961e-71 Wilcoxon rank sum test with continuity correction   two.sided      8714         991
if(DNAPKI==TRUE){
  
  outfilep1p2 = "../../results/after_FC03/RNA_seq_compD_ABD_genes_090222_withDRIP_FC3_DNAPKi.pdf"
}else{
  outfilep1p2 = "../../results/after_FC03/RNA_seq_compD_ABD_genes_090222_withDRIP_FC3.pdf"
}
pdf(outfilep1p2,height=8,width=8)
print(p1)
print(p2)
dev.off()


## PLOT DISTANCE

# chr.to.study <- keep.chr <- glue::glue("chr{c(2,6,9,13,18,1,17,20,'X')}")
keep.chr <- chr.to.study 
# chr.to.study <- keep.chr <- glue::glue("chr{c(1)}")


mes_genes_dist <- list("CompD"=mes_genes.PC1 %>% as_granges(),
                  "NoCompD"=mes_genes.PC1_nocompD%>% as_granges()
)

##Distance to nearest DSB

res_distance <- mes_genes_dist %>% map(function(x){
  sres <- x %>% distanceToNearest(mes_DSB)
  x[queryHits(sres)] %>% mutate(distance = mcols(sres)$distance) %>% as_tibble()
})%>% bind_rows() %>% separate(Type,into = c("Group.PC1","Type"))


res_f <- res_distance %>% drop_na()

require(patchwork)


pb1 <- res_f %>% filter(Type == "Upregulated") %>% ggplot(aes(x=Group.PC1,y=(distance),fill=Group.PC1)) + geom_boxplot(outlier.shape=NA) + geom_jitter(alpha=0.2) +ggtitle("Upregulated Genes")
pb2 <- res_f %>% filter(Type == "Upregulated") %>% ggplot(aes(x=Group.PC1,y=log10(distance),fill=Group.PC1)) + geom_boxplot(outlier.shape=NA) + geom_jitter(alpha=0.2)
pf2 <- pb1/pb2 + plot_layout(guide="collect")
pdf("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/after_FC03/compD_value_vs_DSB_distance.pdf",height=12,width=16)
print(pf2)
dev.off()


## PLOT DISTANCE DNAPKi
mes_genes.PC1.dnapki <- genes_DIVA %>%
  filter(seqnames %in% DNAPKichr) %>%
  filter_by_overlaps(as_granges(compD.bw.DNAPKi) %>% filter(score==1)) %>%
  as_tibble() %>%
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("CompD_{Type}"))


mes_genes.PC1_nocompD.dnapki <- genes_DIVA %>%
  filter(seqnames %in% DNAPKichr) %>%
  filter_by_overlaps(as_granges(compD.bw.DNAPKi) %>% filter(score==0)) %>%
  as_tibble() %>%
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("NoCompD_{Type}")) 

mes_genes.PC1_nocompD.duplicate <- rbind(mes_genes.PC1.dnapki,mes_genes.PC1_nocompD.dnapki) %>% dplyr::count(gene_id) %>% filter(n>1) %>% pull(gene_id)

mes_genes.PC1.dnapki <- mes_genes.PC1.dnapki %>% filter(!gene_id %in% mes_genes.PC1_nocompD.duplicate)
mes_genes.PC1_nocompD.dnapki <- mes_genes.PC1_nocompD.dnapki %>% filter(!gene_id %in% mes_genes.PC1_nocompD.duplicate)


# chr.to.study <- keep.chr <- glue::glue("chr{c(2,6,9,13,18,1,17,20,'X')}")
keep.chr <- DNAPKichr
# chr.to.study <- keep.chr <- glue::glue("chr{c(1)}")


mes_genes_dist <- list("CompD"=mes_genes.PC1.dnapki %>% as_granges(),
                       "NoCompD"=mes_genes.PC1_nocompD.dnapki%>% as_granges()
)

##Distance to nearest DSB

res_distance <- mes_genes_dist %>% map(function(x){
  sres <- x %>% distanceToNearest(mes_DSB)
  x[queryHits(sres)] %>% mutate(distance = mcols(sres)$distance) %>% as_tibble()
})%>% bind_rows() %>% separate(Type,into = c("Group.PC1","Type"))


res_f <- res_distance %>% drop_na()

require(patchwork)


pb1 <- res_f %>% filter(Type == "Upregulated") %>% ggplot(aes(x=Group.PC1,y=(distance),fill=Group.PC1)) + geom_boxplot(outlier.shape=NA) + geom_jitter(alpha=0.2) +ggtitle("Upregulated Genes")
pb2 <- res_f %>% filter(Type == "Upregulated") %>% ggplot(aes(x=Group.PC1,y=log10(distance),fill=Group.PC1)) + geom_boxplot(outlier.shape=NA) + geom_jitter(alpha=0.2)
pf2 <- pb1/pb2 + plot_layout(guide="collect")
pdf("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/after_FC03/compD_value_vs_DSB_distance_DNAPKi.pdf",height=12,width=16)
print(pf2)
dev.off()


