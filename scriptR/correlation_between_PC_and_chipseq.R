require(tidyverse)
require(ranger)
require(corrr)
require(plyranges)
require(rtracklayer)
library(FactoMineR)
# Correlation between PCs and Histone marks
# PCs.ratio.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHTDNAPKi.bw"
PCs.ratio.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipA.bw"
# PCs.ratio.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHT.bw"
# PCs.ratio.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio.bw"
# PCs.ratio.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHT.bw"
PCs.ratio <- import.bw(PCs.ratio.path)
out_name_file <- basename(PCs.ratio.path) %>% str_remove(".bw")
PCs.ratio$name <- str_c("bin",1:length(PCs.ratio))
list_bw <- c(
  "RNA_C1_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYL5LBGX2/HYL5LBGX2_RNA_C1_1_17s002464-1-1_Clouaire_lane117s002464_hg19plusERCC_normalized.bw",
  "RNA_C1_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYL5LBGX2/HYL5LBGX2_RNA_C1_OHT_1_17s002465-1-1_Clouaire_lane117s002465_hg19plusERCC_normalized.bw",
  # "RNA_C2_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYMHNBGX2/HYMHNBGX2_RNA_C1_2_17s002468-1-1_Clouaire_lane117s002468_hg19plusERCC_normalized.bw",
  # "RNA_C2_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYMHNBGX2/HYMHNBGX2_RNA_C1_OHT_2_17s002469-1-1_Clouaire_lane117s002469_hg19plusERCC_normalized.bw",
  # "RNA_C1_pvsmOHT"="/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/BAMCOMPARE/HYL5LBGX2_RNA_C1_bamcompare_logFC_normalized.bw",
  # "RNA_C2_pvsmOHT"="/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/BAMCOMPARE/HYMHNBGX2_RNA_C2_bamcompare_logFC_normalized.bw",
  # "Bru-seq_Induced_30min"="/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/DDR_Ianelli/PROCESSED_08102018/mapping/BIGWIG/SRR5441994_normalized.bw",
  # "Bru-seq_Uninduced"="/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/DDR_Ianelli/PROCESSED_08102018/mapping/BIGWIG/SRR5441998_normalized.bw",
  # "RNA_siSCC1_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/RNASeq_Coline/PROCESSED/mapping/BIGWIG/HVHWHBGXG_RNAseqA_siSCC1_DIvA_20s003791-1-1_Clouaire_lane120s003791_normalized.bw",
  # "RNA_siSCC1_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/RNASeq_Coline/PROCESSED/mapping/BIGWIG/HVHWHBGXG_RNAseqA_siSCC1_OHT_20s003794-1-1_Clouaire_lane120s003794_normalized.bw",
  # "RNA_siWAPL_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/RNASeq_Coline/PROCESSED/mapping/BIGWIG/HVHWHBGXG_RNAseqA_siWAPL_DIvA_20s003792-1-1_Clouaire_lane120s003792_normalized.bw",
  # "RNA_siWAPL_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/RNASeq_Coline/PROCESSED/mapping/BIGWIG/HVHWHBGXG_RNAseqA_siWAPL_OHT_20s003795-1-1_Clouaire_lane120s003795_normalized.bw",
  "h3k79me2_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K79me2/BGI_RUN1_201203/PROCESSED/BIGWIG/h3k79me2_mOHT_normalized_hg19.bw",
  "h3k79me2_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K79me2/BGI_RUN1_201203/PROCESSED/BIGWIG/h3k79me2_pOHT_normalized_hg19.bw",
  "H3K4me3_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K4me3/BGI_RUN3_201412/PROCESSED/BIGWIG/H3K4me3_normalized_hg19.bw",
  "H3K4me3_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K4me3/BGI_RUN3_201412/PROCESSED/BIGWIG/H3K4me3_OHT_normalized_hg19.bw",
  "h2bub_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H2BK120Ub/BGI_RUN1_201203/PROCESSED/BIGWIG/h2bub_mOHT_normalized_hg19.bw",
  "h2bub_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H2BK120Ub/BGI_RUN1_201203/PROCESSED/BIGWIG/h2bub_pOHT_normalized_hg19.bw",
  "H3K9me3_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K9me3/EMBL_H2VFWBGX5_HWNTLBGX3/PROCESSED/mapping/BIGWIG/HWNTLBGX3_H3K9me3_DIVA_17s005723-1-1_Clouaire_lane117s005723_sequence_normalized.bw",
  "H3K9me3_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K9me3/EMBL_H2VFWBGX5_HWNTLBGX3/PROCESSED/mapping/BIGWIG/HWNTLBGX3_H3K9me3_OHT_DIVA_17s005724-1-1_Clouaire_lane117s005724_sequence_normalized.bw",
  "H4K20me3_AS_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H4K20me3/EMBL4_HYJW5BGX5_AS_S/PROCESSED/mapping/BIGWIG/HYJW5BGX5_H4K20me3_AS_DIVA_18s000866-1-1_Clouaire_lane118s000866_sequence_normalized.bw",
  "H4K20me3_AS_pOHT"= "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H4K20me3/EMBL4_HYJW5BGX5_AS_S/PROCESSED/mapping/BIGWIG/HYJW5BGX5_H4K20me3_AS_OHT_18s000867-1-1_Clouaire_lane118s000867_sequence_normalized.bw",
  "H4K20me3_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H4K20me3/EMBL1_HMNY7BGXX_AS/PROCESSED/BIGWIG/HMNY7BGXX_H4K20me3_DIVA_16s000422-1-1_Clouaire_lane116s000422_sequence.nodups.bam_normalized.bw",
  "H4K20me3_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H4K20me3/EMBL1_HMNY7BGXX_AS/PROCESSED/BIGWIG/HMNY7BGXX_H4K20me3_OHT_DI_16s000423-1-1_Clouaire_lane116s000423_sequence.nodups.bam_normalized.bw",
  "h4k20me1_mOHT"= "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H4K20me1/BGI_RUN1_201203/PROCESSED/BIGWIG/h4k20me1_mOHT_normalized_hg19.bw",
  "h4k20me1_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H4K20me1/BGI_RUN1_201203/PROCESSED/BIGWIG/h4k20me1_pOHT_normalized_hg19.bw",
  "H4K20me1Mono_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H4K20me1/BGI_RUN3_201412/PROCESSED/BIGWIG/H4K20me1Mono_normalized_hg19.bw",
  "H4K20me1Mono_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H4K20me1/BGI_RUN3_201412/PROCESSED/BIGWIG/H4K20me1Mono_OHT_normalized_hg19.bw",
  "H1_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H1/PROCESSED/BIGWIG/HN25TBGXX_H1_DIVA_16s000414-1-1_Clouaire_lane116s000414_sequence.nodups.bam_normalized.bw",
  "H1_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H1/PROCESSED/BIGWIG/HN25TBGXX_H1_OHT_DIVA_16s000415-1-1_Clouaire_lane116s000415_sequence.nodups.bam_normalized.bw",
  "H2AZ_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H2AZ/EMBL_HMNY7BGXX/PROCESSED/BIGWIG/HMNY7BGXX_H2AZ_DIVA_16s000418-1-1_Clouaire_lane116s000418_sequence.nodups.bam_normalized.bw",
  "H2AZ_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H2AZ/EMBL_HMNY7BGXX/PROCESSED/BIGWIG/HMNY7BGXX_H2AZ_OHT_DIVA_16s000419-1-1_Clouaire_lane116s000419_sequence.nodups.bam_normalized.bw",
  "h3k36me3_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K36me3/BGI_RUN1_201203/PROCESSED/BIGWIG/h3k36me3_mOHT_normalized_hg19.bw",
  "h3k36me3_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K36me3/BGI_RUN1_201203/PROCESSED/BIGWIG/h3k36me3_pOHT_normalized_hg19.bw",
  "H2AZac_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H2AZac/EMBL_C8LVDACXX/PROCESSED/BIGWIG/H2AZac_DIVA_normalized.bw",
  "H2AZac_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H2AZac/EMBL_C8LVDACXX/PROCESSED/BIGWIG/H2AZac_OHT_DIVA_normalized.bw"
  # "repli-Seq_Rep_12_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Repli-seq/REPLISEQ_Clouaire_H2VJLBGX5_HC2FTBGX9_HWWTVBGX7/PROCESSED_repliseq_26_10_18/mapping/BIGWIG/H2VJLBGX5/H2VJLBGX5_Rep_12_17s005733-1-1_Clouaire_lane117s005733_normalized.bw",
  # "repli-Seq_Rep_10_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Repli-seq/REPLISEQ_Clouaire_H2VJLBGX5_HC2FTBGX9_HWWTVBGX7/PROCESSED_repliseq_26_10_18/mapping/BIGWIG/H2VJLBGX5/H2VJLBGX5_Rep_10_17s005732-1-1_Clouaire_lane117s005732_normalized.bw",
  # "53BP1_S_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/53BP1/EMBL3_HV5MYBGX5_OHT_NT_AS_S/Clouaire_HV5MYBGX5/PROCESSED/mapping/BIGWIG/HV5MYBGX5_53BP1_S_DIVA_18s000856-1-1_Clouaire_lane118s000856_sequence_normalized.bw",
  # "pATM_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/phosphoATM/Clouaire_HLNKYBGXC_pATM/PROCESSED/BIGWIG/Clouaire_HLNKYBGXC_pATMdiva.bw",
  # "GAM_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/GAMMA/U2OS/BGI_RUN1_201203/PROCESSED/ALIGNED/WIGGLE/hg19/GAM.clean_hg19_normalized.bw"
  # "H4K20me2_AS_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H4K20me2/EMBL_HVL7VBGX5_AS/PROCESSED/mapping/BIGWIG/HVL7VBGX5_H4K20me2_AS_DIVA_18s000862-1-1_Clouaire_lane118s000862_sequence_normalized.bw",
  # "H4K20me2_AS_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H4K20me2/EMBL_HVL7VBGX5_AS/PROCESSED/mapping/BIGWIG/HVL7VBGX5_H4K20me2_AS_OHT_18s000863-1-1_Clouaire_lane118s000863_sequence_normalized.bw"
  # "H3K9me3_pvsmOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K9me3/EMBL_H2VFWBGX5_HWNTLBGX3/PROCESSED/mapping/BIGWIG/BamCompare_log2ratio_pOHTvsmOHt_bin50/H3K9me3_EMBL_bin50bp_log2ratio.bw",
  # "h3k79me2_pvsmOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K79me2/BGI_RUN1_201203/PROCESSED/BIGWIG/h3k79me2_BGI_RUN_1_bin50_log2ratio.bw",
  # "H3K4me3_pvsmOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K4me3/BGI_RUN3_201412/PROCESSED/BIGWIG/H3K4me3_BGI_RUN_3_bin50_log2ratio.bw",
  # "h2bub_pvsmOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H2BK120Ub/BGI_RUN1_201203/PROCESSED/BIGWIG/h2bub_BGI_RUN_1_bin50_log2ratio.bw"
)
sub_list_bw <- list_bw[str_detect(names(list_bw),"_pOHT")]
# list_bw <- append(sub_list_bw,list_bw[c("repli-Seq_Rep_12_mOHT","repli-Seq_Rep_10_mOHT","53BP1_S_mOHT","pATM_mOHT")])

# DSB174 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")
mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
  arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 1000000)

if(str_detect(basename(PCs.ratio.path),"OHTDNAPKi")){
  chr.to.study <- glue::glue("chr{c(2,6,9,13,18,1,17,20,'X')}")
}else{
  chr.to.study <- glue::glue("chr{c(1,17,'X')}")
}


PCs.ratio.nogamma <- PCs.ratio %>% filter_by_non_overlaps(gamma_region) %>% filter(seqnames %in% chr.to.study)

full_signal <- lapply(sub_list_bw,function(one_w){
  bw.data <- import.bw(one_w,as="RleList")
  message(one_w)
  res <- PhDfunc::Get1val(Name = "binsPC1",x = PCs.ratio.nogamma,one.w = bw.data)
  
}) %>% bind_rows(.id="wig")

full_signal_tsv <- full_signal %>%
  left_join(as_tibble(PCs.ratio.nogamma) %>% dplyr::select(seqnames,score,name),by = c("rowname"="name")) %>% 
  spread(key = wig,value = value)

colnames(full_signal_tsv) <- str_replace_all(colnames(full_signal_tsv),"-","_")

# my_corr <- full_signal_tsv %>% filter(seqnames %in% c("chr1","chr17"))%>% dplyr::select(-rowname,-seqnames) %>% correlate() %>%  focus(score)
# 
# my_corr %>%
#     ggplot(aes(x=fct_reorder(rowname,score),y=score,fill=score)) +
#     geom_bar(stat="identity") +
#     coord_flip() + theme_classic(base_size=18) +
#     scale_fill_distiller(palette = "Spectral") -> p_corr
# pdf("results/corr_between_pCs_bins_and_histones_no_DSB_chr1_chr17.pdf",height=12,width=8)
# print(p_corr)
# dev.off()
# 

my_corr <- full_signal_tsv %>% dplyr::select(-rowname,-seqnames) %>% correlate() %>%  focus(score)

my_corr %>%
  ggplot(aes(x=fct_reorder(rowname,score),y=score,fill=score)) +
  geom_bar(stat="identity") +
  coord_flip() + theme_classic(base_size=18) +
  scale_fill_distiller(palette = "Spectral") -> p_corr


my_corr <- full_signal_tsv %>% dplyr::select(-rowname,-seqnames) %>% correlate(method = "spearman") %>%  focus(score)

my_corr %>%
  ggplot(aes(x=fct_reorder(rowname,score),y=score,fill=score)) +
  geom_bar(stat="identity") +
  coord_flip() + theme_classic(base_size=18) +
  scale_fill_distiller(palette = "Spectral") -> p_corr_spearman


data_plot <- full_signal_tsv %>% gather(key = Histone,value = value,-rowname:-score)

myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
p <- data_plot %>% ggplot(aes(x=score,y=value)) + geom_bin2d(bins=100)  + facet_wrap(~Histone,scales="free_y",ncol=3) + scale_fill_gradientn(colours = myColor) + theme_minimal(base_size = 8)  + scale_y_log10() 


pdf(glue::glue("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/corr_between_pCs_bins_and_histones_no_DSB_spearman_200_{out_name_file}.pdf"),height=12,width=12)
print(p_corr + ggtitle("pearson"))
print(p_corr_spearman + ggtitle("spearman"))
print(p)
dev.off()


#Correlation between multiple compD and gamma
mes_bw <- c(
  list.files("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/",pattern="PC1_all_chr_log2ratio_100kb_.+\\.bw",full.names =T),
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/GAMMA/U2OS/BGI_RUN1_201203/PROCESSED/ALIGNED/WIGGLE/hg19/GAM.clean_hg19_normalized.bw"
) %>% setNames(basename(.))%>% map(import.bw,as="RleList")



bless80 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
my_dsb <- regioneR::filterChromosomes(bless80,keep.chr = keep.chr)


my_dsb.chr<- my_dsb %>%
  anchor_center() %>% mutate(width = 1000000) 

Random.pos <- regioneR::randomizeRegions(DSB174.chr,per.chromosome = T)
Random.pos$name <- str_c("Random",1:length(Random.pos))

Random80 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/80random.bed")  %>%
  anchor_center() %>% mutate(width = 1000000) %>% regioneR::filterChromosomes(keep.chr = keep.chr)


full_signal <- lapply(names(mes_bw),function(one_n){
  list("DSB"=my_dsb.chr,"Ctrl"=Random.pos,"Random80"=Random80) %>% lapply(function(one_bed){
    PhDfunc::Get1val(Name = one_n,x = one_bed,one.w = mes_bw[[one_n]])
  }) %>% bind_rows(.id="Type")
}) %>% bind_rows()


full_signal.corr <- full_signal %>%
  spread(key=wig,value=value) %>% 
  group_by(Type) %>% nest() %>% 
  mutate(data = map(data,select,-rowname)) %>%
  mutate(data = map(data,correlate)) %>%
  mutate(data = map(data,focus,GAM.clean_hg19_normalized.bw)) %>% unnest(data)

p_corr_gamma <- full_signal.corr %>%
  mutate(rowname = str_remove_all(rowname,"PC1_all_chr_log2ratio_100kb_|.bw")) %>% 
  filter(rowname != "HiC_D_OHTPARPi") %>% 
  mutate(rowname = fct_recode(rowname,`Replicate 1`="OHT_manipA",`Replicate 2`="OHT_manipB",`Replicate 3`="HiC_D_OHT",`DNAPKi`="HiC_D_OHTDNAPKi")) %>% 
  mutate(rowname = factor(rowname,labels = rev(c("Replicate 1","Replicate 2","Replicate 3","DNAPKi")))) %>% 
ggplot(aes(y=rowname,x=Type,fill=GAM.clean_hg19_normalized.bw,label = round(GAM.clean_hg19_normalized.bw,1))) + geom_tile() +
  geom_text() +
  theme_classic() + scale_fill_gradient2(high = "#c0392b",low = "#2980b9",mid = "white",limits = c(-1,1)) +
  ylab("") + xlab("")
pdf("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/correlation_gamma_withDSB_control_1_17_X.pdf",height=6)
print(p_corr_gamma)
dev.off()
