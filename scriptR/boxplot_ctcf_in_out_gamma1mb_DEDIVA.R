# ./HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/scripts/profile_RNA_pol_on_DSB_by_active_gene.R:pdf("../results/Boxplot_CTCF_IN_OUT_gamma1Mb_DEDIVA.PDF")
require(rtracklayer)
require(tidyverse)
require(plyranges)
seqlens <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 %>% seqlengths()
bless80 = read_bed("/mnt/NAS/DATA/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
DSB174 <- read_bed("/mnt/NAS/DATA/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")
ens.genes <- read.table("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed",sep="\t",h=T) %>% GRanges()
ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))
ens.genes$name <- ens.genes$gene_name
ens.extended <- ens.genes  %>% anchor_start() %>% stretch(3000) %>% anchor_end() %>% stretch(3000)

#Gene based on ctcf position
DE_DIVA <- PhDfunc::GetDE_DIvA()
# domains <- bless80 %>% anchor_center() %>% mutate(width=2000000)
domains <- bless80 %>% anchor_center() %>% mutate(width=1000000)

direct.damaged <- ens.extended %>% filter_by_overlaps(bless80)
genes.domains <- ens.genes %>%
  filter(gene_id %in% DE_DIVA$rowname) %>%
  filter_by_overlaps(domains) %>% mutate(Cat = case_when(
    gene_id %in% direct.damaged$gene_id ~ "direct_damaged",
    TRUE ~ "gamma_domain"
  ))
genes.domains.extended.ctcf <- genes.domains  %>% anchor_start() %>% stretch(1000) %>% anchor_end() %>% stretch(1000)
ctcf_peaks <- read_narrowpeaks("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/CTCF/Clouaire_H5HJKBGXC_CTCF/PROCESSED/mapping/MACS/CTCFalfDIvA_vs_INPUT_peaks.narrowPeak")
ctcf_motifs <- read_gff("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/REVISIONS/ctcf_motifs.gff")
ctcf_motifs.in.DIVA <- ctcf_motifs %>% filter_by_overlaps(ctcf_peaks)
#ctcf peak exported for gaelle

# genes.domains.extended.ctcf <- genes.domains.extended.ctcf %>% filter_by_overlaps(ctcf_peaks) 
genes.domains.extended <- promoters(genes.domains,1000,1000)
# genes.domains.extended <- genes.domains %>% anchor_3p() %>% stretch(1000)
ctcf_motifs.in.DIVA.in.gamma <- ctcf_motifs.in.DIVA%>% filter_by_overlaps(domains)


#Split by DSB and keep only 
res <- lapply(bless80$name,function(dsbname){
  message(dsbname)
  DSB.pos <- bless80[bless80$name == dsbname]
  gamma.pos <- domains[domains$name == dsbname]
  
  sres_tot <- c()
  
  sres_ctcf <- ctcf_motifs.in.DIVA.in.gamma %>% filter_by_overlaps(gamma.pos)
  
  c(
    sres_ctcf %>% join_follow(DSB.pos) %>% mutate(Group= ifelse(strand =="+","OUT","IN")),
    sres_ctcf %>% join_precede(DSB.pos) %>% mutate(Group= ifelse(strand =="+","IN","OUT"))
  )
  
  
}) %>% do.call("c",.)


genes.domains.extended.IN <- genes.domains.extended %>% filter_by_overlaps(filter(res,Group =="IN")) %>% .$gene_id
genes.domains.extended.OUT <- genes.domains.extended %>% filter_by_overlaps(filter(res,Group =="OUT")) %>% .$gene_id

genes.domains.extended.IN.OUT <- list(
  "IN"= genes.domains[genes.domains$gene_id %in% genes.domains.extended.IN] %>% filter(Cat =="gamma_domain"),
  "OUT" = genes.domains[genes.domains$gene_id %in% genes.domains.extended.OUT]%>% filter(Cat =="gamma_domain")
)
common_in_IN <- (genes.domains.extended.IN.OUT[[1]]$gene_id %in% genes.domains.extended.IN.OUT[[2]]$gene_id)
common_in_OUT <- (genes.domains.extended.IN.OUT[[2]]$gene_id %in% genes.domains.extended.IN.OUT[[1]]$gene_id)

genes.domains.extended.IN.OUT[["IN"]] <- genes.domains.extended.IN.OUT[["IN"]][!common_in_IN]
genes.domains.extended.IN.OUT[["OUT"]] <- genes.domains.extended.IN.OUT[["OUT"]][!common_in_OUT]
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
#RNA bru seq / RNA seq
wigs_for_ctcf_genes <- c(
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
  "Bru-seq_Induced_150min"="/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/DDR_Ianelli/PROCESSED_08102018/mapping/BIGWIG/SRR5441996_normalized.bw",
  "Bru-seq_Induced_240min"="/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/DDR_Ianelli/PROCESSED_08102018/mapping/BIGWIG/SRR5441997_normalized.bw",
  "Bru-seq_Uninduced"="/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/DDR_Ianelli/PROCESSED_08102018/mapping/BIGWIG/SRR5441998_normalized.bw"
) %>% mclapply(import.bw,as="RleList",mc.cores=4)


#RNA-seq / BRUSEQ
ens.exons.domains <- ensembldb::exons(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,filter = AnnotationFilter::GeneIdFilter(genes.domains$gene_id))
ens.exons.domains <- regioneR::filterChromosomes(ens.exons.domains,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.exons.domains) <- paste0("chr",seqlevels(ens.exons.domains))
ens.exons.domains$name <- ens.exons.domains$gene_id
exons.domains.extended.IN.OUT <- list(
  "IN"= ens.exons.domains[ens.exons.domains$gene_id %in% genes.domains.extended.IN.OUT[["IN"]]$gene_id],
  "OUT" = ens.exons.domains[ens.exons.domains$gene_id %in% genes.domains.extended.IN.OUT[["OUT"]]$gene_id]
)


res.boxplot.exon.RNA <- lapply(wigs_for_ctcf_genes[c("RNA_C1_mOHT","RNA_C1_pOHT","RNA_C1_pvsmOHT","RNA_C2_mOHT","RNA_C2_pOHT","RNA_C2_pvsmOHT")],function(mon_wig){
  lapply(exons.domains.extended.IN.OUT,function(mon_bed){
    Get1valMean(Name = "",x = mon_bed,one.w = mon_wig) %>% dplyr::select(-wig)
  })%>% bind_rows(.id="Type")
  
  
}) %>% bind_rows(.id="wig")%>% mutate(Condition = str_extract(wig,"mOHT|pOHT|Induced_30min|Uninduced")) %>% mutate(wig = str_remove(wig,"_mOHT|_pOHT|_Induced_30min|_Uninduced"))
res.boxplot.exon.RNA <- res.boxplot.exon.RNA %>% group_by(wig,Type,rowname,Condition) %>% summarise(value = mean(value))



require(ggforce)
p1 <- res.boxplot.exon.RNA %>% filter(wig == "RNA_C1") %>% ggplot(aes(x=Type,y=value,fill=Condition)) + geom_boxplot()+ facet_zoom(y = value >= 5 & value <= 10) + scale_y_log10() +labs(title="RNA_C1") + theme_light() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p2 <- res.boxplot.exon.RNA %>% filter(wig == "RNA_C2") %>% ggplot(aes(x=Type,y=value,fill=Condition)) + geom_boxplot()  + facet_zoom(y = value >= 5 & value <= 10) + scale_y_log10()+labs(title="RNA_C2")+ theme_light() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
pf1 <- cowplot::plot_grid(p1,p2,ncol=1)


res.boxplot.exon.RNA %>% filter(wig == "RNA_C1")%>% group_by(Type) %>% nest() %>% 
  mutate(pval = map_dbl(data,function(x){
    x <- x %>% spread(key = Condition,value=value)
    wilcox.test(x$mOHT,x$pOHT,paired=T)$p.value
  })) %>% dplyr::select(-data)


res.boxplot.gene.RNA <- lapply(wigs_for_ctcf_genes[c("BruSeq_30vsU","BruSeq_60vsU","BruSeq_150vsU","Bru-seq_Induced_30min","Bru-seq_Induced_150min","Bru-seq_Induced_240min","Bru-seq_Uninduced")],function(mon_wig){
  lapply(genes.domains.extended.IN.OUT,function(mon_bed){
    Get1valMean(Name = "",x = mon_bed,one.w = mon_wig) %>% dplyr::select(-wig)
  })%>% bind_rows(.id="Type")
  
}) %>% bind_rows(.id="wig")%>% mutate(Condition = str_extract(wig,"mOHT|pOHT|Induced_30min|Induced_150min|Induced_240min|Uninduced")) %>% mutate(wig = str_remove(wig,"_mOHT|_pOHT|_Induced_30min|_Induced_150min|_Induced_240min|_Uninduced"))


p2.1 <- res.boxplot.gene.RNA %>% filter(wig %in% c("Bru-seq"),Condition %in% c("Induced_30min","Uninduced")) %>% ggplot(aes(x=Type,y=value,fill=Condition)) + geom_boxplot() + theme_classic() + facet_wrap(~wig) + scale_y_log10()
p2.2 <- res.boxplot.gene.RNA %>% filter(wig %in% c("Bru-seq"),Condition %in% c("Induced_150min","Uninduced")) %>% ggplot(aes(x=Type,y=value,fill=Condition)) + geom_boxplot() + theme_classic() + facet_wrap(~wig) + scale_y_log10()
p2.3 <- res.boxplot.gene.RNA %>% filter(wig %in% c("Bru-seq"),Condition %in% c("Induced_240min","Uninduced")) %>% ggplot(aes(x=Type,y=value,fill=Condition)) + geom_boxplot() + theme_classic() + facet_wrap(~wig) + scale_y_log10()


res.boxplot.gene.RNA %>% filter(wig %in% c("Bru-seq"),Condition %in% c("Induced_30min","Uninduced"))  %>% group_by(Type) %>% nest() %>% 
  mutate(pval = map_dbl(data,function(x){
    x <- x %>% spread(key = Condition,value=value)
    wilcox.test(x$Induced_30min,x$Uninduced,paired=T)$p.value
    })) %>% dplyr::select(-data)
