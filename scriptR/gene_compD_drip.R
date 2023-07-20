Get1valMean <- function (Name, one.w, x) 
  ##First plot -> 3 ABD with 1,17,X
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

# %>% filter(score == 3)

# 
mes_genes.PC1 <- genes_DIVA %>%
  # filter(seqnames %in% chr.to.study) %>%
  filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 3)) %>%
  as_tibble() %>%
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(TypeD = glue::glue("CompD"))

# mes_genes.PC1 %>% dplyr::select(seqnames,start,end,strand,gene_name,Type,TypeD) %>% write_tsv("../../results/Genes_CompD.bed",col_names = F)
# mes_genes.PC1 %>% dplyr::select(seqnames,start,end,strand,gene_name,Type,TypeD) %>% filter(Type=="Upregulated") %>% write_tsv("../../results/Genes_CompD_Upregulated_0.3.bed",col_names = F)

mes_genes.PC1_nocompD <- genes_DIVA %>%
  filter(seqnames %in% chr.to.study) %>%
  filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 0)) %>%
  as_tibble() %>%
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(TypeD = glue::glue("NoCompD"))
# 
# mes_genes.PC1_nocompD %>% dplyr::select(seqnames,start,end,strand,gene_name,Type,TypeD) %>% write_tsv("../../results/Genes_NoCompD.bed",col_names = F)
# mes_genes.PC1_nocompD %>% dplyr::select(seqnames,start,end,strand,gene_name,Type,TypeD) %>% filter(Type=="Upregulated") %>% write_tsv("../../results/Genes_NoCompD_Upregulated_0.3.bed",col_names = F)


# mes_genes <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD) %>% as_granges()
mes_genes <- genes_DIVA
mes_genes$name <- mes_genes$gene_id


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

res1 <- 
  rbind(
    res %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 3)) %>% as_tibble() %>% mutate(TypeD = glue::glue("CompD")),
    res %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 0)) %>% as_tibble() %>% mutate(TypeD = glue::glue("NoCompD"))
  )

res1.duplicate <- res1 %>% dplyr::count(gene_id,wig,GeneType) %>% filter(n>1) %>% pull(gene_id)

res1 <- res1 %>% filter(!gene_id %in% res1.duplicate)

p1 <- res1 %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() + scale_y_log10() + facet_wrap(~GeneType)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
# p2 <- res1 %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() + scale_y_log10() + facet_wrap(GeneType~Type) + theme_classic()+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
### same but with chr1 only

res2 <- res1 %>% filter(seqnames =="chr1")
p3 <- res2 %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() + scale_y_log10()  + facet_wrap(~GeneType)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
# p4 <- res2 %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() + scale_y_log10() + facet_wrap(GeneType~Type) + theme_classic()+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))




###SAME with DNAPKi


res3 <- 
  rbind(
    res %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw.DNAPKi)%>% filter(score == 1)) %>% as_tibble() %>% mutate(TypeD = glue::glue("CompD")),
    res %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw.DNAPKi)%>% filter(score == 0)) %>% as_tibble() %>% mutate(TypeD = glue::glue("NoCompD"))
  )

res3.duplicate <- res3 %>% dplyr::count(gene_id,wig,GeneType) %>% filter(n>1) %>% pull(gene_id) %>% unique()

res3 <- res3 %>% filter(!gene_id %in% res3.duplicate)

p5 <- res3 %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() + scale_y_log10()  + facet_wrap(~GeneType)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
# p2 <- res1 %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() + scale_y_log10() + facet_wrap(GeneType~Type) + theme_classic()+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))

##PLOTS

pdf("../../results/Drip_seq_compD_ABD_genes.pdf",height=10,width=18)
print(p1)
# print(p2)
print(p3 + ggtitle("chr1"))
# print(p4 + ggtitle("chr1"))

print(p5 + ggtitle(glue::glue("DNAPKi ({paste(DNAPKichr,collapse=',')})")))
dev.off()
#plot au 14/10 en genebody only (et en ayant viré les qDRIP)

res4 <- 
  rbind(
    res %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 3)) %>% as_tibble() %>% mutate(TypeD = glue::glue("CompD")),
    res %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 0)) %>% as_tibble() %>% mutate(TypeD = glue::glue("NoCompD"))
  ) %>% filter(GeneType == "GeneBody")


res4.duplicate <- res4 %>% dplyr::count(gene_id,wig,GeneType) %>% filter(n>1) %>% pull(gene_id) %>% unique()

res4 <- res4 %>% filter(!gene_id %in% res4.duplicate)

p7 <- res4 %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() + scale_y_log10()  + facet_wrap(~GeneType)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
p7.2 <- res4 %>% filter(seqnames =="chr1")%>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() + scale_y_log10()  + facet_wrap(~GeneType)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))

pdf("../../results/Drip_seq_compD_ABD_genes_03022022.pdf",height=6,width=5)
p702 <- res4 %>% filter(wig =="DRIP2_C1_pOHT") %>%  ggplot(aes(x=wig,y=value,fill=TypeD)) +
  geom_boxplot() + scale_y_log10()+coord_cartesian(ylim = c(0.1,2.5))  + facet_wrap(~GeneType)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
print(p702)
dev.off()

##PVALUES
numbers <- res4 %>% dplyr::count(wig,TypeD) %>% spread(key = TypeD,value=n)
pvalues.dat <- res4 %>% group_by(wig) %>% nest() %>% 
  mutate(data = map(data,function(x){
    wilcox.test(value~TypeD,data=x) %>% broom::tidy()
  })) %>% unnest(data) %>% mutate(p.value = format.pval(p.value,3))
pvalues.dat %>% left_join(numbers,by="wig")
#profil gene body en qDRIP de ces gènes

gene_bed <- split(res4,res4$TypeD) %>% map(as_granges) %>% map(unique)

for(wig_n in names(mes_bw.list)){
  # wig_n <- basename(wig) %>% str_remove("_normalized.bw")
  wig <- mes_bw.list[[wig_n]]
  message(wig_n)
  for(my_bed in names(gene_bed)){
    
    bed_n <- glue::glue("gene_in_{my_bed}.bed")
    export.bed(gene_bed[[my_bed]],bed_n)
    message(bed_n)
    outfile <- str_c(my_bed,wig_n,"tab.gz",sep=".")
    outfile <- str_c("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/qDRIP_profil_genes/MATRIX_GENES/",outfile)
    if(!file.exists(outfile)){
      system(str_c("/home/rochevin/.local/bin/computeMatrix scale-regions -S ",wig," -R ",bed_n," -bs 50 -m 10000 -b 3000 -a 3000 -out ",outfile," -p 8 --skipZeros --missingDataAsZero"))
    }
    file.remove(bed_n)
  }
}



mydir <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/qDRIP_profil_genes/MATRIX_GENES/"
files <- list.files(mydir)
dataplot.drip <- mclapply(files,function(onefile){
  fread(cmd=str_c("zcat ",mydir,onefile),skip=1,header=F)  %>% dplyr::select(-V1:-V6) %>% colMeans()  %>%
    enframe() %>% mutate(seq = 1:dplyr::n())%>%
    mutate(wig = str_remove(onefile,".tab.gz")) %>% separate(wig,into = c("Bed","File"),sep="\\.")
},mc.cores=12) %>% bind_rows()

# dataplot.ctcf <- mclapply(files,function(onefile){
#   fread(cmd=str_c("zcat ",mydir,onefile),skip=1,header=F)  %>% dplyr::select(-V1:-V6) %>% colMeans()  %>%
#     enframe() %>% mutate(seq = seq(-3000 - ,3000 +5,5))%>%
#     mutate(wig = str_remove(onefile,".tab.gz")) %>% separate(wig,into = c("Bed","File"),sep="\\.")
# },mc.cores=12) %>% bind_rows()


profile.drip <- dataplot.drip %>% 
  ggplot(aes(seq,value,col=Bed)) +
  geom_line()+
  geom_vline(xintercept = 60,linetype="longdash") +
  geom_vline(xintercept = 260,linetype="longdash") +
  scale_x_continuous(name = 'Position',
                     breaks = c(1,60,260,320),
                     labels = c("TSS-3kb", 'TSS', 'TES', 'TES+3kb')
  ) +
  facet_grid(File~.) +
  theme_classic(base_size = 14) +
  theme(legend.position = "bottom") + scale_color_manual(values=MetBrewer::met.brewer("Thomas", 2))


pdf("../../results/Drip_seq_compD_ABD_genes_14012022.pdf",height=6,width=5)
print(p7)
print(p7.2 + ggtitle("chr1"))
# print(p4 + ggtitle("chr1"))
print(profile.drip)
dev.off()

##AU 20/01/22 
## ADD BIGWIG COMPARE TO THE PLOT (GENE BODY)
##with stranded data
## if gene == +, then plot reverse qDRIP

mes_bwc.list <- c(
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/qDRIP/Clouaire_HWGF7BGXJ_4h_24h_OHT_ALINE/PROCESSED/mapping/BIGWIG/BIGWIG_COMPARE" %>%
    list.files(pattern="_spikeinfactor_subtract_bigwigcompare.bw",full.names=T) %>% str_subset("total",negate=T)%>% str_subset("qDRIPOHT_vs_DIVA_"),
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/qDRIP/Clouaire_HH7TKBGXK_SETX_FTO_FLORIAN/Clouaire_HH7TKBGXK_SETX_FTO/PROCESSED/mapping/BIGWIG/BIGWIG_COMPARE" %>%
    list.files(pattern="_spikeinfactor_subtract_bigwigcompare.bw",full.names=T) %>% str_subset("total",negate=T) %>% str_subset("siFTODIVA_[total|reverse|forward]",negate=T)
) 
mes_bwc <- mes_bwc.list %>%mclapply(import.bw,as="RleList",mc.cores=length(mes_bwc.list))
names(mes_bwc) <- basename(mes_bwc.list)
res.plot.compare <- mclapply(names(mes_bwc),function(one_wig){
  # my_x <- mes_genes %>% promoters(500,500)
  my_x <- mes_genes %>% anchor_start() %>% plyranges::stretch(1000)%>% anchor_end() %>% plyranges::stretch(1000)
  dat1 <- Get1valMean(Name = one_wig,one.w = mes_bwc[[one_wig]],x = my_x)
  
  dat1
},mc.cores=length(mes_bwc)) %>% bind_rows()


res.compare <- as_tibble(mes_genes) %>% left_join(DE_DIVA_INFO,by = c("gene_id"="rowname"))%>% right_join(res.plot.compare,by = c("name"="rowname"))


res5 <- 
  rbind(
    res.compare %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 3)) %>% as_tibble() %>% mutate(TypeD = glue::glue("CompD")),
    res.compare %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 0)) %>% as_tibble() %>% mutate(TypeD = glue::glue("NoCompD"))
  )


res5.duplicate <- res5 %>% dplyr::count(gene_id,wig) %>% filter(n>1) %>% pull(gene_id) %>% unique()

res5 <- res5 %>% filter(!gene_id %in% res5.duplicate)

res5 <- res5 %>%
  mutate(Orientation = str_remove(str_extract(wig,"forward_.+|total_.+|reverse_.+"),"_spikeinfactor_subtract_bigwigcompare.bw")) %>% 
  mutate(wig = str_remove(wig,"_(forward_|total_|reverse_).+"))%>% filter(Orientation == "forward" & strand =="-" | Orientation == "reverse" & strand =="+")


p8 <- res5 %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot()+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
p8.1 <- p8 + coord_cartesian(ylim = c(-2,2)) 
p8.2 <- res5 %>% filter(seqnames =="chr1") %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot()+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2)) +ggtitle("chr1")
p8.3 <- p8.2 + coord_cartesian(ylim = c(-2,2))
##Same thing but with regular bigwig


mes_bw.list2 <- c(
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/qDRIP/Clouaire_HWGF7BGXJ_4h_24h_OHT_ALINE/PROCESSED/mapping/BIGWIG" %>%
    list.files(pattern="_spikeinfactor.bw",full.names=T) %>% str_subset("total",negate=T) %>% str_subset("qDRIPDIVA_|qDRIPOHT_"),
  "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/qDRIP/Clouaire_HH7TKBGXK_SETX_FTO_FLORIAN/Clouaire_HH7TKBGXK_SETX_FTO/PROCESSED/mapping/BIGWIG" %>%
    list.files(pattern="_spikeinfactor.bw",full.names=T) %>% str_subset("total",negate=T) %>% str_subset("siFTOOHT_|siFTODIVA_",negate=T)
) 
mes_bw2 <- mes_bw.list2 %>%mclapply(import.bw,as="RleList",mc.cores=length(mes_bw.list2))
names(mes_bw2) <- basename(mes_bw.list2)
res.plot.2 <- mclapply(names(mes_bw2),function(one_wig){
  # my_x <- mes_genes %>% promoters(500,500)
  my_x <- mes_genes %>% anchor_start() %>% plyranges::stretch(1000)%>% anchor_end() %>% plyranges::stretch(1000)
  dat1 <- Get1valMean(Name = one_wig,one.w = mes_bw2[[one_wig]],x = my_x)
  
  dat1
},mc.cores=length(mes_bw2)) %>% bind_rows()


res.2 <- as_tibble(mes_genes) %>% left_join(DE_DIVA_INFO,by = c("gene_id"="rowname"))%>% right_join(res.plot.2,by = c("name"="rowname"))


res6 <- 
  rbind(
    res.2 %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 3)) %>% as_tibble() %>% mutate(TypeD = glue::glue("CompD")),
    res.2 %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 0)) %>% as_tibble() %>% mutate(TypeD = glue::glue("NoCompD"))
  )


res6.duplicate <- res6 %>% dplyr::count(gene_id,wig) %>% filter(n>1) %>% pull(gene_id) %>% unique()

res6 <- res6 %>% filter(!gene_id %in% res6.duplicate)

res6 <- res6 %>%
  mutate(Orientation = str_remove(str_extract(wig,"forward_.+|total_.+|reverse_.+"),"_spikeinfactor.bw")) %>% 
  mutate(wig = str_remove_all(wig,".+_lane1|_(forward_|total_|reverse_).+"))%>%
  mutate(Condition = str_extract(wig,"OHT|DIVA")) %>% 
  mutate(wig = str_remove(wig,"OHT|DIVA")) %>% 
  filter(Orientation == "forward" & strand =="-" | Orientation == "reverse" & strand =="+")

res6 <- res6 %>% spread(key=Condition,value=value) %>% mutate(diff = OHT-DIVA,ratio=log2(OHT/DIVA))

p9 <- res6 %>% gather(key=CalculationType,value=value,-seqnames:-OHT) %>%  ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot()+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2)) +facet_grid(~CalculationType)
p9.1 <- p9 + coord_cartesian(ylim = c(-2,2)) 
p9.2 <- res6 %>% gather(key=CalculationType,value=value,-seqnames:-OHT)  %>% filter(seqnames =="chr1") %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot()+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2)) +ggtitle("chr1")+facet_grid(~CalculationType)
p9.3 <- p9.2 + coord_cartesian(ylim = c(-2,2))
pdf("../../results/Drip_seq_compD_ABD_genes_20012022.pdf",height=6,width=5)

print(p8.1)
print(p8.3)
print(p9.1)
print(p9.3)
dev.off()



#### PLOT CompD/noCOmpD with RNA-seq


#script au 090222
# mettre a coté compD non compD sur le boxplot 
# 
# le faire sur tous les genes sans compD non compD (pour les gènes up/down/none)
# et faire a coté tous les genes vs up 
# 
# 3eme graphe page1 merger down et none, et virer les diva (only OHT) 
# faut mettre à coter les compD/non compD dans les plots
# virer les plots chromosomes 1 
# virer oht vs DIVA en SETX et en CTRL


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
res_comp <-
  rbind(
    res.rna %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 3)) %>% as_tibble() %>% mutate(TypeD = glue::glue("CompD")),
    res.rna %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 0)) %>% as_tibble() %>% mutate(TypeD = glue::glue("NoCompD"))
  )

res_comp.duplicate <- res_comp %>% dplyr::count(gene_id,wig,GeneType) %>% filter(n>1) %>% pull(gene_id)

res_comp <- res_comp %>% filter(!gene_id %in% res_comp.duplicate)

regroup_data_per_plot <- res_comp %>% mutate(GroupPlot = case_when(str_detect(wig,"STX_vs_C1") ~ "STX_vs_C1",TRUE ~ "OHT_vs_DIVA")) %>%
  mutate(bwType = str_extract(wig,"bigwigCompare|bamCompare"))%>% mutate(wig = str_remove_all(wig,"HYL5LBGX2_|log2.bamCompare|log2.bigwigCompare")) %>%
  group_by(GroupPlot,bwType) %>% nest()

theme_set(theme_classic(base_size=12))
theme_update(axis.text.x = element_text(angle = 90))
res_plot_1 <- map(1:nrow(regroup_data_per_plot),function(i){
  plot_name <- glue::glue("{regroup_data_per_plot$GroupPlot[i]}_{regroup_data_per_plot$bwType[i]}")
  p1 <- regroup_data_per_plot$data[[i]] %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot()  + scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2)) +ggtitle(plot_name) +
    facet_zoom(y = value >= -0.25 & value <= 0.25) 
  p1
})
#plot2 UP/DOWN/NONE ALL GENES (all chromosomes)
require(ggforce)
regroup_data_per_plot_2 <- res.rna %>% filter(str_detect(wig,"OHT")) %>%  mutate(GroupPlot = case_when(str_detect(wig,"STX_vs_C1") ~ "STX_vs_C1",TRUE ~ "OHT_vs_DIVA")) %>%
  mutate(bwType = str_extract(wig,"bigwigCompare|bamCompare"))%>% mutate(wig = str_remove_all(wig,"HYL5LBGX2_|log2.bamCompare|log2.bigwigCompare")) %>%
  group_by(GroupPlot,bwType) %>% nest()

res_plot_2 <- map(1:nrow(regroup_data_per_plot_2),function(i){
  plot_name <- glue::glue("{regroup_data_per_plot_2$GroupPlot[i]}_{regroup_data_per_plot_2$bwType[i]}")
  p1 <- regroup_data_per_plot_2$data[[i]] %>% ggplot(aes(x=wig,y=value,fill=Type)) + geom_boxplot() +
    facet_zoom(y = value >= -0.25 & value <= 0.25) +ggtitle(plot_name)
  
})
#merge NONE/DOWN
regroup_data_per_plot_3 <- res.rna %>% filter(str_detect(wig,"OHT")) %>%  mutate(GroupPlot = case_when(str_detect(wig,"STX_vs_C1") ~ "STX_vs_C1",TRUE ~ "OHT_vs_DIVA")) %>%
  mutate(bwType = str_extract(wig,"bigwigCompare|bamCompare"))%>% mutate(wig = str_remove_all(wig,"HYL5LBGX2_|log2.bamCompare|log2.bigwigCompare")) %>%
  mutate(Type = case_when(Type == "Upregulated" ~ "Upregulated",TRUE ~ "NONE/DOWN")) %>% 
  group_by(GroupPlot,bwType) %>% nest()
res_plot_3 <- map(1:nrow(regroup_data_per_plot_2),function(i){
  plot_name <- glue::glue("{regroup_data_per_plot_3$GroupPlot[i]}_{regroup_data_per_plot_3$bwType[i]}")
  p1 <- regroup_data_per_plot_3$data[[i]] %>% ggplot(aes(x=wig,y=value,fill=Type)) + geom_boxplot()  + scale_fill_manual(values=MetBrewer::met.brewer("Egypt", 2)) +
    facet_zoom(y = value >= -0.25 & value <= 0.25) +ggtitle(plot_name)
  p1
})
#plot3
#It's UP/DOWN/NONE on compD/no compD (so with only chr1,17,X)
res_plot_4 <- map(1:nrow(regroup_data_per_plot),function(i){
  plot_name <- glue::glue("{regroup_data_per_plot$GroupPlot[i]}_{regroup_data_per_plot$bwType[i]}")
  p1 <- regroup_data_per_plot$data[[i]] %>% filter(str_detect(wig,"OHT")) %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() +facet_grid(~Type)  + scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2)) +ggtitle(plot_name)
  p1
})

pdf("../../results/RNA_seq_compD_ABD_genes_090222.pdf",height=8,width=8)
res_plot_1 %>% walk(print)
res_plot_2 %>% walk(print)
res_plot_3 %>% walk(print)
res_plot_4 %>% walk(print)
dev.off()

##refaire graphe5 sans zoom STX_vs_C1_bamCompare

plot_name <- glue::glue("{regroup_data_per_plot_3$GroupPlot[1]}_{regroup_data_per_plot_3$bwType[1]}")
p1 <- regroup_data_per_plot_3$data[[1]] %>% ggplot(aes(x=wig,y=value,fill=Type)) + geom_boxplot()  + scale_fill_manual(values=MetBrewer::met.brewer("Egypt", 2)) +ggtitle(plot_name)

##PVALUES
datapvals <- regroup_data_per_plot_3$data[[1]] %>% filter(str_detect(wig,"OHT")) %>%
  mutate(Type = case_when(Type == "Upregulated" ~ "Upregulated",TRUE ~ "NONE/DOWN"))


numbers <- datapvals %>% group_by(Type) %>% dplyr::count()
pvalues.dat <- wilcox.test(value~Type,data=datapvals) %>% broom::tidy()
pvalues.dat %>% left_join(numbers,by=c("Type"))


## refaire graphe 7 STX_vs_C1_bamCompare mais DOWN/NONE Combinés
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

## refaire graphe 7 STX_vs_C1_pOHT24_bwcompare en DRIP mais DOWN/NONE Combinés -1/+1 genebody

drip_bw_list <- list(
  "HKKWHBGX7_DRIP_STX_vs_CTRL_pOHT24H"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HKKWHBGX7/PROCESSED/mapping/bigwigCompare/HKKWHBGX7_DRIP_STX_vs_CTRL_pOHT24H.bigwigCompare.bw"
  
) 
drip_bw_list <- drip_bw_list %>% map(import.bw,as="RleList")

res.plot.drip_bw_list <- mclapply(names(drip_bw_list),function(one_wig){
  # my_x <- mes_genes %>% promoters(500,500)
  my_x <- mes_genes %>% anchor_start() %>% plyranges::stretch(1000)%>% anchor_end() %>% plyranges::stretch(1000)
  dat1 <- Get1valMean(Name = one_wig,one.w = drip_bw_list[[one_wig]],x = my_x)
  
  dat1
},mc.cores=length(drip_bw_list)) %>% bind_rows()
res.drip <- as_tibble(mes_genes) %>% left_join(DE_DIVA_INFO,by = c("gene_id"="rowname"))%>% right_join(res.plot.drip_bw_list,by = c("name"="rowname"))

res_comp.drip <-
  rbind(
    res.drip %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 3)) %>% as_tibble() %>% mutate(TypeD = glue::glue("CompD")),
    res.drip %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 0)) %>% as_tibble() %>% mutate(TypeD = glue::glue("NoCompD"))
  )

res_comp.drip.duplicate <- res_comp.drip %>% dplyr::count(gene_id,wig) %>% filter(n>1) %>% pull(gene_id)

res_comp.drip <- res_comp.drip %>% filter(!gene_id %in% res_comp.drip.duplicate)

p3 <- res_comp.drip %>%
  mutate(wig = str_remove_all(wig,"HKKWHBGX7_DRIP_")) %>% 
  mutate(Type = case_when(Type == "Upregulated" ~ "Upregulated",TRUE ~ "NONE/DOWN")) %>% 
  ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() +facet_grid(~Type)  + scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2)) +ggtitle("HKKWHBGX7_DRIP_STX_vs_CTRL_pOHT24H")

numbers <- res_comp.drip %>% group_by(TypeD,Type) %>% dplyr::count() %>% spread(key=TypeD,value=n)

pvalues.dat <- res_comp.drip %>% group_by(Type) %>% nest() %>% 
  mutate(data = map(data,function(x){
    wilcox.test(value~TypeD,data=x) %>% broom::tidy()
  })) %>% unnest(data) %>% mutate(p.value = format.pval(p.value,2))
pvalues.dat %>% left_join(numbers,by=c("Type"))


pdf("../../results/RNA_seq_compD_ABD_genes_090222_withDRIP_FC3.pdf",height=8,width=8)
print(p1)
print(p2)
print(p3)
dev.off()

# 
# mes_bw.rna <- list(
#   # "RNA_C1_1_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYL5LBGX2/HYL5LBGX2_RNA_C1_OHT_1_17s002465-1-1_Clouaire_lane117s002465_hg19plusERCC_normalized.bw"
#   # ,"RNA_STX_1_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYL5LBGX2/HYL5LBGX2_RNA_STX_OHT_1_17s002467-1-1_Clouaire_lane117s002467_hg19plusERCC_normalized.bw"
#   # ,"RNA_C1_2_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYMHNBGX2/HYMHNBGX2_RNA_C1_OHT_2_17s002469-1-1_Clouaire_lane117s002469_hg19plusERCC_normalized.bw"
#   # ,"RNA_STX_2_pOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/RNA-Seq/EMBL_201707/PROCESSED/WIGGLE/HYMHNBGX2/HYMHNBGX2_RNA_STX_OHT_2_17s002471-1-1_Clouaire_lane117s002471_hg19plusERCC_normalized.bw"
#   "HYL5LBGX2_RNA_STX_vs_C1_DIVA.log2.bamCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_STX_vs_C1_DIVA.log2.bamCompare.bw"
#   ,"HYL5LBGX2_RNA_STX_vs_C1_DIVA.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_STX_vs_C1_DIVA.log2.bigwigCompare.bw"
#   ,"HYL5LBGX2_RNA_STX_vs_C1_OHT.log2.bamCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_STX_vs_C1_OHT.log2.bamCompare.bw"
#   ,"HYL5LBGX2_RNA_STX_vs_C1_OHT.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_STX_vs_C1_OHT.log2.bigwigCompare.bw"
#   # ,"HYMHNBGX2_RNA_STX_vs_C1_DIVA.log2.bamCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYMHNBGX2_RNA_STX_vs_C1_DIVA.log2.bamCompare.bw"
#   # ,"HYMHNBGX2_RNA_STX_vs_C1_DIVA.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYMHNBGX2_RNA_STX_vs_C1_DIVA.log2.bigwigCompare.bw"
#   # ,"HYMHNBGX2_RNA_STX_vs_C1_OHT.log2.bamCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYMHNBGX2_RNA_STX_vs_C1_OHT.log2.bamCompare.bw"
#   # ,"HYMHNBGX2_RNA_STX_vs_C1_OHT.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYMHNBGX2_RNA_STX_vs_C1_OHT.log2.bigwigCompare.bw"
#   ,"HYL5LBGX2_RNA_OHT_vs_DIVA_STX.log2.bamCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_OHT_vs_DIVA_STX.log2.bamCompare.bw"
#   ,"HYL5LBGX2_RNA_OHT_vs_DIVA_STX.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_OHT_vs_DIVA_STX.log2.bigwigCompare.bw"
#   ,"HYL5LBGX2_RNA_OHT_vs_DIVA_CTRL.log2.bamCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_OHT_vs_DIVA_CTRL.log2.bamCompare.bw"
#   ,"HYL5LBGX2_RNA_OHT_vs_DIVA_CTRL.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_OHT_vs_DIVA_CTRL.log2.bigwigCompare.bw"
#   
# ) 
# mes_bw.rna <- mes_bw.rna %>% mclapply(import.bw,as="RleList",mc.cores=length(mes_bw.rna))
# 
# 
# ens.exons.domains <- ensembldb::exons(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,filter = AnnotationFilter::GeneIdFilter(mes_genes$gene_id))
# ens.exons.domains <- regioneR::filterChromosomes(ens.exons.domains,keep.chr=c(1:22,"X","Y"))
# seqlevels(ens.exons.domains) <- paste0("chr",seqlevels(ens.exons.domains))
# ens.exons.domains$name <- ens.exons.domains$gene_id
# 
# 
# res.plot.rna <- mclapply(names(mes_bw.rna),function(one_wig){
#   # my_x <- mes_genes %>% promoters(500,500)
#   my_x <- mes_genes
#   dat1 <- Get1valMean(Name = one_wig,one.w = mes_bw.rna[[one_wig]],x = my_x)
#   my_x <- ens.exons.domains 
#   dat2 <- Get1valMean(Name = one_wig,one.w = mes_bw.rna[[one_wig]],x = my_x)%>% group_by(wig,rowname) %>% summarise(value = mean(value))
#   rbind(
#     dat1 %>% mutate(GeneType = "GeneBody")
#     # ,dat2 %>% mutate(GeneType = "Exons")
#   )
# },mc.cores=length(mes_bw.rna)) %>% bind_rows()
# 
# res.rna <- as_tibble(mes_genes) %>% left_join(DE_DIVA_INFO,by = c("gene_id"="rowname"))%>% right_join(res.plot.rna,by = c("name"="rowname"))


#script au 080222
# res1 <- 
#   rbind(
#     res.rna %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 3)) %>% as_tibble() %>% mutate(TypeD = glue::glue("CompD")),
#     res.rna %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 0)) %>% as_tibble() %>% mutate(TypeD = glue::glue("NoCompD"))
#   )
# 
# res1.duplicate <- res1 %>% dplyr::count(gene_id,wig,GeneType) %>% filter(n>1) %>% pull(gene_id)
# 
# res1 <- res1 %>% filter(!gene_id %in% res1.duplicate)
# 
# theme_set(theme_classic(base_size=12))
# theme_update(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#  
# regroup_data_per_plot <- res1 %>% mutate(GroupPlot = case_when(str_detect(wig,"STX_vs_C1") ~ "STX_vs_C1",TRUE ~ "OHT_vs_DIVA")) %>% 
#   mutate(bwType = str_extract(wig,"bigwigCompare|bamCompare"))%>% mutate(wig = str_remove_all(wig,"HYL5LBGX2_|log2.bamCompare|log2.bigwigCompare")) %>% 
#   group_by(GroupPlot,bwType) %>% nest()
# 
# res_plot_1 <- map(1:nrow(regroup_data_per_plot),function(i){
#   plot_name <- glue::glue("{regroup_data_per_plot$GroupPlot[i]}_{regroup_data_per_plot$bwType[i]}")
#   p1 <- regroup_data_per_plot$data[[i]] %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot()  + facet_grid(~TypeD)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
#   p2 <- regroup_data_per_plot$data[[i]] %>% ggplot(aes(x=wig,y=value,fill=Type)) + geom_boxplot()
#   p3 <- regroup_data_per_plot$data[[i]] %>% ggplot(aes(x=wig,y=value,fill=Type)) + geom_boxplot() + facet_grid(~TypeD)
#   p1/p2/p3 +  plot_annotation(title = plot_name) +  plot_layout(guides = 'collect')
# })
# 
# res_plot_2 <- map(1:nrow(regroup_data_per_plot),function(i){
#   plot_name <- glue::glue("{regroup_data_per_plot$GroupPlot[i]}_{regroup_data_per_plot$bwType[i]}_chr1")
#   p1 <- regroup_data_per_plot$data[[i]] %>% filter(seqnames =="chr1")%>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot()  + facet_grid(~TypeD)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
#   p2 <- regroup_data_per_plot$data[[i]] %>% filter(seqnames =="chr1")%>% ggplot(aes(x=wig,y=value,fill=Type)) + geom_boxplot()
#   p3 <- regroup_data_per_plot$data[[i]] %>% filter(seqnames =="chr1")%>% ggplot(aes(x=wig,y=value,fill=Type)) + geom_boxplot() + facet_grid(~TypeD)
#   p1/p2/p3 +  plot_annotation(title = plot_name) +  plot_layout(guides = 'collect')
# })
# pdf("../../results/RNA_seq_compD_ABD_genes_080222.pdf",height=16,width=5)
# res_plot_1 %>% walk(print)
# res_plot_2 %>% walk(print)
# dev.off()



# p1.bamcompare <- res1 %>%
#   filter(str_detect(wig,"bamCompare")) %>%
#   ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot()  + facet_grid(~TypeD)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
# p1.bigwigCompare <- res1 %>%
#   filter(str_detect(wig,"bigwigCompare")) %>% mutate(wig = str_remove(wig,"log2.bigwigCompare")) %>%
#   ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot()  + facet_grid(~TypeD)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))


# p2.bamcompare <- res1 %>% filter(seqnames =="chr1") %>%
#   filter(str_detect(wig,"bamCompare")) %>% mutate(wig = str_remove(wig,"log2.bamCompare"))%>% 
#   ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot()  + facet_grid(GeneType~TypeD)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
# p2.bigwigCompare <- res1 %>% filter(seqnames =="chr1") %>%
#   filter(str_detect(wig,"bigwigCompare")) %>% mutate(wig = str_remove(wig,"log2.bigwigCompare")) %>% 
#   ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot()  + facet_grid(GeneType~TypeD)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
# 
# p2.classic <- res1%>% filter(seqnames =="chr1") %>%
#   filter(str_detect(wig,"bamCompare|bigwigCompare",negate = T)) %>% 
#   ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() +scale_y_log10()  + facet_wrap(~GeneType)+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))

p2.l

# pdf("../../results/RNA_seq_compD_ABD_genes_070222_bwcompare_only.pdf",height=8,width=5)
# res_plot_1 %>% walk(print)
# res_plot_2 %>% walk(print)
# dev.off()

# p2 <- res1 %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() + scale_y_log10() + facet_wrap(GeneType~Type) + theme_classic()+ scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2))
### same but with chr1 only


###sort by qDRIP signal and then plot compD


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



#Top/Bottom 20
theme_set(theme_classic(base_size=12))

pDRIP.compd <- qDRIPgenes.compD %>% mutate(Group = case_when(
  Group %in% 1:20 ~ "Bottom20",
  Group %in% 81:100 ~ "Top20",
  TRUE ~ "Other"
)) %>% filter(Group != "Other") %>% dplyr::select(seqnames,gene_name,wig:Group) %>% gather(key=CompD_Source,value=compD,-seqnames:-GeneType,-Group)

plotqDRIP1 <- pDRIP.compd %>% ggplot(aes(x=GeneType,y=compD,fill=Group)) + geom_boxplot() + facet_grid(CompD_Source~wig) + scale_fill_manual(values=MetBrewer::met.brewer("Tsimshian", 2))


#Top/Bottom 10
pDRIP.compd <- qDRIPgenes.compD %>% mutate(Group = case_when(
  Group %in% 1:10 ~ "Bottom10",
  Group %in% 91:100 ~ "Top10",
  TRUE ~ "Other"
)) %>% filter(Group != "Other") %>% dplyr::select(seqnames,gene_name,wig:Group) %>% gather(key=CompD_Source,value=compD,-seqnames:-GeneType,-Group)

plotqDRIP2 <- pDRIP.compd %>% ggplot(aes(x=GeneType,y=compD,fill=Group)) + geom_boxplot() + facet_grid(CompD_Source~wig) + scale_fill_manual(values=MetBrewer::met.brewer("Tsimshian", 2))

#Top/Bottom 100

pDRIP.compd <- qDRIPgenes.compD %>% mutate(Group = case_when(
  Group == 1 ~ "Bottom1",
  Group == 100 ~ "Top1",
  TRUE ~ "Other"
)) %>% filter(Group != "Other") %>% dplyr::select(seqnames,gene_name,wig:Group) %>% gather(key=CompD_Source,value=compD,-seqnames:-GeneType,-Group)

plotqDRIP3 <- pDRIP.compd %>% ggplot(aes(x=GeneType,y=compD,fill=Group)) + geom_boxplot() + facet_grid(CompD_Source~wig) + scale_fill_manual(values=MetBrewer::met.brewer("Tsimshian", 2))


pdf("../../results/compD_by_Drip_Cat.pdf",height=12,width=13)
print(plotqDRIP1 + ggtitle("20"))
print(plotqDRIP2 + ggtitle("10"))
print(plotqDRIP3 + ggtitle("1"))
dev.off()


#plot au 14/10 en genebody only (et en ayant viré les qDRIP) et sur top/bottom 10
#Top/Bottom 10
pDRIP.compd <- qDRIPgenes.compD %>%
  filter(GeneType == "GeneBody") %>% 
  mutate(Group = case_when(
    Group %in% 1:10 ~ "Bottom10",
    Group %in% 91:100 ~ "Top10",
    TRUE ~ "Other"
  )) %>% filter(Group != "Other") %>% dplyr::select(seqnames,gene_name,wig:Group) %>% gather(key=CompD_Source,value=compD,-seqnames:-GeneType,-Group)

plotqDRIP14 <- pDRIP.compd %>% ggplot(aes(x=GeneType,y=compD,fill=Group)) + geom_boxplot() + facet_grid(CompD_Source~wig) + scale_fill_manual(values=MetBrewer::met.brewer("Tsimshian", 2))


pdf("../../results/compD_by_Drip_Cat_14012022.pdf",height=12,width=8)
print(plotqDRIP14 + ggtitle("10"))
dev.off()

##PVALUES

numbers <- pDRIP.compd %>% group_by(CompD_Source,wig,Group) %>% dplyr::count() %>% spread(key=Group,value=n)

pvalues.dat <- pDRIP.compd %>% group_by(CompD_Source,wig) %>% nest() %>% 
  mutate(data = map(data,function(x){
    wilcox.test(compD~Group,data=x) %>% broom::tidy()
  })) %>% unnest(data) %>% mutate(p.value = format.pval(p.value,2))
pvalues.dat %>% left_join(numbers,by=c("wig","CompD_Source"))

##CHR1 ONLY

qDRIPgenes.compD <- qDRIPgenes.compD %>% filter(seqnames =="chr1") %>%  group_by(GeneType,wig) %>% arrange(desc(value)) %>% mutate(Group = ntile(value,100))

#Top/Bottom 10

pDRIP.compd <- qDRIPgenes.compD %>% mutate(Group = case_when(
  Group %in% 1:20 ~ "Bottom20",
  Group %in% 81:100 ~ "Top20",
  TRUE ~ "Other"
)) %>% filter(Group != "Other") %>% dplyr::select(seqnames,gene_name,wig:Group) %>% gather(key=CompD_Source,value=compD,-seqnames:-GeneType,-Group)

plotqDRIP1 <- pDRIP.compd %>% ggplot(aes(x=GeneType,y=compD,fill=Group)) + geom_boxplot() + facet_grid(CompD_Source~wig) + scale_fill_manual(values=MetBrewer::met.brewer("Tsimshian", 2))

pDRIP.compd <- qDRIPgenes.compD %>% mutate(Group = case_when(
  Group %in% 1:10 ~ "Bottom10",
  Group %in% 91:100 ~ "Top10",
  TRUE ~ "Other"
)) %>% filter(Group != "Other") %>% dplyr::select(seqnames,gene_name,wig:Group) %>% gather(key=CompD_Source,value=compD,-seqnames:-GeneType,-Group)

plotqDRIP2 <- pDRIP.compd %>% ggplot(aes(x=GeneType,y=compD,fill=Group)) + geom_boxplot() + facet_grid(CompD_Source~wig) + scale_fill_manual(values=MetBrewer::met.brewer("Tsimshian", 2))




pDRIP.compd <- qDRIPgenes.compD %>% mutate(Group = case_when(
  Group == 1 ~ "Bottom1",
  Group == 100 ~ "Top1",
  TRUE ~ "Other"
)) %>% filter(Group != "Other") %>% dplyr::select(seqnames,gene_name,wig:Group) %>% gather(key=CompD_Source,value=compD,-seqnames:-GeneType,-Group)

plotqDRIP3 <- pDRIP.compd %>% ggplot(aes(x=GeneType,y=compD,fill=Group)) + geom_boxplot() + facet_grid(CompD_Source~wig) + scale_fill_manual(values=MetBrewer::met.brewer("Tsimshian", 2))


pdf("../../results/compD_by_Drip_Cat_chr1.pdf",height=12,width=13)
print(plotqDRIP1 + ggtitle("20"))
print(plotqDRIP2 + ggtitle("10"))
print(plotqDRIP3 + ggtitle("1"))
dev.off()



# Faudrait une estime du siSETX sur tout ces genes en DRIP et en RNA-seq (p/m en SETX & siCTRL + SETX/CTRL en pOHT)
# les bigwigcompare des drip à 24h des drip SETX/CTRL → celui utilisé pr le papier BLM
# niveau là de drip sur les compD non compD (genebody etendu sur le drip on peut -1/+1)
# ensuite UP/DOWN/NONE Donc en résume :
#   un boxplot compD/nonCOMPD
# un boxplot  UP/DOWN/NONE
# et un boxplot combiné 
# en control mettre le mOHT à coté du pOHT
pdfname <- "siSETX_and_RNAseq_on_DIVAgenes_compD_no_compD_Up_none_down"

drip_bw_list <- list(
  "HKKWHBGX7_DRIP_STX_vs_CTRL_pOHT24H"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HKKWHBGX7/PROCESSED/mapping/bigwigCompare/HKKWHBGX7_DRIP_STX_vs_CTRL_pOHT24H.bigwigCompare.bw"
  ,"HKKWHBGX7_DRIP_STX_vs_CTRL_DIVA"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/DRIP-SEQ/DRIP-SEQ_in_siRNA/Clouaire_HKKWHBGX7/PROCESSED/mapping/bigwigCompare/HKKWHBGX7_DRIP_STX_vs_CTRL_DIVA.bigwigCompare.bw"
  
) 
drip_bw_list <- drip_bw_list %>% map(import.bw,as="RleList")

res.plot.drip_bw_list <- mclapply(names(drip_bw_list),function(one_wig){
  # my_x <- mes_genes %>% promoters(500,500)
  my_x <- mes_genes %>% anchor_start() %>% plyranges::stretch(1000)%>% anchor_end() %>% plyranges::stretch(1000)
  dat1 <- Get1valMean(Name = one_wig,one.w = drip_bw_list[[one_wig]],x = my_x)
  
  dat1
},mc.cores=length(drip_bw_list)) %>% bind_rows()



mes_bw.rna.for.specific_need <- list(
  "HYL5LBGX2_RNA_STX_vs_C1_OHT.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYL5LBGX2_RNA_STX_vs_C1_OHT.log2.bigwigCompare.bw",
  "HYMHNBGX2_RNA_OHT_vs_DIVA_CTRL.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYMHNBGX2_RNA_OHT_vs_DIVA_CTRL.log2.bigwigCompare.bw",
  "HYMHNBGX2_RNA_OHT_vs_DIVA_STX.log2.bigwigCompare"="/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/BIGWIG_BAM_COMPARE/HYMHNBGX2_RNA_OHT_vs_DIVA_STX.log2.bigwigCompare.bw"
)
mes_bw.rna.for.specific_need <- mes_bw.rna.for.specific_need %>% mclapply(import.bw,as="RleList",mc.cores=length(mes_bw.rna.for.specific_need))

ens.exons.domains <- ensembldb::exons(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,filter = AnnotationFilter::GeneIdFilter(mes_genes$gene_id))
ens.exons.domains <- regioneR::filterChromosomes(ens.exons.domains,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.exons.domains) <- paste0("chr",seqlevels(ens.exons.domains))
ens.exons.domains$name <- ens.exons.domains$gene_id


res.mes_bw.rna.for.specific_need <- lapply(names(mes_bw.rna.for.specific_need),function(one_wig){
  # my_x <- mes_genes %>% promoters(500,500)
  message(one_wig)
  my_x <- mes_genes
  dat1 <- Get1valMean(Name = one_wig,one.w = mes_bw.rna.for.specific_need[[one_wig]],x = my_x)
  my_x <- ens.exons.domains 
  dat2 <- Get1valMean(Name = one_wig,one.w = mes_bw.rna.for.specific_need[[one_wig]],x = my_x)%>% group_by(wig,rowname) %>% summarise(value = mean(value))
  rbind(
    dat1 %>% mutate(GeneType = "GeneBody")
    ,dat2 %>% mutate(GeneType = "Exons")
  )
}) %>% bind_rows()

res.rna.for.specific_need <- as_tibble(mes_genes) %>% left_join(DE_DIVA_INFO,by = c("gene_id"="rowname"))
res.rna.for.specific_need <- list(
  res.rna.for.specific_need %>% right_join(res.plot.drip_bw_list,by = c("name"="rowname")) ,
  res.rna.for.specific_need %>% right_join(res.mes_bw.rna.for.specific_need,by = c("name"="rowname"))
) %>% bind_rows()


res1 <- 
  rbind(
    res.rna.for.specific_need %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 3)) %>% as_tibble() %>% mutate(TypeD = glue::glue("CompD")),
    res.rna.for.specific_need %>% as_granges() %>% filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 0)) %>% as_tibble() %>% mutate(TypeD = glue::glue("NoCompD"))
  )

res1.duplicate <- res1 %>% dplyr::count(gene_id,wig,GeneType) %>% filter(n>1) %>% pull(gene_id)

res1 <- res1 %>% filter(!gene_id %in% res1.duplicate)


res1.bc.rnaseq <- res1 %>% filter(str_detect(wig,"RNA"))%>% mutate(wig = str_remove(wig,"log2.bigwigCompare"))
res1.bc.dripseq <- res1 %>% filter(str_detect(wig,"DRIP"))


plotboxplotcompD <- function(x){
  x  %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot() + scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2)) + theme_classic(base_size=18) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

plotboxplotupregul <- function(x){
  x  %>% ggplot(aes(x=wig,y=value,fill=Type)) + geom_boxplot() + theme_classic(base_size=18) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

pdf(glue::glue("../../results/{pdfname}_right_plots_for_gaelle_22_11_2022.pdf"),height=12,width=7)
p2 <- res1.bc.dripseq %>% ggplot(aes(x=wig,y=value,fill=TypeD)) + geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=MetBrewer::met.brewer("Thomas", 2)) + theme_classic(base_size=18) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(p2 + coord_cartesian(ylim=c(-0.05,0.05)))
dev.off()

# pvals
res1.bc.dripseq %>% group_by(wig) %>% nest() %>% mutate(data = map(data,function(x){
  tryCatch({
    wilcox.test(value~TypeD,data=x) %>% broom::tidy()
  },error=function(x){
    NA
  })
})) %>% unnest(data) %>% mutate(p.value = format.pval(p.value,3))


pdf(glue::glue("../../results/{pdfname}_right_plots_for_gaelle.pdf"),height=12,width=7)
p2 <- res1.bc.dripseq %>% plotboxplotcompD
print(p2)
dev.off()

pdf(glue::glue("../../results/{pdfname}.pdf"),height=12,width=8)
p1 <- res1.bc.rnaseq %>% plotboxplotcompD
p2 <- res1.bc.dripseq %>% plotboxplotcompD
pf1 <- p1+p2 + plot_layout(guides = "collect")

p1 <- res1.bc.rnaseq %>% plotboxplotupregul
p2 <- res1.bc.dripseq %>% plotboxplotupregul

pf2 <- p1+p2 + plot_layout(guides = "collect")

p1 <- res1.bc.rnaseq %>% plotboxplotcompD + facet_grid(Type~.)
p2 <- res1.bc.dripseq %>% plotboxplotcompD+ facet_grid(Type~.)
pf3 <- p1+p2 + plot_layout(guides = "collect")
print(pf1)
print(pf2)
print(pf3)

dev.off()





## Scatterplot bin100kb compD vs DRIP2_C1_pOHT


compD.signal <- c(
  "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHT.bw",
  "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipA.bw",
  "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipB.bw"
) %>%
  setNames(str_remove_all(basename(.),"PC1_all_chr_log2ratio_100kb_|.bw")) %>% 
  map(import.bw) %>%
  map(filter,seqnames %in% chr.to.study) %>% 
  map(as_tibble) %>% 
  map(dplyr::select,-strand,-width) %>% 
  bind_rows(.id="Name") %>% spread(key=Name,value=score) %>% 
  as_granges() %>% filter_by_non_overlaps(gamma_region)

compD.signal$name <- paste(compD.signal)
compD.signal <- compD.signal %>% sortSeqlevels() %>% sort()
compD.signal$bin <- 1:length(compD.signal)
res.signal <- lapply(names(mes_bw),function(one_wig){
  
  PhDfunc::Get1val(Name = one_wig,one.w = mes_bw[[one_wig]],x = compD.signal)
}) %>% bind_rows()

res.signal <- res.signal %>% spread(key=wig,value=value)


compD.signal <- compD.signal %>% as_tibble() %>% left_join(res.signal,by=c("name"="rowname"))


require(ggforce)
p3 <- compD.signal %>% ggplot() + geom_point(aes(x = .panel_x, y = .panel_y)) +
  facet_matrix(rows = vars(HiC_D_OHT,OHT_manipA,OHT_manipB), cols = vars(DRIP2_C1_mOHT,DRIP2_C1_pOHT,DRIP2_SETX_mOHT,DRIP2_SETX_pOHT)) + theme_classic()
pdf("../../results/scatterplot_drip_vs_compD.pdf",height=8,width=16)
print(p3)
dev.off()


## 