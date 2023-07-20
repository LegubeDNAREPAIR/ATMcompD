
require(tidyverse)
require(plyranges)
require(rtracklayer)
# PC.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHT.bw"
PC.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipB.bw"
# PC.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipA.bw"
# PC.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHTDNAPKi.bw"
# PC.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio.bw"
out_name_file <- basename(PC.path) %>% str_remove(".bw")
cc.cov <- PC.path %>% import.bw(as="RleList")
#Get profiles
#Keep only DSB within a good PCs1 in chromosome
# keep.chr <- res_pca %>% bind_rows() %>% filter(PCs == "Comp.1") %>% dplyr::select(seqnames,PCs,my_cor,p.val) %>% group_by(seqnames) %>% dplyr::slice(1) %>% arrange(desc(abs(my_cor))) %>% ungroup() %>% dplyr::slice(1:5) %>% pull(seqnames)
if(str_detect(basename(PC.path),"OHTDNAPKi")){
    keep.chr <- glue::glue("chr{c(2,6,9,13,18,1,17,20,'X')}")
}else{
    keep.chr <- glue::glue("chr{c(1,17,'X')}")
}



m.w <- "50kb"
sw <- 50000
rds.files <- list.files("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/TRANS/KR/HiTC/observed/",pattern=str_c("_",m.w),full.names=T)


DSB174 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")
bless80 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
my_dsb <- regioneR::filterChromosomes(bless80,keep.chr = keep.chr)


my_dsb.chr<- my_dsb %>%
    anchor_center() %>% mutate(width = 1000000) 

HR = import.bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed") %>% 
    regioneR::filterChromosomes(keep.chr = keep.chr) %>% 
    anchor_center() %>% mutate(width = 1000000)

NHEJ = import.bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed") %>% 
    regioneR::filterChromosomes(keep.chr = keep.chr) %>% 
    anchor_center() %>% mutate(width = 1000000)

Random.pos <- regioneR::randomizeRegions(my_dsb.chr,per.chromosome = T)
Random.pos$name <- str_c("Random",1:length(Random.pos))

Random80 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/80random.bed")  %>%
    anchor_center() %>% mutate(width = 1000000) %>% regioneR::filterChromosomes(keep.chr = keep.chr)

cc.dsb <- list("DSB"=my_dsb.chr,"Ctrl"=Random.pos,"Random80"=Random80,"NHEJ"=NHEJ,"HR"=HR) %>% 
    map(function(one_dsb){
        PhDfunc::Get1val(Name = "cc",one.w = cc.cov,x = one_dsb)
    }) %>% bind_rows(.id = "Name")

p1 <- cc.dsb  %>% 
    ggplot(aes(x=Name,y=value,fill=Name)) + geom_boxplot() + theme_classic()

#special plot for gaelle
p1.bis <- cc.dsb  %>% 
    filter(Name %in% c("Ctrl","DSB","Random80")) %>% 
    ggplot(aes(x=Name,y=value,fill=Name)) + geom_boxplot() + theme_classic()

# cc.dsb.prof <- list("DSB"=DSB174.chr,"Ctrl"=Random.pos,"Random80"=Random80,"NHEJ"=NHEJ,"HR"=HR) %>% 
#     map(function(one_dsb){
#         PhDfunc::computeProfileAugmented(one_dsb,wig = cc.cov,w = 500000,span = 2000,seqlens = seqlens,heatmap = T)
#     }) %>% bind_rows(.id = "Name")
# 
# p2 <- cc.dsb.prof %>%
#     group_by(Windows,Name) %>% summarise(meanVal = mean(value,na.rm = T)) %>% 
#     ggplot(aes(x= Windows,y=meanVal,col=Name)) + geom_line() + theme_classic()

pdf(glue::glue("../../results/PC1_signal_over_DSB_in_same_chr_{out_name_file}.pdf"),height=4,width=8)
print(p1)
# print(p2)
dev.off()

pdf(glue::glue("../../results/PC1_signal_over_80bless_in_same_chr_{out_name_file}.pdf"),height=4,width=8)
print(p1.bis)
# print(p2)
dev.off()

cc.dsb %>% group_by(Name) %>% nest() %>% mutate(pval = map_dbl(data,function(x){wilcox.test(x$value)$p.value}))


DE_DIVA <- PhDfunc::GetDE_DIvA()

DE_DIVA <- DE_DIVA %>% 
  mutate(Type = case_when(
    # logFC < 0 & FILTER.P == 1 ~ "Downregulated",
    logFC < 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Downregulated",
    # logFC < 0 & FILTER.FC == 1 ~ "Downregulated",
    # logFC < -0.3 ~ "Downregulated",
    # logFC < -0.15 ~ "Downregulated",
    logFC > 0 & FILTER.FC == 1 & FILTER.P == 1~ "Upregulated",
    # logFC > 0 & FILTER.P == 1~ "Upregulated",
    # logFC > 0 & FILTER.FC == 1~ "Upregulated",
    # logFC > 0.15 ~ "Upregulated",
    # logFC > 0.3 ~ "Upregulated",
    TRUE ~ "None"
  )) 

ens.genes <- read.table("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed",sep="\t",h=T) %>% GRanges()
ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))
mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
    arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
# gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 1000000)
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 2000000)

# %>% sample(300)

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


# mes_genes <- list("up.genes"=up.genes,"down.genes"=down.genes,"no.genes"=no.genes)
mes_genes <- list("up.genes"=up.genes,"no.genes"=no.genes)

# PC.path <- 
PC.path.multi <- c("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipB.bw",
                   "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHT.bw",
                   "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipA.bw")
names(PC.path.multi) <- PC.path.multi %>% map(function(i){basename(i) %>% str_remove(".bw")})
keep.chr <- glue::glue("chr{c(1,17,'X')}")

up.genes <- ens.genes %>% filter(gene_id %in% filter(DE_DIVA,Type =="Upregulated")$rowname)%>%
  filter_by_non_overlaps(gamma_region) %>% 
  regioneR::filterChromosomes(keep.chr = keep.chr)
down.genes <- ens.genes %>% filter(gene_id %in% filter(DE_DIVA,Type =="Downregulated")$rowname)%>%
  filter_by_non_overlaps(gamma_region) %>% 
  regioneR::filterChromosomes(keep.chr = keep.chr)
# no.genes <- ens.genes %>% filter(gene_id %in% filter(DE_DIVA,Type =="Upregulated",FILTER.FC != 1)$rowname) %>%
#     filter_by_non_overlaps(gamma_region) %>% 
#     regioneR::filterChromosomes(keep.chr = keep.chr)


no.genes <- ens.genes %>% filter(gene_id %in% filter(DE_DIVA,Type !="Upregulated")$rowname) %>%
  filter_by_non_overlaps(gamma_region) %>% 
  regioneR::filterChromosomes(keep.chr = keep.chr)

ccdsb.pre <- PC.path.multi %>%  map(function(one_compD){
  cc.cov <- one_compD %>% import.bw(as="RleList")
  cc.dsb <- mes_genes %>% 
    map(function(one_dsb){
      one_dsb$name <- one_dsb$gene_id
      Get1val(Name = "cc",one.w = cc.cov,x = one_dsb)
    }) %>% bind_rows(.id = "Name")
  
}) %>% bind_rows(.id = "TypeofCompD")

cc.dsb <- ccdsb.pre %>% group_by(Name,rowname) %>% summarise(value = mean(value))

p3 <- cc.dsb  %>% 
    mutate(Name = factor(Name,levels = c("no.genes","up.genes","down.genes"))) %>% 
    ggplot(aes(x=Name,y=value,fill=Name)) + geom_boxplot(outlier.shape = NA) + theme_classic()

pdf(glue::glue("../../results/after_FC03/PC1_signal_over_DE_DIVA_Genes_in_same_chr_inAverage_A_B_D.pdf"),height=4,width=8)
# print(p3)
print(p3 + coord_cartesian(ylim = c(-0.025,0.025)))
dev.off()


## WITH DNAPKI
keep.chr <- glue::glue("chr{c(2,6,9,13,18,1,17,20,'X')}")

PC.path.multi <- c("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHTDNAPKi.bw")
names(PC.path.multi) <- PC.path.multi %>% map(function(i){basename(i) %>% str_remove(".bw")})

up.genes <- ens.genes %>% filter(gene_id %in% filter(DE_DIVA,Type =="Upregulated")$rowname)%>%
  filter_by_non_overlaps(gamma_region) %>% 
  regioneR::filterChromosomes(keep.chr = keep.chr)
down.genes <- ens.genes %>% filter(gene_id %in% filter(DE_DIVA,Type =="Downregulated")$rowname)%>%
  filter_by_non_overlaps(gamma_region) %>% 
  regioneR::filterChromosomes(keep.chr = keep.chr)

ccdsb.pre <- PC.path.multi %>%  map(function(one_compD){
  cc.cov <- one_compD %>% import.bw(as="RleList")
  cc.dsb <- mes_genes %>% 
    map(function(one_dsb){
      one_dsb$name <- one_dsb$gene_id
      Get1val(Name = "cc",one.w = cc.cov,x = one_dsb)
    }) %>% bind_rows(.id = "Name")
  
}) %>% bind_rows(.id = "TypeofCompD")

cc.dsb <- ccdsb.pre %>% group_by(Name,rowname) %>% summarise(value = mean(value))

p3 <- cc.dsb  %>% 
  mutate(Name = factor(Name,levels = c("no.genes","up.genes","down.genes"))) %>% 
  ggplot(aes(x=Name,y=value,fill=Name)) + geom_boxplot(outlier.shape = NA) + theme_classic()

pdf(glue::glue("../../results/after_FC03/PC1_signal_over_DE_DIVA_Genes_in_same_chr_{paste(names(PC.path.multi),sep='_')}.pdf"),height=4,width=8)
# print(p3)
print(p3 + coord_cartesian(ylim = c(-0.03,0.03)))
dev.off()

#Only chr 1 / 17

up.genes <- ens.genes %>% filter(gene_id %in% filter(DE_DIVA,Type =="Upregulated")$rowname)%>%
    filter_by_non_overlaps(gamma_region) %>% 
    regioneR::filterChromosomes(keep.chr = c("chr1","chr17"))
down.genes <- ens.genes %>% filter(gene_id %in% filter(DE_DIVA,Type =="Downregulated")$rowname)%>%
    filter_by_non_overlaps(gamma_region) %>% 
    regioneR::filterChromosomes(keep.chr = c("chr1","chr17"))
no.genes <- ens.genes %>% filter(gene_id %in% filter(DE_DIVA,Type =="None",FILTER.FC != 1)$rowname) %>%
    filter_by_non_overlaps(gamma_region) %>% 
    regioneR::filterChromosomes(keep.chr = c("chr1","chr17")) %>% sample(300)
no.genes <- ens.genes %>% filter(gene_id %in% filter(DE_DIVA,Type =="None",FILTER.FC != 1)$rowname) %>%
    filter_by_non_overlaps(gamma_region) %>% 
    regioneR::filterChromosomes(keep.chr = c("chr1","chr17")) %>% sample(300)

mes_genes <- list("up.genes"=up.genes,"down.genes"=down.genes,"no.genes"=no.genes)



cc.dsb <- mes_genes %>% 
    map(function(one_dsb){
        one_dsb$name <- one_dsb$gene_id
        Get1val(Name = "cc",one.w = cc.cov,x = one_dsb)
    }) %>% bind_rows(.id = "Name")

p3 <- cc.dsb  %>% 
    mutate(Name = factor(Name,levels = c("no.genes","up.genes","down.genes"))) %>% 
    ggplot(aes(x=Name,y=value,fill=Name)) + geom_boxplot() + theme_classic()

pdf(glue::glue("../../results/PC1_signal_over_DE_DIVA_Genes_in_chr1_chr17_{out_name_file}.pdf"),height=4,width=8)
print(p3)
dev.off()

##MODIF AU 27/10/21 TO CHECK IF GAMMA_DOMAIN 2Mb change anything
cc.dsb %>% dplyr::select(-wig) %>% rename(compD="value")%>% write_tsv(glue::glue("../../results/genes_compD_value_used_to_produce_PC1_signal_over_DE_DIVA_Genes_in_chr1_chr17_{out_name_file}.tsv"))
pdf(glue::glue("../../results/PC1_signal_over_DE_DIVA_Genes_in_chr1_chr17_{out_name_file}_with2mbgammadomains.pdf"),height=4,width=8)
print(p3)
dev.off()

cc.dsb %>% group_by(Name) %>% nest() %>% mutate(pval = map_dbl(data,function(x){wilcox.test(x$value)$p.value}))



cc.dsb %>% group_by(Name) %>% nest() %>% mutate(pval = map_dbl(data,function(x){wilcox.test(x$value)$p.value}))


ExtractSubMatrixforGenes <- function(file, bed1,bed2) {
    c(my.ranges,intdata) %<-% readRDS(file)
    
    val1 <- intdata %>% sum()
    overlaps1 <- findOverlaps(bed1,my.ranges)
    overlaps2 <- findOverlaps(bed2,my.ranges)
    
    resSub <- intdata[subjectHits(overlaps1),subjectHits(overlaps2)] %>% colMeans()
    # res <- my.ranges[-subjectHits(overlaps)]
    res <- my.ranges[subjectHits(overlaps2)]
    res$SumRes <- resSub
    overlaps.l <- as(findOverlaps(bed2,res), "List")
    res <- overlaps.l %>% lapply(function(x){
        res[x]$SumRes %>% mean()
    }) %>% unlist()
    bed2$score <- res
    return(list(bed2,val1))
}

feed_bed <- mes_genes %>% map(as_tibble) %>% bind_rows(.id="Type") %>% as_granges()
gamma_region80 <- bless80 %>% anchor_center() %>% mutate(width = 1000000)
total_res <- rds.files %>% parallel::mclapply(ExtractSubMatrixforGenes,bed1=gamma_region80,bed2=feed_bed,mc.cores=length(rds.files))

names(total_res) <- rds.files %>% basename() %>% str_remove_all("HTC_observed_KR_|_full_matrix_[0-9]+kb.rds")
total_nbread <- total_res %>% map_dbl(2) %>% enframe()
total_table  <- total_res %>% map(1) %>% map(as_tibble) %>% bind_rows(.id = "name") %>% 
    left_join(total_nbread,by = "name")


total_table <- total_table %>% as_tibble() %>% 
    mutate(Condition = str_remove(name,"HiC_D_")) %>% 
    dplyr::select(-width,-name) %>% 
    pivot_wider(names_from = Condition, 
                values_from = c("value", "score")) %>%
    mutate(score_OHT = score_OHT * (value_DIvA/value_OHT)) %>%
    mutate(score_OHTATMi = score_OHTATMi * (value_DIvA/value_OHTATMi)) %>%
    mutate(score_OHTDNAPKi = score_OHTDNAPKi * (value_DIvA/value_OHTDNAPKi)) %>%
    mutate(score_OHTPARPi = score_OHTPARPi * (value_DIvA/value_OHTPARPi)) %>%
    dplyr::select(-value_DIvA,-value_OHT,-value_OHTDNAPKi,-value_OHTATMi,-value_OHTPARPi)
# Compute ratio

total_table <- total_table %>% 
    gather(key = Condition,value = value,-seqnames:-seq_coord_system) %>% 
    mutate(value = value +1) %>% 
    spread(key = Condition,value = value) %>%
    gather(key = Condition,value = value,-seqnames:-score_DIvA) %>% 
    mutate(ratio = log2(value/score_DIvA))


p.ratio <- total_table %>%
    drop_na()%>% ggplot(aes(x=Type,y=ratio,fill=Type)) + geom_boxplot() + facet_wrap(~Condition) + theme_classic(base_size=22)
ggsave(glue::glue("results/ratio_interaction_signal_over_genes_chr1_chr17_{m.w}_manipD.pdf"),p.ratio,width=12,height=12)





#On TP53 genes
require(motifmatchr)
require(TFBSTools)
require(JASPAR2018)
searchMotif <- function(motif_name = "TP53",mybed,extend = 3000){
    mybed <- mybed %>% promoters(upstream = extend,downstream = extend)
    opts <- list()
    opts[["species"]] <- 9606
    opts[["name"]] <- motif_name
    opts[["collection"]] <- "CORE"
    
    PFMatrixList <- getMatrixSet(JASPAR2018, opts)
    seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, mybed)
    
    motif_ix <- matchMotifs(PFMatrixList, seq,out = c("scores"),p.cutoff = 1e-03) 
    motif_id <- colnames(motif_ix)
    res <- motifMatches(motif_ix)
    counts <- motifCounts(motif_ix) %>% as.vector()
    return(counts)
}
mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
    arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 1000000)
p53_target <- read_csv("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/p53_target/p53_targets_up_regulated_in_at_least_one_dataset.csv")
# p53_target_bis <- read_csv("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/p53_target/p53_targets_up_regulated_in_U20S_at_least_one_dataset.csv")

p53_target.genes <- ens.genes %>% filter(gene_id %in% p53_target$`Ensembl ID`) %>% 
    filter_by_non_overlaps(gamma_region) %>% regioneR::filterChromosomes(keep.chr = keep.chr)

p53_direct_target.genes <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/direct_target_p53_hg19.bed") %>% 
    filter_by_non_overlaps(gamma_region) %>% regioneR::filterChromosomes(keep.chr = keep.chr) %>% mutate(gene_id = name)

# p53_target_bis.genes <- ens.genes %>% filter(gene_id %in% p53_target_bis$`Ensembl ID`)
p53_motifs.genes <- ens.genes %>% 
    filter_by_non_overlaps(gamma_region) %>% 
    mutate(occ_motif = searchMotif(mybed=.)) %>% regioneR::filterChromosomes(keep.chr = keep.chr)

p53_motifs.genes <- p53_motifs.genes %>% filter(occ_motif > 10)

# filter.genes <- ens.genes %>%
#     filter(!gene_id %in% p53_target.genes$gene_id) %>%
#     # filter(!gene_id %in% p53_target_bis.genes$gene_id) %>%
#     filter(!gene_id %in% p53_motifs.genes$gene_id) %>%
#     filter(gene_biotype == "protein_coding") %>% 
#     mutate(occ_motif = searchMotif(mybed=.)) %>% regioneR::filterChromosomes(keep.chr = keep.chr) %>% filter_by_non_overlaps(gamma_region)
# control.genes <- regioneR::resampleRegions(A=p53_target.genes,universe = (filter.genes %>% filter(occ_motif <= 7)))
set.seed(1234)
control.genes <- ens.genes %>% 
    filter(!gene_id %in%c(p53_motifs.genes$gene_id,p53_direct_target.genes$gene_id,p53_target.genes$gene_id)) %>% 
    filter_by_non_overlaps(gamma_region) %>% 
    filter(gene_id %in% DE_DIVA$rowname) %>% sample(1000)
cc.dsb <- list(p53_motifs.genes,
               p53_target.genes,
               p53_direct_target.genes,
               # DSB174.chr %>% mutate(gene_id = name),
               control.genes) %>%
    setNames(c(glue::glue("p53 MOTIF genes ({length(p53_motifs.genes)})"),
               glue::glue("p53 KNOWN genes ({length(p53_target.genes)})"),
               glue::glue("p53 DIRECT TARGET genes ({length(p53_direct_target.genes)})"),
               # glue::glue("DSB174 ({length(DSB174.chr)})"),
               glue::glue("Control genes ({length(control.genes)})"))) %>% 
    map(function(one_dsb){
        one_dsb$name <- one_dsb$gene_id
        Get1val(Name = "cc",one.w = cc.cov,x = one_dsb)
    }) %>% bind_rows(.id = "Name")

require(ggforce)
require(ggtext)
p4 <- cc.dsb  %>% 
    # mutate(Name = factor(Name,levels = c("no.genes","up.genes","down.genes"))) %>% 
    ggplot(aes(x=Name,y=value,fill=Name)) + geom_boxplot() + theme_classic() +
    facet_zoom(ylim = c(-0.025,0.025))+theme(axis.text.x = element_markdown(angle=45,hjust=1),strip.text = element_markdown(),legend.position="none")

LOG2FC.genes <- ens.genes %>% as_tibble() %>% dplyr::select(gene_id,gene_name) %>% 
    right_join(DE_DIVA,by = c("gene_id"="rowname"))


p.4.5 <- list(p53_motifs.genes %>% mutate(gene_id = gene_name),
              p53_target.genes %>% mutate(gene_id = gene_name),
              p53_direct_target.genes,
              control.genes %>% mutate(gene_id = gene_name)) %>% map(as_tibble) %>%
    setNames(c(glue::glue("p53 MOTIF genes"),
               glue::glue("p53 KNOWN genes"),
               glue::glue("p53 DIRECT TARGET genes"),
               glue::glue("Control genes"))) %>%
    map(left_join,LOG2FC.genes,by = c("gene_id"="gene_name")) %>%
    map(as_granges) %>%
    map(function(one_dsb){
        one_dsb$name <- one_dsb$gene_id
        res <- Get1val(Name = "cc",one.w = cc.cov,x = one_dsb) %>% dplyr::select(-wig)
        one_dsb %>% as_tibble() %>% left_join(res,by = c("name"="rowname"))
    }) %>% bind_rows(.id = "Name")


p.4.5 <- p.4.5  %>% 
    filter(!is.na(logFC)) %>% 
    filter(abs(logFC)>0.5) %>% mutate(Type = ifelse(logFC>0,"Upregulated","Downregulated")) %>% 
    # mutate(Name = factor(Name,levels = c("no.genes","up.genes","down.genes"))) %>% 
    ggplot(aes(x=Name,y=value,fill=Type)) + geom_boxplot() + theme_classic() +
    facet_zoom(ylim = c(-0.025,0.025))+theme(axis.text.x = element_markdown(angle=45,hjust=1),strip.text = element_markdown(),legend.position="bottom")
pdf(glue::glue("../results/PC1_signal_over_p53genes_{out_name_file}.pdf"),height=4,width=8)
print(p4)
print(p.4.5)
dev.off()

#Try with different control -> genes subsampled other than p53 targets
set.seed(1234)
control.genes <- ens.genes %>% 
    filter(!gene_id %in%c(p53_motifs.genes$gene_id,p53_direct_target.genes$gene_id,p53_target.genes$gene_id)) %>% 
    filter_by_non_overlaps(gamma_region) %>% 
    filter(gene_id %in% DE_DIVA$rowname) %>% sample(length(p53_direct_target.genes))
# cc.dsb <- list(p53_motifs.genes %>% mutate(gene_id = gene_name),
#                p53_target.genes %>% mutate(gene_id = gene_name),
#                p53_direct_target.genes,
#                control.genes %>% mutate(gene_id = gene_name)) %>% map(as_tibble) %>% 
#     setNames(c(glue::glue("p53 MOTIF genes"),
#                glue::glue("p53 KNOWN genes"),
#                glue::glue("p53 DIRECT TARGET genes"),
#                glue::glue("Control genes"))) %>% 
#     map(left_join,LOG2FC.genes,by = c("gene_id"="gene_name")) %>%
#     map(as_granges) %>% 
#     map(function(one_dsb){
#         one_dsb$name <- one_dsb$gene_id
#         res <- Get1val(Name = "cc",one.w = cc.cov,x = one_dsb) %>% dplyr::select(-wig)
#         one_dsb %>% as_tibble() %>% left_join(res,by = c("name"="rowname"))
#     }) %>% bind_rows(.id = "Name")

cc.dsb <- list(
               p53_target.genes %>% mutate(gene_id = gene_name)) %>% map(as_tibble) %>% 
    setNames(c(
               glue::glue("p53 KNOWN genes"))) %>% 
    map(left_join,LOG2FC.genes,by = c("gene_id"="gene_name")) %>%
    map(as_granges) %>% 
    map(function(one_dsb){
        one_dsb$name <- one_dsb$gene_id
        res <- Get1val(Name = "cc",one.w = cc.cov,x = one_dsb) %>% dplyr::select(-wig)
        one_dsb %>% as_tibble() %>% left_join(res,by = c("name"="rowname"))
    }) %>% bind_rows(.id = "Name")

p.4.5 <- cc.dsb  %>% 
    filter(!is.na(logFC)) %>% 
    mutate(Type = case_when(
        logFC < 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Downregulated",
        logFC > 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Upregulated",
        TRUE ~ "None"
    )) %>% 
    # mutate(Name = factor(Name,levels = c("no.genes","up.genes","down.genes"))) %>% 
    ggplot(aes(x=Name,y=value,fill=Type)) + geom_boxplot() + theme_classic() +
    facet_zoom(ylim = c(-0.025,0.025))+theme(axis.text.x = element_markdown(angle=45,hjust=1),strip.text = element_markdown(),legend.position="bottom")
pdf(glue::glue("../results/PC1_signal_over_p53genes_{out_name_file}_with_ctrl_same_size.pdf"),height=4,width=8)
print(p.4.5)
dev.off()

cc.dsb  %>% 
    filter(!is.na(logFC)) %>% 
    mutate(Type = case_when(
        logFC < 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Downregulated",
        logFC > 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Upregulated",
        TRUE ~ "None"
    )) %>% group_by(Type) %>% nest() %>% mutate(pval = map_dbl(data,function(x){wilcox.test(x$value)$p.value}))


feed_bed <- p53_target.genes %>% mutate(gene_id = gene_name) %>% as_tibble() %>% left_join(dplyr::select(LOG2FC.genes,gene_name,logFC,FILTER.FC,FILTER.P),by = c("gene_id"="gene_name"))  %>% 
    filter(!is.na(logFC)) %>% 
    mutate(Type = case_when(
        logFC < 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Downregulated",
        logFC > 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Upregulated",
        TRUE ~ "None"
    )) %>% as_granges()

gamma_region80 <- bless80 %>% anchor_center() %>% mutate(width = 1000000)
total_res <- rds.files %>% parallel::mclapply(ExtractSubMatrixforGenes,bed1=gamma_region80,bed2=feed_bed,mc.cores=length(rds.files))

names(total_res) <- rds.files %>% basename() %>% str_remove_all("HTC_observed_KR_|_full_matrix_[0-9]+kb.rds")
total_nbread <- total_res %>% map_dbl(2) %>% enframe()
total_table  <- total_res %>% map(1) %>% map(as_tibble) %>% bind_rows(.id = "name") %>% 
    left_join(total_nbread,by = "name")


total_table <- total_table %>% as_tibble() %>% 
    mutate(Condition = str_remove(name,"HiC_D_")) %>% 
    dplyr::select(-width,-name) %>% 
    pivot_wider(names_from = Condition, 
                values_from = c("value", "score")) %>%
    mutate(score_OHT = score_OHT * (value_DIvA/value_OHT)) %>%
    mutate(score_OHTATMi = score_OHTATMi * (value_DIvA/value_OHTATMi)) %>%
    mutate(score_OHTDNAPKi = score_OHTDNAPKi * (value_DIvA/value_OHTDNAPKi)) %>%
    mutate(score_OHTPARPi = score_OHTPARPi * (value_DIvA/value_OHTPARPi)) %>%
    dplyr::select(-value_DIvA,-value_OHT,-value_OHTDNAPKi,-value_OHTATMi,-value_OHTPARPi)
# Compute ratio

total_table <- total_table %>% 
    gather(key = Condition,value = value,-seqnames:-Type) %>% 
    mutate(value = value +1) %>% 
    spread(key = Condition,value = value) %>%
    gather(key = Condition,value = value,-seqnames:-score_DIvA) %>% 
    mutate(ratio = log2(value/score_DIvA))


p.ratio <- total_table %>%
    drop_na()%>% ggplot(aes(x=Type,y=ratio,fill=Type)) + geom_boxplot() + facet_wrap(~Condition) + theme_classic(base_size=22)
ggsave(glue::glue("results/ratio_interaction_signal_over_p53genes_manipD.pdf"),p.ratio,width=12,height=12)


#with voie nfkb

nfkb_path <- read_tsv("../../data/voie_nfkb.tsv",col_names = F) %>% pull(1)
library("biomaRt")

#use ensembl mart
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

my_refseq_loci <- getBM(attributes=c("refseq_mrna","external_gene_name","ensembl_gene_id"),
                        filters = "refseq_mrna",
                        values = nfkb_path,
                        mart = ensembl)

nfkb_genes <- ens.genes %>% 
    filter(gene_name %in%my_refseq_loci$external_gene_name) %>%
    # filter(gene_id %in%my_refseq_loci$ensembl_gene_id) %>% 
    filter_by_non_overlaps(gamma_region) %>% 
    regioneR::filterChromosomes(keep.chr = keep.chr)

cc.dsb <- list(
    nfkb_genes %>% mutate(gene_id = gene_name)) %>% map(as_tibble) %>% 
    setNames(c(
        glue::glue("NFKB KNOWN genes"))) %>% 
    map(left_join,LOG2FC.genes,by = c("gene_id"="gene_name")) %>%
    map(as_granges) %>% 
    map(function(one_dsb){
        one_dsb$name <- one_dsb$gene_id
        res <- Get1val(Name = "cc",one.w = cc.cov,x = one_dsb) %>% dplyr::select(-wig)
        one_dsb %>% as_tibble() %>% left_join(res,by = c("name"="rowname"))
    }) %>% bind_rows(.id = "Name")

p.4.5 <- cc.dsb  %>% 
    filter(!is.na(logFC)) %>% 
    mutate(Type = case_when(
        logFC < 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Downregulated",
        logFC > 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Upregulated",
        TRUE ~ "None"
    )) %>% 
    # mutate(Name = factor(Name,levels = c("no.genes","up.genes","down.genes"))) %>% 
    ggplot(aes(x=Name,y=value,fill=Type)) + geom_boxplot() + theme_classic() +
    facet_zoom(ylim = c(-0.025,0.025))+theme(axis.text.x = element_markdown(angle=45,hjust=1),strip.text = element_markdown(),legend.position="bottom")
pdf(glue::glue("../../results/PC1_signal_over_nfkb_genes_{out_name_file}_with_ctrl_same_size.pdf"),height=4,width=8)
print(p.4.5)
dev.off()


cc.dsb  %>% 
    filter(!is.na(logFC)) %>% 
    mutate(Type = case_when(
        logFC < 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Downregulated",
        logFC > 0 & FILTER.FC == 1 & FILTER.P == 1 ~ "Upregulated",
        TRUE ~ "None"
    )) %>% group_by(Type) %>% nest() %>% mutate(pval = map_dbl(data,function(x){wilcox.test(x$value)$p.value}))

#Take DE GENES in p53
DE_DIVA.GR <- ens.genes %>% as_tibble() %>% dplyr::select(seqnames,start,end,gene_id) %>% right_join(DE_DIVA,by = c("gene_id"="rowname")) %>% drop_na()%>% as_granges()
DE_DIVA.GR <- DE_DIVA.GR %>%
    as_tibble() %>% 
    mutate(p53_Group = ifelse(gene_id%in% p53_target.genes$gene_id,"p53 (KNOWN)","no p53")) %>% 
    mutate(color_p53 = ifelse(p53_Group == "p53","#8e44ad","black")) %>% 
    mutate(color_Type = case_when(
        Type == "Upregulated" ~ "#27ae60",
        Type == "Downregulated" ~ "#c0392b",
        TRUE ~ "#2c3e50"        
    )) %>% 
    mutate(ManualGroups = glue::glue("<b style='color:{color_p53}'>{p53_Group}</b><br><b style='color:{color_Type}'>{Type}</b>")) %>% 
    group_by(ManualGroups) %>% mutate(ManualGroups = glue::glue("{ManualGroups} ({n()})")) %>% as_granges()


cc.dsb <- split(DE_DIVA.GR,DE_DIVA.GR$ManualGroups) %>% 
    lapply(function(one_dsb){
        one_dsb$name <- one_dsb$gene_id
        Get1val(Name = "cc",one.w = cc.cov,x = one_dsb)
    }) %>% bind_rows(.id = "Name")

require(ggforce)
p5 <- cc.dsb  %>% 
    # mutate(Name = factor(Name,levels = c("no.genes","up.genes","down.genes"))) %>% 
    ggplot(aes(x=Name,y=value,fill=Name)) + geom_boxplot() + theme_classic() +
    facet_zoom(ylim = c(-0.05,0.05)) +theme(axis.text.x = element_markdown(angle=45,hjust=1),strip.text = element_markdown(),legend.position="none")
#same but with p53 motifs
DE_DIVA.GR <- ens.genes %>% as_tibble() %>% dplyr::select(seqnames,start,end,gene_id) %>% right_join(DE_DIVA,by = c("gene_id"="rowname")) %>% drop_na()%>% as_granges()
DE_DIVA.GR <- DE_DIVA.GR %>%
    as_tibble() %>% 
    mutate(p53_Group = ifelse(gene_id%in% p53_motifs.genes$gene_id,"p53 (MOTIF)","no p53")) %>% 
    mutate(color_p53 = ifelse(p53_Group == "p53","#8e44ad","black")) %>% 
    mutate(color_Type = case_when(
        Type == "Upregulated" ~ "#27ae60",
        Type == "Downregulated" ~ "#c0392b",
        TRUE ~ "#2c3e50"        
    )) %>% 
    mutate(ManualGroups = glue::glue("<b style='color:{color_p53}'>{p53_Group}</b><br><b style='color:{color_Type}'>{Type}</b>")) %>% 
    group_by(ManualGroups) %>% mutate(ManualGroups = glue::glue("{ManualGroups} ({n()})")) %>% as_granges()


cc.dsb <- split(DE_DIVA.GR,DE_DIVA.GR$ManualGroups) %>% 
    lapply(function(one_dsb){
        one_dsb$name <- one_dsb$gene_id
        Get1val(Name = "cc",one.w = cc.cov,x = one_dsb)
    }) %>% bind_rows(.id = "Name")

require(ggforce)
p6 <- cc.dsb  %>% 
    # mutate(Name = factor(Name,levels = c("no.genes","up.genes","down.genes"))) %>% 
    ggplot(aes(x=Name,y=value,fill=Name)) + geom_boxplot() + theme_classic() +
    facet_zoom(ylim = c(-0.05,0.05)) +theme(axis.text.x = element_markdown(angle=45,hjust=1),strip.text = element_markdown(),legend.position="none")

pdf(glue::glue("../results/PC1_signal_over_p53genes_{out_name_file}.pdf"),height=4,width=8)
print(p4)
print(p5)
print(p6)
dev.off()

#Take RNA seq and 3 category by log2FC signal by bin, and after that plot the PCs1
require(ggforce)
RNA_C1_pvsmOHT="/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/BAMCOMPARE/HYL5LBGX2_RNA_C1_bamcompare_logFC_normalized.bw" %>% import.bw(as="RleList")
PCs.ratio <- import.bw(PC.path)
PCs.ratio$name <- str_c("bin",1:length(PCs.ratio))
DSB174 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")
mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
    arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 1000000)
PCs.ratio.nogamma <- PCs.ratio %>% filter_by_non_overlaps(gamma_region)
res <- Get1val(Name = "binsPC1",x = PCs.ratio.nogamma,one.w = RNA_C1_pvsmOHT)
res <- res %>%
    left_join(as_tibble(PCs.ratio.nogamma) %>% dplyr::select(seqnames,score,name),by = c("rowname"="name")) %>%
    mutate(Group_PC1 = dplyr::ntile(value,5)) %>% group_by(Group_PC1) %>% mutate(meanval = round(mean(value,5),3))%>% mutate(Group_PC1_b = glue::glue("D Comp. {Group_PC1}"))
p.RNA.bc <- res %>% ggplot(aes(x=fct_reorder(Group_PC1_b,Group_PC1),y=score,fill=meanval)) + geom_boxplot() +scale_fill_viridis_c() + facet_zoom(ylim = c(-0.02,0.02)) + theme(axis.text.x = element_text(angle=45,hjust=1))

pdf(glue::glue("../results/RNAC1_bamcompare_logFC_by_compDgroups_{out_name_file}.pdf"),height=4,width=10)
print(p.RNA.bc)
dev.off()

# se mettre ds des conditions RNA-seq donc filtrer sur les position ou y'a des exons


library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
Exon <- exons(edb)
sequence <- seqlevels(Exon)
newStyle <- mapSeqlevels(sequence,"UCSC")
newStyle <- newStyle[complete.cases(newStyle)] # removing NA cases.
## rename the seqlevels
Exon <- renameSeqlevels(Exon,newStyle)
Exon.PCs.ratio <- subsetByOverlaps(Exon,PCs.ratio.nogamma,type="within") %>% plyranges::reduce_ranges()
Exon.PCs.ratio$name <- str_c("Exon",1:length(Exon.PCs.ratio))
# filtrer bin aussi pour bin negatifs bins positifs 
res <- Get1val(Name = "RNA_C1_pvsmOHT",x = Exon.PCs.ratio,one.w = RNA_C1_pvsmOHT)
res2 <- Get1val(Name = "CompD",x = Exon.PCs.ratio,one.w = coverage(PCs.ratio.nogamma,weight="score"))
res.f <- rbind(res,res2) %>%  spread(key = wig,value = value) %>% 
    left_join(as_tibble(Exon.PCs.ratio) ,by = c("rowname"="name")) %>%
    mutate(Group_PC1 = dplyr::ntile(CompD,5)) %>% group_by(Group_PC1) %>% mutate(meanval = round(mean(CompD,5),3))%>% mutate(Group_PC1_b = glue::glue("D Comp. {Group_PC1}"))
p.RNA.bc <- res.f %>% ggplot(aes(x=fct_reorder(Group_PC1_b,Group_PC1),y=RNA_C1_pvsmOHT,fill=meanval)) + geom_boxplot() +scale_fill_viridis_c() + facet_zoom(ylim = c(-0.02,0.02)) + theme(axis.text.x = element_text(angle=45,hjust=1))

pdf(glue::glue("../results/RNAC1_bamcompare_logFC_by_compDgroups_EXONS_{out_name_file}.pdf"),height=4,width=10)
print(p.RNA.bc)
dev.off()







#With Hic count directly








