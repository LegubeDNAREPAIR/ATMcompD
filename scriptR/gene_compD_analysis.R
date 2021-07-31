# gènes qui dans les 3 conditions A B ET D OHT pour le 1,17,X et 20 et qui sont positifs dans les 3 conditions :
#   - bins qui ont un niveau positif de façon reproductible dans les 3 manips, les gènes associé et combien y'en a
# - go pathway sur les gènes et les motifs
# - revoir leur expression m et plus oht
# 
# - pareil sur dnapki mais en prenant tous les chromosomes (1,2,3,6,7,8,9,13,17,18,20,21,X)
# 
# - refaire fig 6 B et D refaire sur les nouveaux set de chromosome et en A/B/D
# 
# 
# - figure 6E et supp p53 NFKPB en prenant les 1,2,3,6,7,8,9,13,17,18,20,21,X sur DNAPKi et sur manip A/B/D en OHT et Parpi 
require(tidyverse)
require(plyranges)
require(BSgenome.Hsapiens.UCSC.hg19)
require(TFBSTools)
require(motifmatchr)
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


#Test if genes overlap more often than expected with compD 
require(regioneR)
A <- genes_DIVA %>% as_tibble() %>% 
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% as_granges() %>% filter(Type =="Upregulated",seqnames %in% chr.to.study)
B <- compD.bw %>% as_granges()
pt <- permTest(A=A, ntimes=50, randomize.function=resampleRegions, universe=filter(genes_DIVA,seqnames %in% chr.to.study),
               evaluate.function=numOverlaps, B=B, verbose=F)

#Do the enrichGO analysis
mes_genes.PC1 <- genes_DIVA %>% filter_by_overlaps(as_granges(compD.bw)) %>% 
  as_tibble() %>% 
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname"))

require(clusterProfiler)
formula_res_1 <- enrichGO(mes_genes.PC1$gene_id,
                              OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                              keyType       = 'ENSEMBL',
                              ont           = "BP",
                              universe = genes_DIVA$gene_id,
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.1,
                              qvalueCutoff  = 0.1,
                              readable=TRUE)
p1 <- formula_res_1 %>% dotplot(showCategory=100)
formula_res_1@result %>% as_tibble() %>% write_tsv(glue::glue("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/EnrichGO_on_DE_genes_by_compD_in3Replicates.tsv"))

formula_res_2 <- compareCluster(gene_id~Type, data=mes_genes.PC1, fun="enrichGO",
                              OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                              keyType       = 'ENSEMBL',
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              universe = genes_DIVA$gene_id,
                              pvalueCutoff  = 0.1,
                              qvalueCutoff  = 0.1,
                              readable=TRUE)


p2 <- formula_res_2 %>% dotplot(showCategory=30)
ggsave(glue::glue("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/EnrichGO_on_DE_genes_by_compD_in3Replicates_updown.pdf"),p2,width = 12,height=10)
formula_res_2@compareClusterResult %>% as_tibble() %>% write_tsv(glue::glue("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/EnrichGO_on_DE_genes_by_compD_in3Replicates_updown.tsv"))


#Overlap between transloc and A/B/D genes

data_files_breakpoints <- c(
  "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/ALL-CANCER-INTERCHROM.csv",
  "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/BRCA Cancer-interchrom.csv"
) %>% setNames(str_remove(basename(.),".csv"))%>% map(read_csv)

data_files_breakpoints.GR <- data_files_breakpoints %>% 
  map(dplyr::select,ID,chrA,posA,chrB,posB) %>% 
  map(mutate,chrA=str_c("chr",chrA),chrB=str_c("chr",chrB)) %>% 
  map(pivot_longer,
      chrA:posB,
      c(".value", "Type"),names_pattern = "(chr|pos)([A-B]{1})") %>% 
  map(dplyr::rename,"seqnames"=chr) %>% 
  map(dplyr::rename,"start"=pos) %>% 
  map(dplyr::rename,"name"=ID) %>% 
  map(mutate,width=1) %>% 
  map(drop_na) %>% 
  map(as_granges)


mes_genes.PC1 <- genes_DIVA %>%
  # filter(seqnames %in% chr.to.study) %>% 
  filter_by_overlaps(as_granges(compD.bw)) %>%
  as_tibble() %>% 
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("CompD_{Type}"))

mes_genes.PC1_nocompD <- genes_DIVA %>%
  filter(seqnames %in% chr.to.study) %>%
  filter_by_non_overlaps(as_granges(compD.bw)) %>%
  as_tibble() %>% 
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("NoCompD_{Type}"))


A.list <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD)%>% as_granges() %>% split(.,.$Type)

require(regioneR)

for(ni in names(data_files_breakpoints)){
  my_transloc.GR <- data_files_breakpoints.GR[[ni]]
  out_file <- glue::glue("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PermTest_overlap_between_compD_ABD_genes_and_transloc_{str_replace(ni,' ','_')}.pdf")
  pdf(out_file)
  for(nii in names(A.list)){
    A <- A.list[[nii]]
    title_name <- glue::glue("Genes {nii} ({length(A)})")
    pt <- permTest(A=A, ntimes=1000, randomize.function=resampleRegions, universe=filter(genes_DIVA,seqnames %in% chr.to.study),
                   evaluate.function=numOverlaps, B=my_transloc.GR, verbose=F)
    p <- as.ggplot(expression(plot(pt))) + ggtitle(title_name)
    print(p)
  }
  dev.off()
}


# A <- mes_genes.PC1 %>% as_granges() %>% filter(Type =="Upregulated")
# A.list <- append(mes_genes.PC1_nocompD %>% as_granges() %>% split(.,.$Type),
#                  mes_genes.PC1 %>% as_granges() %>% split(.,.$Type))
# B <- data_files_breakpoints.GR[[1]]
# 
# pt.res <- A.list %>% lapply(function(A){
#   permTest(A=A, ntimes=1000, randomize.function=resampleRegions, universe=filter(genes_DIVA,seqnames %in% chr.to.study),
#            evaluate.function=numOverlaps, B=B, verbose=F)
# })
# 
# # pdf("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PermTest_overlap_between_compD_ABD_genes_and_transloc_ALL-CANCER-INTERHCROM.pdf")
# pdf("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PermTest_overlap_between_compD_ABD_genes_and_transloc_BRCA-CANCER-INTERHCROM.pdf")
# as.ggplot(expression(plot(pt.res[["Upregulated"]]))) + ggtitle("Genes in compD (Upregulated)")
# as.ggplot(expression(plot(pt.res[["Downregulated"]]))) + ggtitle("Genes in compD (Downregulated)")
# as.ggplot(expression(plot(pt.res[["None"]]))) + ggtitle("Genes in compD (None)")
# dev.off()


seqlens <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 %>% seqlengths()
bin_genome <- tileGenome(seqlens[chr.to.study],tilewidth = 1000,cut.last.tile.in.chrom = T)
#ALL cancer
bin_genome <- bin_genome %>% mutate(count_transloc = count_overlaps(bin_genome,B))


pt <- permTest(A=A, randomize.function=resampleRegions, universe=filter(genes_DIVA,seqnames %in% chr.to.study),
               evaluate.function=meanInRegions, x=bin_genome,ntimes=1000)

#test with dnapki
compD.dnapi <- "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHTDNAPKi.bw" %>% import.bw() %>% filter(score>0) %>% 
  filter(seqnames %in% glue::glue("chr{c(1,2,3,6,7,8,9,13,17,18,20,21,'X')}"))
mes_genes.PC1 <- genes_DIVA %>% filter_by_overlaps(compD.dnapi) %>% 
  as_tibble() %>% 
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("CompD_{Type}"))

mes_genes.PC1_nocompD <- genes_DIVA %>%
  filter(seqnames %in% glue::glue("chr{c(1,2,3,6,7,8,9,13,17,18,20,21,'X')}")) %>%
  filter_by_non_overlaps(as_granges(compD.dnapi)) %>%
  as_tibble() %>% 
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("NoCompD_{Type}"))


A.list <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD)%>% as_granges() %>% split(.,.$Type)


for(ni in names(data_files_breakpoints)){
  my_transloc.GR <- data_files_breakpoints.GR[[ni]]
  out_file <- glue::glue("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PermTest_overlap_between_compD_ABD_genes_and_transloc_{str_replace(ni,' ','_')}_DNAPKi.pdf")
  pdf(out_file)
  for(nii in names(A.list)){
    A <- A.list[[nii]]
    title_name <- glue::glue("Genes {nii} ({length(A)})")
    pt <- permTest(A=A, ntimes=1000, randomize.function=resampleRegions, universe=filter(genes_DIVA,seqnames %in% glue::glue("chr{c(1,2,3,6,7,8,9,13,17,18,20,21,'X')}")),
                   evaluate.function=numOverlaps, B=my_transloc.GR, verbose=F)
    p <- as.ggplot(expression(plot(pt))) + ggtitle(title_name)
    print(p)
  }
  dev.off()
}



# 
# A <- mes_genes.PC1 %>% as_granges()
# B <- data_files_breakpoints.GR[[1]]
# pt <- permTest(A=A, ntimes=1000, randomize.function=resampleRegions, universe=filter(genes_DIVA,seqnames %in% glue::glue("chr{c(1,2,3,6,7,8,9,13,17,18,20,21,'X')}")),
#                evaluate.function=numOverlaps, B=B, verbose=F)
# 
# pdf("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PermTest_overlap_between_compD_ABD_genes_and_transloc_ALL-CANCER-INTERHCROM_DNAPKi.pdf")
# plot(pt)
# dev.off()

#motifs
#from recherche_motif_dnapki

#Gene compD vs gene noCompD
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
  mutate(score = rowSums(across(starts_with("score"))))


mes_genes.PC1 <- c(
  genes_DIVA %>% filter_by_overlaps(as_granges(filter(compD.bw,score==3))) %>% mutate(group = "compD"),
  genes_DIVA %>% filter_by_overlaps(as_granges(filter(compD.bw,score==0))) %>% mutate(group = "NocompD")
) %>% 
  as_tibble() %>% 
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname"))


MotifList <- readJASPARMatrix("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
#matchMotifs function
up.distance <- 1000
down.distance <- 0
up.compD <- mes_genes.PC1 %>% filter(group=="compD")
down.compD <- mes_genes.PC1 %>% filter(group=="NocompD")
# down.compD <- cc.dsb %>% filter(Group.PC1=="Everything_else")
genes.up <- mes_genes.PC1 %>% as_granges() %>% 
  promoters(up.distance,down.distance)
# genes.up$Group.PC1 <- c(rep("TopCompD",nrow(up.compD)),rep("Everything_else",nrow(down.compD)))


genes.up.seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19,genes.up)
names(genes.up.seq) <- genes.up$gene_name
# genes.up.seq[genes.up$group =="CompD"] %>% writeXStringSet(glue::glue("../results/DNA_sequences_in_compD_up_genes_{out_name_file}_{up.distance}_{down.distance}.fa"))
# genes.up.seq[!genes.up$group =="CompD"] %>% writeXStringSet(glue::glue("../results/DNA_sequences_in_NOcompD_up_genes_{out_name_file}_{up.distance}_{down.distance}.fa"))
ix2 = motifmatchr::matchMotifs(MotifList,genes.up.seq,out = "scores")

#count number of motifs in each sequence
mat2 = motifmatchr::motifCounts(ix2)
hihi = motifmatchr::motifScores(ix2)
colnames(mat2) <- names(MotifList)
rownames(mat2) <- paste(genes.up)

countMat <- mat2 %>% as.matrix() %>% reshape2::melt()
genes.up.count <- genes.up %>% as_tibble() %>% mutate(pos = paste(genes.up)) %>% right_join(countMat, by = c("pos"="Var1")) %>% dplyr::rename(Motif="Var2",count=value)



#group by Type
genes.up.count.group <- genes.up.count %>% group_by(Type,group,Motif) %>% summarise(count = sum(count))
genes.up.count.group12 <- genes.up.count.group %>% spread(key = group,value = count) %>% nest()
genes.up.count.group12 <- genes.up.count.group12 %>% mutate(enrichment = map(data,function(x){
  total <- colSums(x[,2:3])
  x$stat.fisher.universe <- apply(x[,2:3],1,function(z){fisher.test(matrix(c(z[[1]],z[[2]],total-z),byrow=T,ncol=2),alternative = "greater")$p.value})
  x$adj.fisher.universe <- p.adjust(x$stat.fisher.universe,method = "BH")
  return(x)
}))%>% dplyr::select(-data) %>% unnest(enrichment)

#Ungrouped
genes.up.count.group <- genes.up.count %>% group_by(group,Motif) %>% summarise(count = sum(count))

genes.up.count.group12 <- genes.up.count.group %>% spread(key = group,value = count)
total <- colSums(genes.up.count.group12[,2:3])
genes.up.count.group12$stat.fisher.universe <- apply(genes.up.count.group12[,2:3],1,function(z){fisher.test(matrix(c(z[[1]],z[[2]],total-z),byrow=T,ncol=2),alternative = "greater")$p.value})
genes.up.count.group12$adj.fisher.universe <- p.adjust(genes.up.count.group12$stat.fisher.universe,method = "BH")

#Final table
genes.up.count.final <- genes.up.count.group12 %>%
  group_by(Type) %>%
  # mutate(total_Low = sum(`Everything_else`),total_High = sum(`TopCompD`)) %>% mutate(Enrichment = (`TopCompD`/total_High)/(`Everything_else`/total_Low)) %>%
  mutate(total_Low = sum(`NocompD`),total_High = sum(`compD`)) %>% mutate(Enrichment = (`compD`/total_High)/(`NocompD`/total_Low)) %>%
  filter(adj.fisher.universe < 0.4) %>%
  dplyr::select(-total_Low,-total_High) %>%
  dplyr::select(Motif,Enrichment,everything())
theme_set(theme_classic(base_size=18))
genes.up.count.final %>%
  filter(is.finite(Enrichment)) %>% 
  ggplot(aes(x=fct_reorder(Motif,Enrichment),y=Enrichment,col=-log10(adj.fisher.universe),size=Enrichment)) + geom_point() +coord_flip() +
  facet_wrap(~Type) +
  xlab("Motifs")+scale_color_viridis_c() +
  # ggtitle(glue::glue("Motif enrichment of top100 compD vs top1 compD genes\nupstream: {up.distance}bp and downstream: {down.distance}bp"))-> pmotifs
  ggtitle(glue::glue("Motif enrichment of compD vs nocompD genes\nupstream: {up.distance}bp and downstream: {down.distance}bp"))-> pmotifs
ggsave(glue::glue("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/Motifs_enrichment_in_compDvsnocompD_genes_in3Rep_nb_genes_{nrow(up.compD)}_{nrow(down.compD)}.pdf"),pmotifs,height=10,width=10)


#Expression levels
set.seed(1234)
mes_genes.PC1.exp <- c(
  genes_DIVA %>% filter_by_overlaps(as_granges(filter(compD.bw,score==3))) %>% mutate(group = "compD"),
  genes_DIVA %>% filter_by_overlaps(as_granges(filter(compD.bw,score==0))) %>% mutate(group = "NocompD"),
  genes_DIVA %>% filter_by_non_overlaps(as_granges(compD.bw)) %>% mutate(group = "random") %>% sample(1000)
) %>% 
  as_tibble() %>% 
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname"))
p.boxplot <- mes_genes.PC1.exp %>% ggplot(aes(x=group,y=logFC,fill=group)) + geom_boxplot()
ggsave(glue::glue("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/boxplot_compDvsnocompD_genes_in3Rep_nb_genes_{nrow(up.compD)}_{nrow(down.compD)}.pdf"),p.boxplot,height=6,width=6)


##Recherche of certain motifs

motif_name_list <- c("C12orf32",
                     "H2AFY2",
                     "TP53",
                     "ZNF384",
                     "RAD51AP1",
                     "MTMR15",
                     "RAD17",
                     "RBBP8",
                     "MCM2",
                     "LIG1",
                     "MSH2",
                     "XPC",
                     "SMARCAD1",
                     "TOPBP1",
                     "POLB",
                     "PSMD4",
                     "RBBP4",
                     "COPS2",
                     "WDR3",
                     "SGOL2",
                     "ZNF829",
                     "UTP3",
                     "ZMYND11",
                     "TRMT6",
                     "ZNF281",
                     "PRPF3",
                     "TEX10",
                     "TWISTNB",
                     "RRP1",
                     "PPAN",
                     "CCT5",
                     "CDC73",
                     "GTF3C2",
                     "PES1",
                     "IFIT3",
                     "RFC4",
                     "RAD18",
                     "CBX5",
                     "NUCKS1",
                     "RFC3",
                     "RNF168",
                     "BOD1L",
                     "TCOF1",
                     "UNG",
                     "MED8",
                     "C18orf25",
                     "CHTF18",
                     "RING1",
                     "RAD23B",
                     "HDAC1",
                     "PNKP",
                     "FANCE",
                     "BAP1",
                     "MLH1",
                     "TRIM28",
                     "ERCC5",
                     "CLSPN",
                     "DAXX",
                     "ATF2",
                     "NPM1",
                     "TIPIN",
                     "FANCD2",
                     "XRCC2",
                     "MCM10",
                     "C9orf80",
                     "POLR2F",
                     "GLRX2",
                     "AHNAK",
                     "USP39",
                     "FBXW7",
                     "ZBTB10",
                     "DTL",
                     "UTP14A",
                     "HTATSF1",
                     "DDX21",
                     "ZC3H11A",
                     "HSPA4",
                     "BCL7C",
                     "BUD31",
                     "Prpf40a",
                     "USP6",
                     "UBE2R2",
                     "CPSF4",
                     "RNF138",
                     "PGAM2",
                     "NSA2",
                     "BCL7A",
                     "HP1BP3",
                     "BATF3",
                     "LIMCH1",
                     "ZC3H8",
                     "PCGF4",
                     "DFFA",
                     "RNF113A",
                     "FLJ23518",
                     "ZRANB2",
                     "MAPK14",
                     "MAGEB3",
                     "AFF4",
                     "PTENP1",
                     "FOS",
                     "RDBP",
                     "DPF2",
                     "POU2F2",
                     "SF3B3",
                     "PIK3R2",
                     "AHDC1",
                     "GTF2F1",
                     "CENPM",
                     "PRPF31",
                     "ZSCAN5A",
                     "MTF2",
                     "TUT1",
                     "ZMAT1",
                     "THOC5",
                     "PPM1G",
                     "PTPN6",
                     "UBA2",
                     "RACGAP1",
                     "WBP11",
                     "TCEA1",
                     "TAF15",
                     "PRPF6",
                     "CPSF2",
                     "PRPF38A",
                     "DNTTIP2",
                     "NOP56",
                     "GATAD1",
                     "MKI67IP",
                     "FBL",
                     "RRS1",
                     "GLTSCR2",
                     "XPO5",
                     "WDR46",
                     "DMAP1",
                     "RRP9",
                     "ISG20L2",
                     "CMAS",
                     "MAD1L1",
                     "TAF6",
                     "IMP4",
                     "TXN",
                     "RRP15",
                     "ERI1",
                     "PSMD11",
                     "PHF10",
                     "TRIM21",
                     "ATP6V1C1",
                     "SAE1",
                     "C18orf25",
                     "CENPJ",
                     "UHRF2",
                     "HNRNPK",
                     "EWSR1",
                     "NFIA",
                     "ATXN3",
                     "WEE1",
                     "WHSC1L1",
                     "CCDC99",
                     "MCM4",
                     "UBA1",
                     "SART1",
                     "C19orf6",
                     "CCNF",
                     "ITGB7",
                     "OTUD7B",
                     "SLC4A1AP",
                     "WHSC1",
                     "ORC3",
                     "PBX2",
                     "HHLA2",
                     "YY1",
                     "ATF7",
                     "CCNA2",
                     "MYBL2",
                     "PTEN",
                     "PSMB6",
                     "RASAL3",
                     "TTN")

my_motifs <- MotifList[grepl(str_c(motif_name_list,collapse="|"),names(MotifList))]

up.distance <- 1000
down.distance <- 0
genes.up <- genes_DIVA %>% 
  promoters(up.distance,down.distance)

genes.up.seq <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19,genes.up)
names(genes.up.seq) <- genes.up$gene_name
ix2 = motifmatchr::matchMotifs(my_motifs,genes.up.seq,out = "scores")

#count number of motifs in each sequence
mat2 = motifmatchr::motifCounts(ix2)
colnames(mat2) <- names(my_motifs)
rownames(mat2) <- paste(genes.up)

countMat <- mat2 %>% rowSums() %>% enframe()
genes.up.count <- genes.up %>% as_tibble() %>% mutate(pos = paste(genes.up)) %>% right_join(countMat, by = c("pos"="name")) %>% dplyr::rename(count=value) %>% 
  # filter(seqnames %in% glue::glue("chr{c(1,2,3,6,7,8,9,13,17,18,20,21,'X')}"))
  filter(seqnames %in% glue::glue("chr{c(1,17)}"))

genes.up.count <- genes.up.count %>% mutate(Cat = ifelse(count>20,"TF","noTF")) %>% split(.,.$Cat)

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



cc.cov <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipA.bw" %>% import.bw(as="RleList")
cc.dsb <- genes.up.count %>%
  map(as_granges) %>% 
  map(function(one_dsb){
    one_dsb$name <- one_dsb$gene_id
    Get1val(Name = "cc",one.w = cc.cov,x = one_dsb)
  }) %>% bind_rows(.id = "Name")


p1 <- cc.dsb %>% ggplot(aes(x=Name,y=value,fill=Name)) + geom_boxplot()


p2 <- genes.up.count %>% bind_rows() %>% left_join(dplyr::select(cc.dsb,rowname,value),by = c("gene_id"="rowname")) %>% 
  ggplot(aes(x=count,y=value)) + geom_point()
pdf("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/TF_list_compD_value.pdf",width=8)
print(p1)
print(p2)
dev.off()
