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
    # logFC < 0 & FILTER.FC == 1 ~ "Downregulated",
    logFC < -0.3 ~ "Downregulated",
    # logFC > 0 & FILTER.FC == 1 & FILTER.P == 1~ "Upregulated",
    # logFC > 0 & FILTER.FC == 1~ "Upregulated",
    logFC > 0.3 ~ "Upregulated",
    TRUE ~ "None"
  )) 

ens.genes <- read.table("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed",sep="\t",h=T) %>% GRanges()
ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))
mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
  arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
# gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 1000000)
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
  map(filter,seqnames %in% chr.to.study) %>% 
  map(as_tibble) %>% 
  map(dplyr::select,-strand,-width) %>% 
  purrr::reduce(left_join,by=c("seqnames","start","end")) %>% 
  mutate_at(vars(starts_with("score")),function(x){
    as.numeric(x>0)
  }) %>% 
  mutate(score = rowSums(across(starts_with("score")))) 
# %>% filter(score == 3)


#Test if genes overlap more often than expected with compD 

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
  map(as_granges) %>% filter(score ==3)


mes_genes.PC1 <- genes_DIVA %>%
  filter(seqnames %in% chr.to.study) %>%
  filter_by_overlaps(as_granges(compD.bw)) %>%
  as_tibble() %>%
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("CompD_{Type}"))

# mes_genes.PC1 <- genes_DIVA %>%
#   filter(seqnames %in% chr.to.study) %>%
#   filter_by_overlaps(as_granges(compD.bw)%>% filter(score == 3)) %>%
#   as_tibble() %>% 
#   left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("CompD_{Type}"))

mes_genes.PC1_nocompD <- genes_DIVA %>%
  filter(seqnames %in% chr.to.study) %>%
  filter_by_non_overlaps(as_granges(compD.bw)) %>%
  as_tibble() %>%
  left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("NoCompD_{Type}"))
# mes_genes.PC1_nocompD <- genes_DIVA %>%
#   filter(seqnames %in% chr.to.study) %>%
#   filter_by_overlaps(as_granges(compD.bw) %>% filter(score < 3)) %>%
#   as_tibble() %>%
#   left_join(DE_DIVA_INFO,by = c("gene_id"="rowname")) %>% mutate(Type = glue::glue("NoCompD_{Type}"))


A.list <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD)%>% as_granges() %>% split(.,.$Type)

require(regioneR)

# for(ni in names(data_files_breakpoints)){
#   my_transloc.GR <- data_files_breakpoints.GR[[ni]]
#   out_file <- glue::glue("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PermTest_overlap_between_compD_ABD_genes_and_transloc_{str_replace(ni,' ','_')}.pdf")
#   pdf(out_file)
#   for(nii in names(A.list)){
#     A <- A.list[[nii]]
#     title_name <- glue::glue("Genes {nii} ({length(A)})")
#     pt <- permTest(A=A, ntimes=1000, randomize.function=resampleRegions, universe=filter(genes_DIVA,seqnames %in% chr.to.study),
#                    evaluate.function=numOverlaps, B=my_transloc.GR, verbose=F)
#     p <- as.ggplot(expression(plot(pt))) + ggtitle(title_name)
#     print(p)
#   }
#   dev.off()
# }
# 
# 
# A.list.compD.nocompD <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD) %>% separate(Type,into=c("compD","Type")) %>% 
#   group_by(compD) %>% nest()
# my_transloc.GR <- data_files_breakpoints.GR[[1]]
# 
# A.list.compD.nocompD.res <- A.list.compD.nocompD %>% mutate(permtest = map(data,
#                                                                            function(x){
#                                                                              A <-as_granges(x)
#                                                                              set.seed(1234)
#                                                                              pt <- permTest(A=A, ntimes=1000, randomize.function=resampleRegions, universe=filter(genes_DIVA,seqnames %in% glue::glue("chr{c(1,17,'X')}")),
#                                                                                             evaluate.function=numOverlaps, B=my_transloc.GR, verbose=F)
#                                                                              
#                                                                              tibble(
#                                                                                n = length(A),
#                                                                                permuted = mean(pt$numOverlaps$permuted),
#                                                                                observed = pt$numOverlaps$observed,
#                                                                                pval = format.pval(pt$numOverlaps$pval,2)
#                                                                              )
#                                                                              
#                                                                            })) %>% dplyr::select(-data) %>% unnest()
# my.p.bar <- A.list.compD.nocompD.res %>% 
#   mutate(observed = observed/permuted) %>% 
#   mutate(permuted = permuted/permuted) %>% 
#   gather(key=Type,value = value,-compD,-pval,-n) %>% 
#   mutate(mycolor = case_when(
#     Type == "permuted" ~ "#7f8c8d",
#     Type == "observed" & value > 1 ~ "#27ae60",
#     Type == "observed" & value < 1 ~ "#27ae60"
#   )) %>% 
#   mutate(compD = glue::glue("{compD} ({n})")) %>% 
#   ggplot(aes(x=compD,y=value,fill=mycolor)) + geom_bar(stat="identity",position = "dodge",col="black") + theme_classic() + theme(axis.text.x = element_markdown()) +
#   scale_fill_identity() + scale_y_continuous(labels = scales::percent,limits = c(0,1.2))
# ggsave("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/barplot_permutatio_observed_permuted_percentage_manipABD_compD_nocompD.pdf",my.p.bar,width=4,height=4)
# 


A.list.compD.nocompD <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD) %>% separate(Type,into=c("compD","Type")) %>% 
  group_by(compD,Type) %>% nest()

my_transloc.GR <- data_files_breakpoints.GR[[1]]

A.list.compD.nocompD.res <- A.list.compD.nocompD %>% mutate(permtest = map(data,
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
ggsave("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/barplot_permutatio_observed_permuted_percentage_manipABD_compD_nocompD_bygenetype.pdf",my.p.bar,width=8,height=4)
A.list.compD.nocompD.res %>% write_tsv("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/barplot_permutatio_observed_permuted_percentage_manipABD_compD_nocompD_bygenetype.tsv")
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



A.list.compD.nocompD <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD) %>% separate(Type,into=c("compD","Type")) %>% 
  group_by(compD) %>% nest()

my_transloc.GR <- data_files_breakpoints.GR[[1]]

A.list.compD.nocompD.res <- A.list.compD.nocompD %>% mutate(permtest = map(data,
                                               function(x){
                                                 A <-as_granges(x)
                                                 set.seed(1234)
                                                 pt <- permTest(A=A, ntimes=1000, randomize.function=resampleRegions, universe=filter(genes_DIVA,seqnames %in% glue::glue("chr{c(1,2,3,6,7,8,9,13,17,18,20,21,'X')}")),
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
  gather(key=Type,value = value,-compD,-pval,-n) %>% 
  mutate(mycolor = case_when(
    Type == "permuted" ~ "#7f8c8d",
    Type == "observed" & value > 1 ~ "#27ae60",
    Type == "observed" & value < 1 ~ "#27ae60"
  )) %>% 
  mutate(compD = glue::glue("{compD} ({n})")) %>% 
  ggplot(aes(x=compD,y=value,fill=mycolor)) + geom_bar(stat="identity",position = "dodge",col="black") + theme_classic() + theme(axis.text.x = element_markdown()) +
    scale_fill_identity() + scale_y_continuous(labels = scales::percent,limits = c(0,1.2))
ggsave("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/barplot_permutatio_observed_permuted_percentage_DNAPKi_compD_nocompD.pdf",my.p.bar,width=4,height=4)

#same but with upregulated toossa

A.list.compD.nocompD <- rbind(mes_genes.PC1,mes_genes.PC1_nocompD) %>% separate(Type,into=c("compD","Type")) %>% 
  group_by(compD,Type) %>% nest()

my_transloc.GR <- data_files_breakpoints.GR[[1]]

A.list.compD.nocompD.res <- A.list.compD.nocompD %>% mutate(permtest = map(data,
                                                                           function(x){
                                                                             A <-as_granges(x)
                                                                             pt <- permTest(A=A, ntimes=1000, randomize.function=resampleRegions, universe=filter(genes_DIVA,seqnames %in% glue::glue("chr{c(1,2,3,6,7,8,9,13,17,18,20,21,'X')}")),
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
  scale_fill_identity() + scale_y_continuous(labels = scales::percent,limits = c(0,1.5)) + facet_wrap(~compD,scales="free_x")
ggsave("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/barplot_permutatio_observed_permuted_percentage_DNAPKi_compD_nocompD_bygenetype.pdf",my.p.bar,width=8,height=4)
