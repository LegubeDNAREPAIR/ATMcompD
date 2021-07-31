require(tidyverse)
require(ranger)
require(corrr)
require(plyranges)
require(rtracklayer)
library(FactoMineR)
# Correlation between PCs and Histone marks
PCs.ratio.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHTDNAPKi.bw"
# PCs.ratio.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio.bw"
# PCs.ratio.path <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_HiC_D_OHT.bw"
PCs.ratio <- import.bw(PCs.ratio.path)
cc.cov <- PC.path %>% import.bw(as="RleList")
out_name_file <- basename(PCs.ratio.path) %>% str_remove(".bw")
PCs.ratio$name <- str_c("bin",1:length(PCs.ratio))


mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
  arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 1000000)

if(length(seqlevels(PCs.ratio))>4){
  chr.to.study <- glue::glue("chr{c(2,6,9,13,18,1,17,20,'X')}")
}else{
  chr.to.study <- seqlevels(PCs.ratio)
}
PCs.ratio.nogamma <- PCs.ratio %>% filter_by_non_overlaps(gamma_region) %>% filter(seqnames %in% chr.to.study)


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


data_files_breakpoints.GR[["gamma_region"]] <- gamma_region %>% anchor_center() %>% mutate(width=100000)

data_files_breakpoints.GR[["random_region"]] <- tileGenome(seqlens[1:23],tilewidth = 100000,cut.last.tile.in.chrom = T) %>% 
  filter_by_non_overlaps(Reduce("c",data_files_breakpoints.GR)) %>% sample(1000) %>% mutate(name = glue::glue("Rdm{1:length(.)}"))

cc.dsb <- data_files_breakpoints.GR%>% 
  map(function(one_dsb){
    PhDfunc::Get1val(Name = "cc",one.w = cc.cov[chr.to.study],x = filter(one_dsb,seqnames %in% chr.to.study))
  }) %>% bind_rows(.id = "Type")


cc.dsb %>% ggplot(aes(x=Type,y=value,fill=wig)) + geom_boxplot() + theme(legend.position="none")



# create bigwig


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

sizebin <- 500000
sizebintxt <- glue::glue("{sizebin/1e03}kb")
seqlens <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 %>% seqlengths()
bin_genome <- tileGenome(seqlens[1:23],tilewidth = sizebin,cut.last.tile.in.chrom = T)


#ALL cancer
bin_genome <- bin_genome %>% mutate(count_transloc = count_overlaps(bin_genome,data_files_breakpoints.GR[[1]]))
bin_genome_bw <- coverage(bin_genome,weight = "count_transloc")
bin_genome_noDSB <- bin_genome %>% filter_by_non_overlaps(gamma_region)
bin_genome_bw_noDSB <- coverage(bin_genome_noDSB,weight = "count_transloc")


bin_genome_bw_noDSB %>% export.bw(glue::glue("../../results/ALL-CANCER-INTERCHROM_noDSB200_{sizebintxt}.bw"))
bin_genome_bw %>% export.bw(glue::glue("../../results/ALL-CANCER-INTERCHROM_{sizebintxt}.bw"))
# BRCA cancer
bin_genome <- bin_genome %>% mutate(count_transloc = count_overlaps(bin_genome,data_files_breakpoints.GR[[2]]))
bin_genome_bw <- coverage(bin_genome,weight = "count_transloc")
bin_genome_noDSB <- bin_genome %>% filter_by_non_overlaps(gamma_region)
bin_genome_bw_noDSB <- coverage(bin_genome_noDSB,weight = "count_transloc")


bin_genome_bw_noDSB %>% export.bw(glue::glue("../../results/BRCA-CANCER-INTERCHROM_noDSB200_{sizebintxt}.bw"))
bin_genome_bw %>% export.bw(glue::glue("../../results/BRCA-CANCER-INTERCHROM_{sizebintxt}.bw"))




