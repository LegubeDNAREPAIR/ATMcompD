# Les DSB switch elle du comp A a du comp B ? sur les 80 pour la composante a 100kb
# 
# Dans les gènes du D compartment
# Les gènes de la DDR up et down-régulés 
# Gènes de la DDR up-régulés mais qui ne vont pas dans le D, sur le même chromosome 
# Correlation distance DSB et niveau de compD
require(tidyverse)
require(plyranges)

bless80 <- read_bed("/mnt/NAS/DATA/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed") %>% sortSeqlevels() %>% sort()


reskb <- "100kb"
AB_GR <- c(list.files("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/AB_comp/",pattern=reskb,full.names = T),
           list.files("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC/compartmentAB/",pattern=reskb,full.names = T)) %>%
  setNames(basename(.)) %>% 
  map(import.bedGraph) %>%  
  map(function(x){
    x %>% mutate(score = ifelse(seqnames =="chr22",score*-1,score)) # reverse signal on chromosome 22 for every compartment
  })


res_AB <- lapply(AB_GR,function(one.GR){
  bless80 %>% plyranges::select(-score)%>% join_overlap_inner(one.GR) %>% as_tibble() %>% dplyr::select(name,score)
}) %>% bind_rows(.id="Condition") %>% mutate(Condition = str_remove(str_remove(Condition,"ABcomp_|compartmentAB_"),glue::glue("_{reskb}.bedGraph"))) %>% spread(key=Condition,value=score)
res_AB %>% write_tsv("../../results/compAB_values_80DSB.tsv")


#Classer la matrice du clustering par rapport a leur niveau de compartiment celles qui clusters sont celles dans le A


# saddle plot + heatmap françois + le tableau reverse le chromosome 22 peu importe le type de DSB