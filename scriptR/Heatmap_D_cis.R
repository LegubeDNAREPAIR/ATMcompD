#From /home/rochevin/Documents/PROJET_INGE/HiC_Coline/REVISIONS/ORIGINAL_SCRIPTS/Heatmap_HiC_data_DSB.R
require(HiTC)
require(tidyverse)
require(plyranges)
require(reshape2)
require(keras)
require(rtracklayer)
require(cowplot)
require(BSgenome.Hsapiens.UCSC.hg19.masked)
require(regioneR)
require(ggtext)
seqlens <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 %>% seqlengths()
#Functions

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

get_trans_matrix_Hic <- function(file,bed,binsize){
  c(my.ranges,intdata) %<-% file
  
  val1 <- intdata %>% sum
  
  overlaps <- as(findOverlaps(bed,my.ranges), "List")
  overlaps.bin <- lapply(overlaps,length) %>% unlist()
  overlaps <- overlaps[overlaps.bin == binsize]
  
  res <- lapply(1:length(overlaps),function(i){
    x <- overlaps[[i]]
    if(as.character(strand(bed[i])) == "-"){
      x <- rev(x)
    }
    intdata[x,x] %>% round()
  })
  return(list(res,val1))
}
get_sub_matrix_Hic <- function(file,bed,binsize){
  HTC <- loadRData(file)
  rows.chr <- x_intervals(HTC)
  val1 <- HTC %>% intdata %>% sum
  # HTC <- normPerReads(HTC)
  
  overlaps <- as(findOverlaps(bed,rows.chr), "List")
  overlaps.bin <- lapply(overlaps,length) %>% unlist()
  overlaps <- overlaps[overlaps.bin == binsize]
  
  res <- lapply(1:length(overlaps),function(i){
    x <- overlaps[[i]]
    if(as.character(strand(bed[i])) == "-"){
      x <- rev(x)
    }
    intdata(HTC)[x,x]
  })
  return(list(res,val1))
}
get_nb_read_by_chr <- function(file){
  HTC <- loadRData(file)
  HTC %>% intdata %>% sum
}


plot_telquel <-function(data.plot,my.quantile = 0.90,facet=TRUE,window = 1000000,fixed = F,fixed.limite = 1000,DSB = "DSB"){
  if(fixed ==T){
    limite <- fixed.limite
    
  }else{
    limite <- data.plot %>% pull(value) %>% quantile(my.quantile) %>% as.numeric()
  }
  
  debut <- data.plot$Var1 %>% levels %>% as.numeric() %>% min()
  milieu <- data.plot$Var1 %>% levels %>% as.numeric() %>% median()
  fin <- data.plot$Var1 %>% levels %>% as.numeric() %>% max()
  p.ratio <- data.plot %>%
    mutate(value = ifelse(value > limite,limite,value)) %>%
    # mutate(value = log10(value)) %>%
    ggplot(aes(x=Var1,y=Var2,fill=value)) + geom_tile() +
    scale_fill_gradient2(low = "white",high = "#EA2027",limits = c(0,limite)) +
    scale_x_discrete(name = 'Position',
                     breaks = c(debut,milieu,fin),
                     labels = c(-window,DSB,window)
    ) +
    scale_y_discrete(name = 'Position',
                     breaks = c(debut,milieu,fin),
                     labels = c(-window,DSB,window)
    ) +
    # geom_vline(xintercept = milieu,linetype ="dashed",col="red")+
    theme_classic() +
    theme(legend.position="bottom")
  if(facet == TRUE){
    p.ratio + facet_wrap(~Condition,ncol=2)
  }else{
    p.ratio
  }
}

#ADD binned 5kb genomic region with no damage and lot of damage (excluding damages)
# ens.genes <- read.table("/home/rochevin/Documents/PROJET_THESE/DSB_EFFECT_ON_GENE/data/EnsDB.Hsapiens.v75.genes.bed",sep="\t",h=T) %>% GRanges()
# ens.genes <- regioneR::filterChromosomes(ens.genes,keep.chr=c(1:22,"X","Y"))
# seqlevels(ens.genes) <- paste0("chr",seqlevels(ens.genes))
# ens.genes.ext <- ens.genes %>% anchor_5p() %>% stretch(3000)
# prom.genes <- ens.genes.ext %>% promoters(1,1)
# prom.genes$start_tss <- start(prom.genes)
DSB174 <- read_bed("/mnt/NAS/DATA/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")
nb_of_damage <- DSB174 %>% as_tibble() %>% dplyr::count(seqnames)
most_damage_chr <- nb_of_damage %>% arrange(desc(n)) %>% dplyr::slice(1:2) %>% pull(seqnames)
less_damage_chr <- nb_of_damage %>% arrange(n) %>% dplyr::slice(1) %>% pull(seqnames)

extend.DSB <- DSB174 %>% anchor_center() %>% mutate(width=4000000)
MaskedRegions <- getMask(BSgenome.Hsapiens.UCSC.hg19.masked)

RdmTiles <- . %>% tileGenome(tilewidth = 5000,cut.last.tile.in.chrom=T) %>%
  filter_by_non_overlaps(extend.DSB) %>%
  filter_by_non_overlaps(MaskedRegions) %>%
  sample(5000) %>% anchor_center() %>% mutate(width=1)
# less_damage_bins <- seqlens[less_damage_chr] %>% RdmTiles
# most_damage_bins <- seqlens[most_damage_chr] %>% RdmTiles
# less_damage_bins_sec <- seqlens["chr13"] %>% RdmTiles
# most_damage_bins_sec <- seqlens["chr17"] %>% RdmTiles
require(regioneR)
RdmRegioneR <- . %>% tileGenome(tilewidth = 5000,cut.last.tile.in.chrom=T) %>%
  sample(5000) %>%
  randomizeRegions(mask = c(extend.DSB,MaskedRegions),allow.overlaps=F, per.chromosome=T)

RdmRegioneR20kb <- . %>% tileGenome(tilewidth = 5000,cut.last.tile.in.chrom=T) %>%
  sample(20000) %>%
  randomizeRegions(mask = c(extend.DSB,MaskedRegions),allow.overlaps=F, per.chromosome=T)

set.seed(174)
less_damage_bins <- seqlens[less_damage_chr] %>% RdmRegioneR
set.seed(174)
most_damage_bins <- seqlens[most_damage_chr] %>% RdmRegioneR
set.seed(174)
less_damage_bins_sec <- seqlens["chr13"] %>% RdmRegioneR
set.seed(174)
most_damage_bins_sec <- seqlens["chr1"] %>% RdmRegioneR
set.seed(174)
full_chr_random_bins <- seqlens[seqlevels(read_bed("/mnt/NAS/DATA/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed"))] %>% RdmRegioneR20kb

##USING DSB
extend.DSB <- DSB174 %>% anchor_center() %>% mutate(width=4000000)
# My80TADs <- "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC/TAD/TopDom/TAD_topdom_DIvA_manipA_hg19_50kb_ws10.bed" %>% read_bed() %>% 
#     filter_by_non_overlaps(extend.DSB) %>% .[sample(1:length(.),80)] %>% anchor_end() %>% mutate(width=1)
MES_DSB <- list(
  # "DSB174" = read_bed("/mnt/NAS/DATA/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed"),
  # "HR" = read_bed("/mnt/NAS/DATA/AsiSI/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed"),
  # "NHEJ" = read_bed("/mnt/NAS/DATA/AsiSI/BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed"),
  # "RANDOM" = read_bed("/mnt/NAS/DATA/AsiSI/30random.bed"),
  # # "most_damage_bins" = regioneR::filterChromosomes(most_damage_bins,keep.chr=paste0("chr",c(1:22,"X"))),
  # # ,"less_damage_bins" = less_damage_bins
  "most_damage_bins_chr1" = regioneR::filterChromosomes(most_damage_bins_sec,keep.chr=paste0("chr",c(1:22,"X"))),
  "less_damage_bins_chr13" = regioneR::filterChromosomes(less_damage_bins_sec,keep.chr=paste0("chr",c(1:22,"X"))),
  # "full_chr_random_bins" = regioneR::filterChromosomes(full_chr_random_bins,keep.chr=paste0("chr",c(1:22,"X"))),
  "80bless" = read_bed("/mnt/NAS/DATA/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
  # "My80TADs"=read_bed("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/REVISIONS/My80RandomTAds_undamaged_50kb_ws10.bed")
  # "DOM" = read_bed("/mnt/NAS/DATA/AsiSI/100bestCleavedSites_hg19.bed")
  # "DOM" = read_bed( "/home/rochevin/Documents/PROJET_INGE/4CSeq_PROCESS_VINCENT/AFTER_EVA_PROCESSING/PROCESSED_21_01_19/201901_Vincent_C-HiC_subDomains/scripts_Marion/DOMAINES2Mbp_hg19.bed" )
)


mes_windows <- c(
  # "100kb"=100000,
  "50kb"=50000
  # "10kb"=10000,
  # "5kb"=5000
)

mes_limites <- list(
  # "100kb"=c(10000000),
  "50kb"=c(5000000)
  # "10kb"=c(1000000),
  # "5kb"=c(500000,250000)
)





ctrlCOND <- "DIvA"
Replace_DiVA <- FALSE

theme_set(theme_classic(base_size=18))
theme_update(legend.position="top",axis.title = element_blank(),strip.background = element_blank(),axis.line = element_blank()) 



for(obserOE in c("observed")){
  if(obserOE == "OE"){
    my.lim <- 0.05
  }else{
    # my.lim <- ifelse(Replace_DiVA,0.1,0.3)
    my.lim <- 0.1
  }
  for(DSB.n in names(MES_DSB)){
    message(DSB.n)
    DSB174 <- MES_DSB[[DSB.n]]
    
    
    
    for(m.w in names(mes_windows)){
      message(m.w)
      sw <- mes_windows[[m.w]]
      
      if(obserOE == "observed"){
        HiC.files <- list.files("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/HTC/observed",pattern=str_c("_",m.w),full.names=T)
        # HiC.files <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/REVISIONS/matrix_generation/test/observed",pattern=m.w,full.names=T)
      }else{
        HiC.files <- list.files("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/HTC/OE",pattern=str_c("_",m.w),full.names=T)
      }
      
      
      
      for(window in mes_limites[[m.w]]){
        message(window)
        bin <- (window/sw)+1
        
        res.DSB <- mclapply(unique(seqnames(DSB174)),function(chrom){
          message(chrom)
          f <- HiC.files[str_detect(HiC.files,str_c("_",chrom,"_"))]
          if(ctrlCOND == "DIvA" & Replace_DiVA){
            New_DiVa <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/HTC",pattern=m.w,full.names=T) %>% 
              str_subset(glue::glue("DIvA_manipA_{chrom}_"))
            f[str_detect(f,ctrlCOND)] <- New_DiVa
          }
          Replicate <- f %>% map(.%>% basename() %>% str_remove_all("HTC_observed_KR_|_manip[AB]|HTC_observed_KR_HiC_D_|_chr[0-9A-Z]+_[0-9]+kb.RData")) %>% unlist()
          names(f) <- Replicate
          # Replicate <- f %>% map(.%>% basename() %>% str_extract("manip[A-Z]+|rep[0-9]+"))
          diva.chr <- DSB174 %>%filter(seqnames ==chrom) %>% 
            anchor_center() %>% mutate(width = window)
          if(length(diva.chr)<2)
            return(NULL)
          
          f_ctrl <- Replicate %>% str_subset(glue::glue("{ctrlCOND}$"))
          f_exp <- Replicate %>% str_subset(glue::glue("{ctrlCOND}$"),negate = T)
          
          
          c(res_CTRL,val_CTRL) %<-% get_sub_matrix_Hic(file = f[f_ctrl],bed = diva.chr,binsize = bin)
          res_CTRL <- Reduce("+",res_CTRL)
          res <- lapply(f[f_exp],function(s){
            c(sres_exp,sval_exp) %<-% get_sub_matrix_Hic(file = s,bed = diva.chr,binsize = bin)
            
            sres_exp <- lapply(sres_exp,function(x){
              x * (val_CTRL/sval_exp)
            })
            Reduce("+",sres_exp)
          })
          
          
          
          res[["CTRL"]] <- res_CTRL
          
          res
          
          
          
        },mc.cores=12)
        Replicate <- res.DSB %>% map(names) %>% unlist() %>% unique()
        data.plot <- lapply(Replicate,function(one_exp){
          exp_count <- lapply(res.DSB,function(x){
            x[[one_exp]]
          }) %>% plyr::compact()
          Reduce("+",exp_count) %>% as.matrix %>% melt
        }) %>% setNames(Replicate) %>% bind_rows(.id = "Condition")
        
        data.plot <- data.plot %>% 
          mutate(Var1 = forcats::lvls_revalue(Var1, as.character(1:bin))) %>%
          mutate(Var2 = forcats::lvls_revalue(Var2, as.character(1:bin)))
        
        
        data.ratio <- data.plot %>%
          spread(key = Condition,value = value)%>%
          gather(key = "Condition",value = value,-Var1,-Var2,-CTRL) %>% 
          # mutate(OHT = OHT * (487/393)) %>% 
          mutate(ratio = log2(value/CTRL)) %>%
          mutate(nlratio = (value/CTRL)) %>% 
          mutate(Condition = glue::glue("{Condition}/{ctrlCOND}"))
        
        
        
        debut <- data.ratio$Var1 %>% levels %>% as.numeric() %>% min()
        milieu <- data.ratio$Var1 %>% levels %>% as.numeric() %>% median()
        fin <- data.ratio$Var1 %>% levels %>% as.numeric() %>% max()
        nb_annot <- scales::number(window/2, scale = 1e-6,accuracy=0.01,suffix = "Mb")
        p.ratio <- data.ratio %>%
          mutate(ratio = ifelse(ratio > my.lim,my.lim,ratio)) %>%
          mutate(ratio = ifelse(ratio < -my.lim,-my.lim,ratio)) %>%
          ggplot(aes(x=Var1,y=Var2,fill=ratio)) + geom_tile() +
          scale_fill_gradient2(low = "#2980b9",high = "#f1c40f",midpoint = 0,mid = "black",limits = c(-my.lim,my.lim)) +
          scale_x_discrete(breaks = c(debut,milieu,fin),
                           labels = c(str_c("-",nb_annot),"0",nb_annot)
          ) +
          scale_y_discrete(breaks = c(debut,milieu,fin),
                           labels = c(str_c("-",nb_annot),"0",nb_annot)
          )  +
          facet_wrap(~Condition,scales="free") +
          guides(fill = guide_colourbar(title.position = "top",title = "log2(Fold Change)",title.hjust=0.5,barwidth =unit(20,"lines"),barheight =unit(1,"lines")))
        my_out_file <- ifelse(Replace_DiVA,
                              glue::glue("{DSB.n}_{m.w}_{window}_{obserOE}_Hich_average_ratio_with_OLD_DIVA.pdf"),
                              ifelse(ctrlCOND == "DIvA",glue::glue("{DSB.n}_{m.w}_{window}_{obserOE}_Hich_average_ratio.pdf"),
                                     glue::glue("{DSB.n}_{m.w}_{window}_{obserOE}_Hich_average_ratio_{ctrlCOND}.pdf")))
        pdf(my_out_file,width=12,height=12)
        print(p.ratio)
        dev.off()
        
      }
      
      
    }
  }
}
