require(BSgenome.Hsapiens.UCSC.hg19)
require(diffHic)
require(edgeR)
require(tidyverse)
require(plyranges)
require(data.table)
# require(normOffsets)
require(ggside)
require(Matrix)
require(patchwork)


GetIndex <- function(df,mes_chrs,indices_ranges){
  
  row.regions <- tibble(seqnames = mes_chrs[1],start=df[,1]+1) %>% 
    left_join(indices_ranges,by = c("seqnames","start")) %>% pull(idx)# interaction start
  col.regions <- tibble(seqnames = mes_chrs[2],start=df[,2]+1) %>% 
    left_join(indices_ranges,by = c("seqnames","start")) %>% pull(idx)# interaction end
  return(tibble(row.regions,col.regions))
}
process_HiC <- function(list_files,indices_ranges=NULL,filters=NULL,filters2=NULL,binSize=1000000){
  
  mes_chrs <- str_extract_all(list_files[1],"chr[A-Z0-9]+")[[1]]
  message(glue::glue("{mes_chrs[1]}:{mes_chrs[2]}"))
  chrom_length <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[mes_chrs]
  
  # first_pos <- indices_ranges %>% filter(seqnames == mes_chrs[[1]]) %>% pull(idx)
  # second_pos <- indices_ranges %>% filter(seqnames == mes_chrs[[2]]) %>% pull(idx)
  res_hic <- list_files %>% map(function(x){
    Condition <- x %>% basename() %>% str_remove(filetypename) %>% str_remove(glue::glue("_chr[A-Z0-9]+_chr[A-Z0-9]+_{reskb}.txt.gz"))
    one_dat <- as.matrix(fread(cmd=glue::glue("zcat {x}"),sep='\t',header=F))
    one_dat[,3]=round(one_dat[,3],2)
    max_cols <- colMaxs(one_dat[,-3]) %>% which.max()
    if(max_cols!=which.max(chrom_length)){
      mes_chrs <- rev(mes_chrs)
      # chrom_length <- rev(chrom_length)
    }
    
    # binendi=ceiling(chrom_length/binSize)
    # one_dat=rbind(one_dat,c((binendi-1)*binSize,0))
    
    hic.gi <- GetIndex(one_dat,mes_chrs=mes_chrs,indices_ranges=indices_ranges)%>% mutate(value = one_dat[,3]) %>% mutate(cond = Condition)#replace position with index on my.ranges
    # data.Mati=sparseMatrix(i=(one_dat[,1]/binSize)+1,j=(one_dat[,2]/binSize)+1,x=one_dat[,3])
    hic.gi <- hic.gi %>% drop_na()
    
    
    return(hic.gi)
  }) %>% bind_rows() 
  
  res_hic <- res_hic %>% pivot_wider(names_from = cond, values_from = value)
  
  
  names_cond <- res_hic %>% dplyr::select(-1:-2) %>% colnames()
  names_cond <- rep(0,length(names_cond)) %>% setNames(names_cond)%>% as.list()
  
  res_hic <- res_hic %>% replace_na(names_cond)
  # %>% purrr::reduce(full_join,by=c("row.regions","col.regions")) %>% drop_na()
  
  
  
  
  #Create group of interactions using index in col and row, from myranges
  gi <- GInteractions(res_hic$row.regions, res_hic$col.regions,my.ranges)
  
  # Finally, create the InteractionSet object
  iset <- InteractionSet(as.matrix(dplyr::select(res_hic,c(-1,-2))), gi)
  total_reads <- colSums(assay(iset))
  # If filters and filters2 exist, get bins which correspond to theses coordinates, and merged then
  if(!is.null(filters)){
    if(!is.null(filters2)){
      out <- linkOverlaps(iset, filters,filters2) # extract bin corresponding to filters and filters2
    }else{
      out <- linkOverlaps(iset, filters) # extract bin corresponding to filters 
      filters2 <- filters
    }
    if(nrow(out)==0)
      return(NULL)
    iset <- iset[out$query] # subtrat matrix
    out <- linkOverlaps(iset, filters,filters2) # re-do the overlap to get the new coordinates
    res_hic <- cbind(out[,-1],assay(iset[out$query])) %>% as_tibble() %>% group_by(subject1,subject2) %>% summarise_all(sum) %>% ungroup()
    gi <- GInteractions(filters[res_hic$subject1],filters2[res_hic$subject2])
    iset <- InteractionSet(as.matrix(dplyr::select(res_hic,c(-1,-2))), gi)
  }
  iset$totals <- total_reads
  interactions(iset) <- as(interactions(iset), "ReverseStrictGInteractions")
  colnames(iset) <- names(names_cond)
  metadata(iset)$width <- median(width(regions(iset)))
  names(assays(iset)) <- "counts"
  
  return(iset)
}

Get1val <- function(my.wigs,one.w,x){
  lapply(split(x,droplevels(seqnames(x))),function(zz){
    message(unique(as.character(seqnames(zz))))
    cov <- one.w[[unique(as.character(seqnames(zz)))]]
    score <- Views( cov, start = start(zz), end = end(zz) ) %>% sum()
    tibble(wig = my.wigs,value = score,rowname = zz$name)
  }) %>% bind_rows()
}

binSize <- 100000
n_percentile <- 21
OEtype <- "observed"
Chr.V <- glue::glue("chr{c(1,17,'X')}") %>% as.character()
# Chr.V <- paste0("chr",c(1:22,"X"))
# full_combi <- expand_grid(chr1=Chr.V,chr2=Chr.V)
full_combi <- combn(Chr.V,2) %>% t() %>% as.matrix() %>% as_tibble()
# full_combi <- rbind(tibble(V1=Chr.V,V2=Chr.V),full_combi)
full_combi <- tibble(V1=Chr.V,V2=Chr.V)
SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
# Build the bins on chromosomes of interest

my.ranges <- tileGenome(SeqInfo[Chr.V],tilewidth=binSize, cut.last.tile.in.chrom=TRUE) %>%
  sortSeqlevels() %>% sort()

mes_DSB <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/DF_allAsiSI_ordBLESSpOHT_fragPE_Rmdups_500bp_24012019.tsv" %>% read_tsv() %>% 
  arrange(desc(value)) %>% dplyr::slice(1:200) %>% as_granges()
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 2000000)

my.ranges <- my.ranges %>% filter_by_non_overlaps(mes_DSB)


indices_ranges <- my.ranges%>% mutate(idx = 1:length(my.ranges)) %>% as_tibble() %>% dplyr::select(seqnames,start,idx)

reskb <- glue::glue("{binSize/10e2}kb")
# AB_comp <- glue::glue("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC/compartmentAB/compartmentAB_DIvA_manipA_{reskb}.bedGraph") %>% read_bed_graph()
Dcomp_GR <- "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/PC1_all_chr_log2ratio_100kb_OHT_manipA.bw" %>% import.bw()

# AB_values <- AB_values %>% mutate(AB = ifelse(score<0,"B","A")) %>% split(.,.$AB) %>% 
#   lapply(function(x){
#     x %>% reduce()
#   }) %>% as("GRangesList") %>% unlist() %>% sortSeqlevels() %>% sort()

for(OEtype in c("observed","OE")){
  
  
  
  filetypename <- glue::glue("dump_{OEtype}_KR_")
  
  
  
  
  
  # HiC.files.manipD <- list.files("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/TRANS/KR",pattern=glue::glue("dump_{OEtype}_.+_{reskb}.txt.gz"),full.names=T)
  HiC.files.manipAline <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Aline/dump/TRANS/KR",pattern=glue::glue("dump_{OEtype}_.+_{reskb}.txt.gz"),full.names=T)
  HiC.files.manipAB <- list.files(glue::glue("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/matrix_generation/matrix_generation/data_retry/KR"),pattern=glue::glue("dump_{OEtype}_.+_{reskb}.txt.gz"),full.names=T)
  
  HiC.files <- c(HiC.files.manipAline,HiC.files.manipAB)
  
  conditionList <-c("DIvA","OHT")
  ctrlCOND <- conditionList[1]
  
  
  # DSB174 <- read_bed("/mnt/NAS/DATA/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")
  # 
  # DSB_n_chr <- DSB174 %>% as_tibble() %>% dplyr::count(seqnames)
  # most_DSB <- DSB_n_chr %>% pull(n) %>% which.max
  # less_DSB <- DSB_n_chr %>% pull(n) %>% which.min
  # mes_chrs <- DSB_n_chr %>% slice(most_DSB,less_DSB) %>% pull(seqnames) %>% as.character()
  # 
  
  full_res <- apply(full_combi,1,function(mes_chr){
    list_files <- HiC.files %>% str_subset(glue::glue("{mes_chr[1]}_{mes_chr[2]}_"))
    process_HiC(list_files,indices_ranges=indices_ranges,binSize=binSize)
  })
  
  
  full_res <- full_res %>%  plyr::compact()
  
  total_counts <- full_res %>% map(function(x){x$totals}) %>% Reduce("+",.)
  
  full_res <- full_res %>% map(function(x){
    x$totals <- total_counts
    return(x)
  })
  
  full_res <- do.call("rbind",full_res)
  colnames(full_res) <- str_replace(colnames(full_res),"DIVA","DIvA")
  
  for(n_percentile in c(5,11,21,31)){
    AB_values <- Dcomp_GR %>% as_tibble() %>% filter(score != 0) %>% group_by(seqnames) %>% nest() %>% mutate(data = map(data,function(x){
      x$quartile <- with(x, cut(score,
                                breaks=quantile(score, probs=seq(0,1, by=1/n_percentile), na.rm=TRUE),
                                include.lowest=TRUE))
      x %>% mutate(quartile = as.numeric(quartile))
    })) %>% unnest(data) %>% as_granges()
    
    
    DSB_pos_1 <- anchors(full_res)[[1]] %>% plyranges::join_overlap_left(AB_values)
    DSB_pos_2 <- anchors(full_res)[[2]] %>% plyranges::join_overlap_left(AB_values)
    
    data.OE.AB <- tibble(seqnames1 = as.vector(seqnames(DSB_pos_1)),seqnames2 = as.vector(seqnames(DSB_pos_2)),POS1=start(DSB_pos_1),POS2=start(DSB_pos_2),percentile_n_1=DSB_pos_1$quartile,percentile_n_2=DSB_pos_2$quartile) %>% 
      cbind(assay(full_res)) %>% 
      gather(key=Condition,value = score,-POS1,-POS2,-seqnames1,-seqnames2,-percentile_n_1,-percentile_n_2) %>% 
      dplyr::rename(Experiment = "Condition") %>% 
      mutate(Condition=str_extract(Experiment,paste(glue::glue("{conditionList}"),collapse="|"))) %>% 
      mutate(Experiment=str_remove(Experiment,paste(glue::glue("{conditionList}"),collapse="|"))) %>% 
      mutate(Condition = ifelse(Condition == ctrlCOND,"CTRL",Condition)) 
    
    
    data.OE.AB.cis <- data.OE.AB %>% filter(seqnames1==seqnames2)
    data.OE.AB.dup <- data.OE.AB.cis %>% filter(POS1!=POS2) %>%
      dplyr::rename(POS1="POS2",POS2="POS1")
    data.OE.AB.cis <- data.OE.AB.cis %>% rbind(data.OE.AB.dup)
    rm(data.OE.AB.dup)
    
    
    # p1 <- data.OE.AB.cis %>%
    #   filter(Condition == "CTRL") %>%
    #   ggplot(aes(x=POS1,y=POS2,fill=log2(score))) +
    #   geom_tile() +
    #   facet_wrap(Condition~seqnames1,scales="free") +
    #   scale_fill_gradient2(low="blue",midpoint=0,mid = "white",high="red",expand=F) +
    #   coord_cartesian(expand=0) +
    #   theme_classic() +
    #   theme(axis.title.x=element_blank(),
    #         axis.title.y=element_blank())
    
    
    
    ##PERCENTILES OF PC1
    data.OE.AB.summarised.ratio <- data.OE.AB.cis %>% spread(key=Condition,value=score) %>% gather(key=Condition,value=score,-seqnames1:-Experiment,-CTRL)
    data.OE.AB.summarised.ratio <- data.OE.AB.summarised.ratio %>% mutate(ratio = log2((score+1)/(CTRL+1)))
    data.OE.AB.summarised.ratio <- data.OE.AB.summarised.ratio %>% drop_na() %>%
      group_by(percentile_n_1,percentile_n_2,Condition,Experiment) %>%
      summarise(score =mean(ratio,na.rm=T)) %>%
      mutate(Experiment=str_remove(Experiment,"^_|_$")) %>% 
      mutate(Condition = glue::glue("{Condition}/CTRL"))
    data.OE.AB.summarised <- data.OE.AB.cis %>% drop_na() %>%  group_by(percentile_n_1,percentile_n_2,Condition,Experiment) %>% summarise(score =log2( mean(score)))
    
    limite <- data.OE.AB.summarised.ratio %>% pull(score) %>% quantile(0.999) %>% as.numeric()
    p1.ratio <- data.OE.AB.summarised.ratio  %>% 
      mutate(score = ifelse(score > limite,limite,score)) %>%
      mutate(score = ifelse(score < -limite,-limite,score)) %>%
      ggplot(aes(x=percentile_n_1,y=percentile_n_2,fill=score)) + 
      geom_tile() +
      facet_grid(Experiment~Condition) +
      scale_fill_gradient2(low = "#2980b9",high = "#f1c40f",midpoint = 0,mid = "black",expand=F) +
      coord_cartesian(expand=0) +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="bottom") 
    p2.ratio <- map(data.OE.AB.summarised.ratio$Condition %>% unique,function(x){
      AB_values %>% as_tibble() %>% mutate(Condition = x)
    }) %>% bind_rows() %>% ggplot(aes(x=quartile,y=score)) + geom_bar(stat="summary",fun="mean") +facet_grid(~Condition) + theme_classic() + coord_cartesian(expand=F)
    
    p1.ratio <- p2.ratio+ p1.ratio + plot_layout(height = c(0.1,0.9))
    
    
    
    # trymoica <- data.OE.AB.summarised %>% filter(Condition == "CTRL") %>% dplyr::select(-Condition)
    # trymoica <- sparseMatrix(i=trymoica[[1]],j=trymoica[[2]],x=trymoica[[3]])
    # 
    # limite <- 1
    # data.OE.AB.summarised <- data.OE.AB.summarised %>%
    #   mutate(score = ifelse(score > limite,limite,score)) %>%
    #   mutate(score = ifelse(score < -limite,-limite,score))
    
    p1.summarised <- data.OE.AB.summarised %>%
      mutate(Experiment=str_remove(Experiment,"^_|_$")) %>% 
      ggplot(aes(x=percentile_n_1,y=percentile_n_2,fill=score)) + 
      geom_tile() +
      facet_grid(Experiment~Condition) +
      scale_fill_gradient2(low="blue",midpoint=0,mid = "white",high="red",expand=F) +
      coord_cartesian(expand=0) +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="bottom") 
    
    p2.summarised <- map(conditionList,function(x){
      AB_values %>% as_tibble() %>% mutate(Condition = x)
    }) %>% bind_rows() %>% ggplot(aes(x=quartile,y=score)) + geom_bar(stat="summary",fun="mean") +facet_grid(~Condition) + theme_classic() + coord_cartesian(expand=F)

    p1.summarised <- p2.summarised+ p1.summarised + plot_layout(height = c(0.1,0.9))
    
    cc <- data.OE.AB.summarised %>% spread(key=Condition,value=score) %>% gather(key = Condition,value = score,-percentile_n_1:-percentile_n_2,-CTRL,-Experiment) %>% 
      mutate(diff = (score-CTRL)) %>%
      mutate(Condition = glue::glue("{Condition}/CTRL"))
    
    
    limite <- cc %>% pull(diff) %>% quantile(0.999) %>% as.numeric()
    p3.summarised <- cc %>%
      mutate(diff = ifelse(diff > limite,limite,diff)) %>%
      mutate(diff = ifelse(diff < -limite,-limite,diff)) %>%
      ggplot(aes(x=percentile_n_1,y=percentile_n_2,fill=diff)) + 
      geom_tile() +
      facet_grid(Experiment~Condition) +
      scale_fill_gradient2(low = "#2980b9",high = "#f1c40f",midpoint = 0,mid = "black",expand=F) +
      coord_cartesian(expand=0) +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="bottom") 
    
    # p2.summarised <- AB_values %>% as_tibble() %>% ggplot(aes(x=quartile,y=score)) + geom_bar(stat="summary",fun="mean") + theme_classic() + coord_cartesian(expand=F)
    
    p2.summarised <- map(cc$Condition %>% unique,function(x){
      AB_values %>% as_tibble() %>% mutate(Condition = x)
    }) %>% bind_rows() %>% ggplot(aes(x=quartile,y=score)) + geom_bar(stat="summary",fun="mean") +facet_grid(~Condition) + theme_classic() + coord_cartesian(expand=F)
    
    
    p3.summarised <- p2.summarised+ p3.summarised + plot_layout(height = c(0.1,0.9))
    
    # #ratio
    # data.ratio.AB <-  tibble(seqnames1 = as.vector(seqnames(DSB_pos_1)),seqnames2 = as.vector(seqnames(DSB_pos_2)),POS1=start(DSB_pos_1),POS2=start(DSB_pos_2),percentile_n_1=DSB_pos_1$quartile,percentile_n_2=DSB_pos_2$quartile) %>%
    #   cbind((assay(full_res)*((1000000/full_res$totals)))) %>%
    #   gather(key=Condition,value = score,-seqnames1:-percentile_n_2) %>%
    #   mutate(Condition=str_extract(Condition,paste(glue::glue("{conditionList}$"),collapse="|"))) %>%
    #   mutate(Condition = ifelse(Condition == ctrlCOND,"CTRL",Condition)) %>%
    #   mutate(score = score+1) %>%
    #   # group_by(pos1,pos2,Condition) %>% summarise(score = mean(score,na.omit=T)) %>%
    #   spread(key = Condition,value = score) %>%
    #   gather(key = Condition,value = score,-seqnames1:-percentile_n_2,-CTRL) %>%
    #   mutate(ratio = log2(score/CTRL)) %>%
    #   mutate(Condition = glue::glue("{Condition}/CTRL"))
    # 
    # 
    # data.OE.AB.cis.ratio <- data.ratio.AB %>% filter(seqnames1==seqnames2)
    # data.OE.AB.dup <- data.OE.AB.cis.ratio %>% filter(POS1!=POS2) %>%
    #   dplyr::rename(POS1="POS2",POS2="POS1")
    # data.OE.AB.cis.ratio <- data.OE.AB.cis.ratio %>% rbind(data.OE.AB.dup)
    # rm(data.OE.AB.dup)
    
    # p1 <- data.ratio.AB %>%
    #   filter(Condition == "OHT/CTRL") %>%
    #   ggplot(aes(x=POS1,y=POS2,fill=log2(score))) +
    #   geom_tile() +
    #   facet_wrap(Condition~seqnames1,scales="free") +
    #   scale_fill_gradient2(low="blue",midpoint=0,mid = "white",high="red",expand=F) +
    #   coord_cartesian(expand=0) +
    #   theme_classic() +
    #   theme(axis.title.x=element_blank(),
    #         axis.title.y=element_blank())
    
    # data.OE.AB.summarised.ratio <- data.OE.AB.cis.ratio %>% drop_na() %>%  group_by(percentile_n_1,percentile_n_2,Condition) %>% summarise(ratio =log2( mean(ratio)))
    
    
    # trymoica <- data.OE.AB.summarised %>% filter(Condition == "CTRL") %>% dplyr::select(-Condition)
    # trymoica <- sparseMatrix(i=trymoica[[1]],j=trymoica[[2]],x=trymoica[[3]])
    #
    # limite <- 0.75
    # data.OE.AB.summarised.ratio <- data.OE.AB.summarised %>%
    #   mutate(score = ifelse(score > limite,limite,score)) %>%
    #   mutate(score = ifelse(score < -limite,-limite,score))
    
    # p3.summarised <- data.OE.AB.summarised.ratio %>%
    #   ggplot(aes(x=percentile_n_1,y=percentile_n_2,fill=score)) +
    #   geom_tile() +
    #   facet_grid(~Condition) +
    #   scale_fill_gradient2(low="blue",midpoint=0,mid = "white",high="red",expand=F) +
    #   coord_cartesian(expand=0) +
    #   theme_classic() +
    #   theme(axis.title.x=element_blank(),
    #         axis.title.y=element_blank())
    # 
    # p2.summarised <- AB_values %>% as_tibble() %>% ggplot(aes(x=quartile,y=score)) + geom_bar(stat="summary",fun="mean") + theme_classic() + coord_flip(expand=F)
    # 
    # p2.summarised+ p1.summarised + plot_layout(widths = c(0.1,0.9))
    
    
    pdf(glue::glue("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/results/heatmap_HiC_D/SADDLE_PLOTS/Heatmap_compD_saddle_{reskb}_percentile_{n_percentile}_{OEtype}.pdf"),width=5,height=23)
    # print(p1)
    print(p1.ratio)
    print(p1.summarised)
    print(p3.summarised)
    # print(p2.summarised)
    dev.off()
    
  }
  
  
}
