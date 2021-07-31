
#from reproduce_françois_plot
require(HiTC)
require(tidyverse)
require(plyranges)
require(reshape2)
require(keras)
require(rtracklayer)
require(cowplot)
require(edgeR)
require(ggside)
seqlens <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19 %>% seqlengths()

#Comp analysis
#Make compA/B based on H3K9me3
Get1val <- function(one.w,my.wigs,x){
  lapply(split(x,droplevels(seqnames(x))),function(zz){
    message(unique(as.character(seqnames(zz))))
    cov <- one.w[[unique(as.character(seqnames(zz)))]]
    score <- Views( cov, start = start(zz), end = end(zz) ) %>% sum()
    tibble(wig = my.wigs,value = score,rowname = zz$name)
  }) %>% bind_rows()
}
compAB <- "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC/compartmentAB/compartmentAB_DIvA_manipA_100kb.bedGraph" %>% 
  read_bed_graph()
# compAB$name <- paste(compAB)
# 
# 
# mes_bw <- list("H3K9me3_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K9me3/EMBL_H2VFWBGX5_HWNTLBGX3/PROCESSED/mapping/BIGWIG/HWNTLBGX3_H3K9me3_DIVA_17s005723-1-1_Clouaire_lane117s005723_sequence_normalized.bw","h3k79me2_mOHT"="/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/HISTONES/H3K79me2/BGI_RUN1_201203/PROCESSED/BIGWIG/h3k79me2_mOHT_normalized_hg19.bw")
# mes_bw <- mes_bw %>% map(import.bw,as="RleList")
# 
# 
# res_bw <- names(mes_bw) %>% 
#   map(function(x){
#     Get1val(mes_bw[[x]],x,compAB)
#   }) %>% bind_rows() %>% spread(key = wig,value =value)
# 
# compAB <- compAB %>% as_tibble() %>% left_join(res_bw,by = c("name"="rowname")) %>% dplyr::select(-name)
# 
# compAB_corr <- compAB %>% dplyr::select(seqnames,score,h3k79me2_mOHT,H3K9me3_mOHT) %>% 
#   group_by(seqnames) %>% nest() %>% 
#   mutate(corr_data = map(data,correlate)) %>%
#   mutate(corr_data = map(corr_data,focus,score)) %>%
#   dplyr::select(-data) %>% unnest(corr_data) 
# 
# compAB_corr %>% spread(key = rowname,value = score) %>% write_tsv("../../results/correlation_PC1_diva_with_histone_marks.tsv")





#Functions

get_all_against_all_matrix_Hic <- function(file,bed,add_for_log = 0){
  c(my.ranges,intdata) %<-% readRDS(file)
  
  val1 <- intdata %>% sum
  overlaps <- as(findOverlaps(bed,my.ranges), "List")
  
  
  res <- mclapply(1:length(overlaps),function(i){
    x <- overlaps[[i]]
    sres <- lapply(1:length(overlaps),function(j){
      y <- overlaps[[j]]
      intdata[x,y] %>% sum() %>% round()
    })
    names(sres) <- bed$name
    
    sres <- sres %>% melt() %>% as_tibble() %>% mutate(bin1 = bed[i]$name)
    colnames(sres) <- c("score","bin2","bin1")
    sres
  },mc.cores=12) %>% bind_rows()
  
  # overlaps <- findOverlaps(bed,my.ranges)%>% as_tibble
  # dfbed <- bed %>% as_tibble() %>%
  #     mutate(index = 1:length(bed)) %>%
  #     left_join(overlaps,by = c("index"="queryHits")) %>% dplyr::select(name,subjectHits) %>% mutate(binName = str_c(name,subjectHits,sep="_"))
  # submat <- intdata[dfbed$subjectHits,dfbed$subjectHits] %>% as.matrix() %>% as.data.frame()
  # colnames(submat) <- rownames(submat) <- dfbed$binName
  #
  # res <- submat %>% tibble::rownames_to_column("bin1") %>% gather(key = bin2,value = score,-bin1) %>%
  #     mutate_if(is_character,str_remove,"_[0-9]+") %>% group_by(bin1,bin2) %>% summarise(score = sum(score,na.rm = T)) %>% ungroup() %>%
  #     mutate(score = score + add_for_log)
  
  
  
  
  return(list(res,val1))
}

process_bin_for_comp <- function(file,bed){
  message(file)
  c(my.ranges,intdata) %<-% readRDS(file)
  # ov <- findOverlaps(my.ranges,bed[bed$name %in% c("bin10004","bin10001")])
  ov <- findOverlaps(my.ranges,bed)
  val1 <- intdata %>% sum
  intdata <- intdata[queryHits(ov),queryHits(ov)]
  dimnames(intdata) = list(bed$name,bed$name)
  df <- as.data.frame(summary(intdata))
  df$bin1 <- rownames(intdata)[df$i]
  df$bin2 <- colnames(intdata)[df$j]
  df <- df %>% dplyr::select(-i,-j) %>% 
    # filter(seqnames.x!=seqnames.y) %>% 
    mutate(x = x * (1000000/val1)) 
  return(df)
}

DSB174 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")
bless80 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
HR <-  read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/BLESS_HR_JunFragPE_Rmdups_pm500bp.bed")
NHEJ <-  read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/BLESS_NHEJ_JunFragPE_Rmdups_pm500bp.bed")
uncut <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/ASIsites_hg19_without_bless80.bed")

compAB_DIVA <- compAB %>% as_granges() %>% filter_by_non_overlaps(DSB174 %>% anchor_center() %>% mutate(width = 2000000)) %>% 
  filter(score != 0) %>% 
  mutate(comp = ifelse(score > 0 ,"A","B"))

compAB_DIVA$name <- glue::glue("bin{1:length(compAB_DIVA)}")

compAB_DIVA <- split(compAB_DIVA,compAB_DIVA$comp) %>% lapply(sample,5000) %>% as("GRangesList") %>% unlist()

compAB_DIVA_infos <- compAB_DIVA %>% as_tibble() %>% dplyr::select(seqnames,comp,name)

m.w <- "100kb"
rds.files <- list.files("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/TRANS/KR/HiTC/observed",pattern=str_c("_",m.w),full.names=T)

Replicate <- rds.files %>% basename() %>% str_remove_all("HTC_observed_KR_HiC_D_|_full_matrix_[0-9]+kb.rds|HTC_observed_KR_|_manipA")
names(rds.files) <- Replicate

# res_DIVA <- process_bin_for_comp(file = rds.files[["DIvA"]],bed=compAB_DIVA)

res <- lapply(rds.files[Replicate[1:2]],function(file){
  
  process_bin_for_comp(file = file,bed=compAB_DIVA)
}) %>% setNames(Replicate[1:2]) %>% bind_rows(.id = "Type")


# res[["CTRL"]] <- res_DIVA
# res <- res %>%  reduce(left_join, by = c("bin1"="bin2"))

res <- res %>% 
  spread(key = Type,value = x) %>%
  drop_na() %>% 
  gather(key = Condition,value = score,-bin2,-bin1,-DIvA) %>% 
  mutate(ratio = log2(score/DIvA)) %>%
  mutate(Condition = glue::glue("{Condition}/DIvA"))

res <- res %>% 
  left_join(compAB_DIVA_infos,by = c("bin1"="name")) %>% 
  left_join(compAB_DIVA_infos,by = c("bin2"="name"))

res <- res %>% filter(seqnames.x!=seqnames.y)

res <- res %>% unite(compAB,comp.x,comp.y)

res %>% group_by(compAB,Condition) %>% summarise(mean(ratio,na.rm=T))

p1 <- res %>% ggplot(aes(x=compAB,y=ratio,fill=Condition)) + geom_boxplot() + facet_zoom(ylim = c(-1,1))

# cc <- read_tsv("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/DSB_DSB_interaction_1mbEach.tsv.gz")
for(ctrlCOND in c("DIvA","OHT")){
  for(m.w in names(mes_windows)){
    message(m.w)
    sw <- mes_windows[[m.w]]
    rds.files <- list.files("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/TRANS/KR/HiTC/observed",pattern=str_c("_",m.w),full.names=T)
    # rds.files <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/matrix_generation/data_retry/KR/HiTC/observed",pattern=str_c("_",m.w),full.names=T) %>% str_subset("manipA")
    Replicate <- rds.files %>% basename() %>% str_remove_all("HTC_observed_KR_HiC_D_|_full_matrix_[0-9]+kb.rds|HTC_observed_KR_|_manipA")
    names(rds.files) <- Replicate
    
    for(window in mes_limites[[m.w]]){
      message(window)
      
      
      diva.chr <- bless80 %>%
        anchor_center() %>% mutate(width = window)
      if(window/1e6 >= 1){
        window_format <- scales::number(window/2,scale=1e-6,suffix="Mb")
      }else{
        window_format <- scales::number(window/2,scale=1e-3,suffix="Kb")
      }
      
      
      
      res <- lapply(rds.files,function(i){
        get_all_against_all_matrix_Hic(file = i,bed = diva.chr,add_for_log = 1)
      })
      
      f_ctrl <- Replicate %>% str_subset(glue::glue("{ctrlCOND}$"))
      # f_exp <- Replicate %>% str_subset(glue::glue("{ctrlCOND}$"),negate=T)
      
      res_norm <- lapply(Replicate,function(i){
        # res[[i]][[1]] %>% mutate(val1 = res[[f_ctrl]][[2]]) %>%  mutate(val2 = res[[i]][[2]])  %>% mutate(Condition = i) %>%
        #   mutate(score = score * (val1/val2)) %>% dplyr::select(-val1,-val2)
        res[[i]][[1]] %>% mutate(val2 = res[[i]][[2]])  %>% mutate(Condition = i) %>%
          mutate(score = score * (1000000/val2)) %>% dplyr::select(-val2)%>% mutate(Condition = ifelse(Condition == ctrlCOND,"CTRL",Condition))
      })%>% setNames(Replicate)
      
      
      # 
      # c(res1,val1) %<-% get_all_against_all_matrix_Hic(file = rds.files[[1]],bed = diva.chr,add_for_log = 1)
      # c(res2,val2) %<-% get_all_against_all_matrix_Hic(file = rds.files[[2]],bed = diva.chr,add_for_log = 1)
      # c(res3,val3) %<-% get_all_against_all_matrix_Hic(file = rds.files[[3]],bed = diva.chr,add_for_log = 1)
      # c(res4,val4) %<-% get_all_against_all_matrix_Hic(file = rds.files[[4]],bed = diva.chr,add_for_log = 1)
      # c(res5,val5) %<-% get_all_against_all_matrix_Hic(file = rds.files[[5]],bed = diva.chr,add_for_log = 1)
      
      # res1 <- res[[1]][[1]] %>% mutate(val1 = res[[1]][[2]]) %>%  mutate(val2 = res[[2]][[2]]) %>% mutate(Condition = Replicate[[1]])
      # res2 <- res[[2]][[1]] %>% mutate(val1 = res[[1]][[2]]) %>%  mutate(val2 = res[[2]][[2]])  %>% mutate(Condition = Replicate[[2]]) %>%
      #   mutate(score = score * (val1/val2))
      # res3 <- res3 %>% mutate(val1 = val1) %>%  mutate(val3 = val3)  %>% mutate(Condition = Replicate[[2]]) %>%
      #   mutate(score = score * (val1/val3))
      # res4 <- res4 %>% mutate(val1 = val1) %>%  mutate(val4 = val4)  %>% mutate(Condition = Replicate[[2]]) %>%
      #   mutate(score = score * (val1/val4))
      # res5 <- res5 %>% mutate(val1 = val1) %>%  mutate(val5 = val5)  %>% mutate(Condition = Replicate[[2]]) %>%
      #   mutate(score = score * (val1/val5))
      
      # data.ratio.AB <- rbind(
      #   res[[f_ctrl]][[1]] %>% mutate(Condition = "CTRL"),
      #   res_norm %>% bind_rows()
      #   ) %>%
      data.ratio.AB <- res_norm %>% bind_rows() %>% 
        mutate(score = score+1) %>% 
        spread(key = Condition,value = score) %>%
        gather(key = Condition,value = score,-bin2,-bin1,-CTRL) %>% 
        mutate(ratio = log2(score/CTRL)) %>%
        mutate(Condition = glue::glue("{Condition}/{ctrlCOND}"))
      bless40 <- bless80 %>% anchor_center() %>% mutate(width=80000)
      gamma.val <- Get1val(x = bless40,my.wigs = names(mes.bw1)[1],one.w = mes.bw1[[1]]) %>% arrange(desc(value)) %>% dplyr::select(value,rowname)
      p53BP1.val <- Get1val(x = bless40,my.wigs = names(mes.bw2)[1],one.w = mes.bw2[[1]]) %>% arrange(desc(value)) %>% dplyr::select(value,rowname)
      bless40 <- bless80 %>% anchor_center() %>% mutate(width=2000)
      pol2.val <- Get1val(x = bless40,my.wigs = names(mes.bw3)[1],one.w = mes.bw3[[1]]) %>% arrange(desc(value)) %>% dplyr::select(value,rowname)
      #for chromosomic order
      
      order_chr_bless <- bless80 %>% sortSeqlevels %>% sort() %>% mutate(pos = 1:n()) %>%  as_tibble() %>% dplyr::select(name,pos,seqnames)
      chromo_color <- tibble(seqnames = unique(order_chr_bless$seqnames)) %>% 
        mutate(col_seq = pals::kelly(dplyr::n()))
      order_chr_bless <-  order_chr_bless%>% 
        left_join(chromo_color,by = "seqnames")
      
      
      
      data.ratio.AB.trans <- data.ratio.AB %>% 
        left_join(order_chr_bless,by = c("bin1"="name")) %>% dplyr::rename(pos_bin1 = "pos")%>% dplyr::rename(chr_bin1 = "seqnames") %>% dplyr::rename(col_seq_bin1 = "col_seq") %>% 
        left_join(order_chr_bless,by = c("bin2"="name")) %>% dplyr::rename(pos_bin2 = "pos")%>% dplyr::rename(chr_bin2 = "seqnames") %>% dplyr::rename(col_seq_bin2 = "col_seq") %>% 
        mutate(bin1 = fct_reorder(bin1,pos_bin1,.desc = T))  %>% 
        mutate(bin2 = fct_reorder(bin2,pos_bin2,.desc = T))
      
      
      p0 <- plot_diff(data.ratio.AB.trans,facet=T,fixed=T,fixed.limite=lim.Plot[[m.w]],my.quantile=0.99) +
        geom_ysidetile(aes(x = "chromosomes", yfill = `chr_bin2`)) +
        geom_xsidetile(aes(y = "chromosomes", xfill = `chr_bin1`),show.legend =F) +
        scale_xfill_manual(values = chromo_color$col_seq,breaks = chromo_color$seqnames) +
        scale_yfill_manual("chromosomes",values = chromo_color$col_seq,breaks = chromo_color$seqnames) +
        guides(fill = guide_colourbar(title.position = "top",title = "log2(Fold Change)",title.hjust=0.5,barwidth =unit(20,"lines"),barheight =unit(1,"lines"))) +
        ggtitle("ORDERED BY CHROM")
      #for gamma
      data.ratio.AB.gamma <- data.ratio.AB %>%
        left_join(gamma.val,by = c("bin1"="rowname")) %>% dplyr::rename(value_gamma_1 = "value") %>% 
        left_join(gamma.val,by = c("bin2"="rowname")) %>% dplyr::rename(value_gamma_2 = "value") %>% 
        mutate(bin1 = fct_reorder(bin1,value_gamma_1)) %>%
        mutate(bin2 = fct_reorder(bin2,value_gamma_2))
      
      p1 <- plot_diff(data.ratio.AB.gamma,facet=T,fixed=T,fixed.limite=lim.Plot[[m.w]],my.quantile=0.99) +
        geom_ysidetile(aes(x = "GAMMA", yfill = log10(value_gamma_2)),show.legend = F)  +
        geom_xsidetile(aes(y = "GAMMA", xfill = log10(value_gamma_1)),show.legend = F)  +
        scale_xfill_gradient(low = "white",high = "red") +
        scale_yfill_gradient(low = "white",high = "red") +
        guides(fill = guide_colourbar(title.position = "top",title = "log2(Fold Change)",title.hjust=0.5,barwidth =unit(20,"lines"),barheight =unit(1,"lines"))) +
        ggtitle("ORDERED BY GAMMA")
      
      #for 53BP1
      data.ratio.AB.53BP1 <- data.ratio.AB %>%
        left_join(p53BP1.val,by = c("bin1"="rowname")) %>% dplyr::rename(value_53BP1_1 = "value") %>% 
        left_join(p53BP1.val,by = c("bin2"="rowname")) %>% dplyr::rename(value_53BP1_2 = "value") %>% 
        mutate(bin1 = fct_reorder(bin1,value_53BP1_1)) %>%
        mutate(bin2 = fct_reorder(bin2,value_53BP1_2))
      
      p12 <- plot_diff(data.ratio.AB.53BP1,facet=T,fixed=T,fixed.limite=lim.Plot[[m.w]],my.quantile=0.99) +
        geom_ysidetile(aes(x = "53BP1", yfill = log10(value_53BP1_2)),show.legend = F)  +
        geom_xsidetile(aes(y = "53BP1", xfill = log10(value_53BP1_1)),show.legend = F)  +
        scale_xfill_gradient(low = "white",high = "red") +
        scale_yfill_gradient(low = "white",high = "red") +
        guides(fill = guide_colourbar(title.position = "top",title = "log2(Fold Change)",title.hjust=0.5,barwidth =unit(20,"lines"),barheight =unit(1,"lines"))) +
        ggtitle("ORDERED BY 53BP1")
      
      
      #Order by HR NHEJ
      orderHR <- bless80 %>% as_tibble() %>% mutate(Type = case_when(
        name %in% HR$name ~ 1,
        name %in% NHEJ$name ~ 2,
        TRUE ~ 0
      )) %>% arrange(Type) %>% filter(Type != 0) %>% 
        mutate(name = glue::glue("<i style='color:{ifelse(name %in% HR$name,'#c0392b','#2980b9')}'>{name}</i>"))
      data.ratio.AB.HR <- data.ratio.AB  %>% 
        mutate(bin1 = glue::glue("<i style='color:{ifelse(bin1 %in% HR$name,'#c0392b','#2980b9')}'>{bin1}</i>")) %>% 
        mutate(bin2 = glue::glue("<i style='color:{ifelse(bin2 %in% HR$name,'#c0392b','#2980b9')}'>{bin2}</i>"))%>% filter(bin1 %in% orderHR$name,bin2 %in% orderHR$name) %>% mutate(bin1 = fct_relevel(bin1,orderHR$name)) %>%
        mutate(bin2 = fct_relevel(bin2,orderHR$name))
      
      p2 <- plot_diff(data.ratio.AB.HR,facet=T,fixed=T,fixed.limite=lim.Plot[[m.w]]) + theme(axis.text.x= element_markdown(),axis.text.y= element_markdown()) +
        guides(fill = guide_colourbar(title.position = "top",title = "log2(Fold Change)",title.hjust=0.5,barwidth =unit(20,"lines"),barheight =unit(1,"lines")))+
        ggtitle("ORDERED BY HR/NHEJ")
      
      
      
      pdf(glue::glue("redo_plot_francois_ratio_pm{window_format}_bin{m.w}_Hich_average_ratio_{ctrlCOND}.pdf"),width=24,height=24)
      print(p0)
      print(p1)
      print(p12)
      print(p2)
      dev.off()
      
      #ADD UNCUT
      set.seed(174)
      diva.chr <- uncut %>%
        anchor_center() %>% mutate(width = window) %>% sample(80)
      
      
      res.uncut <- lapply(rds.files,function(i){
        get_all_against_all_matrix_Hic(file = i,bed = diva.chr,add_for_log = 1)
      })
      
      f_ctrl <- Replicate %>% str_subset(glue::glue("{ctrlCOND}$"))
      # f_exp <- Replicate %>% str_subset(glue::glue("{ctrlCOND}$"),negate=T)
      
      res_norm.uncut <- lapply(Replicate,function(i){
        # res.uncut[[i]][[1]] %>% mutate(val1 = res[[f_ctrl]][[2]]) %>%  mutate(val2 = res.uncut[[i]][[2]])  %>% mutate(Condition = i) %>%
        #   mutate(score = score * (val1/val2)) %>% dplyr::select(-val1,-val2)
        res.uncut[[i]][[1]] %>% mutate(val2 = res.uncut[[i]][[2]])  %>% mutate(Condition = i) %>%
          mutate(score = score * (1000000/val2)) %>% dplyr::select(-val2)%>% mutate(Condition = ifelse(Condition == ctrlCOND,"CTRL",Condition))
      })%>% setNames(Replicate)
      
      
      
      
      # data.ratio.AB.uncut <- rbind(
      #   res.uncut[[f_ctrl]][[1]] %>% mutate(Condition = "CTRL"),
      #   res_norm.uncut %>% bind_rows()
      # ) %>%
      data.ratio.AB.uncut <- res_norm.uncut %>% bind_rows() %>%
        mutate(score = score+1) %>%
        spread(key = Condition,value = score) %>%
        gather(key = Condition,value = score,-bin2,-bin1,-CTRL) %>%
        mutate(ratio = log2(score/CTRL)) %>%
        mutate(Condition = glue::glue("{Condition}/{ctrlCOND}"))
      #ADD RANDOM80
      diva.chr <- read_bed("/mnt/NAS/DATA/AsiSI/80random.bed") %>%
        anchor_center() %>% mutate(width = window)
      
      
      res.random <- lapply(rds.files,function(i){
        get_all_against_all_matrix_Hic(file = i,bed = diva.chr,add_for_log = 1)
      })
      
      f_ctrl <- Replicate %>% str_subset(glue::glue("{ctrlCOND}$"))
      # f_exp <- Replicate %>% str_subset(glue::glue("{ctrlCOND}$"),negate=T)
      
      res_norm.random <- lapply(Replicate,function(i){
        # res.random[[i]][[1]] %>% mutate(val1 = res[[f_ctrl]][[2]]) %>%  mutate(val2 = res.random[[i]][[2]])  %>% mutate(Condition = i) %>%
        #   mutate(score = score * (val1/val2)) %>% dplyr::select(-val1,-val2)
        res.random[[i]][[1]] %>% mutate(val2 = res.random[[i]][[2]])  %>% mutate(Condition = i) %>%
          mutate(score = score * (1000000/val2)) %>% dplyr::select(-val2)%>% mutate(Condition = ifelse(Condition == ctrlCOND,"CTRL",Condition))
      })%>% setNames(Replicate)
      
      
      
      
      # data.ratio.AB.random <- rbind(
      #   res.random[[f_ctrl]][[1]] %>% mutate(Condition = "CTRL"),
      #   res_norm.random %>% bind_rows()
      # ) %>%
      data.ratio.AB.random <- res_norm.random %>% bind_rows() %>%
        mutate(score = score+1) %>%
        spread(key = Condition,value = score) %>%
        gather(key = Condition,value = score,-bin2,-bin1,-CTRL) %>%
        mutate(ratio = log2(score/CTRL)) %>%
        mutate(Condition = glue::glue("{Condition}/{ctrlCOND}"))
      
      datacutUNCUT <-  rbind(
        data.ratio.AB %>% mutate(Type = "80bless"),
        data.ratio.AB.uncut %>% mutate(Type = "uncut"),
        data.ratio.AB.random %>% mutate(Type = "random")
      )
      
      res.boxplot1 <- datacutUNCUT %>% filter(bin1 != bin2) %>% ggplot(aes(x=Type,y=ratio,fill=Condition)) + geom_boxplot(outlier.shape = NA) +geom_hline(yintercept = 0,linetype="dashed",col="red")
      
      pdf(glue::glue("boxplot_pm{window_format}_bin{m.w}_Ratio_D_normalized_rpkm_{ctrlCOND}.pdf"),width=10,height=6)
      print(res.boxplot1+coord_cartesian(ylim=c(-lim.Plot[[m.w]],lim.Plot[[m.w]])) + ggtitle("80Bless vs uncut (sampled 80)"))
      dev.off()
    }
    
    
    
    
  }
}