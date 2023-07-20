#from process_diffHic_from_juicer.Rmd & gene_compD_drip.R
require(BSgenome.Hsapiens.UCSC.hg19)
require(diffHic)
require(edgeR)
require(tidyverse)
require(plyranges)
require(data.table)
require(normOffsets)
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

###SET REGIONS OF INTEREST
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
gamma_region <- mes_DSB %>% anchor_center() %>% mutate(width = 2000000)
# chr.to.study <- glue::glue("chr{c(2,6,9,13,18,1,17,20,'X')}")
chr.to.study <- glue::glue("chr{c(1)}")
# chr.to.study <- glue::glue("chr{c(1,17,'X')}")
# chr.to.study <- glue::glue("chr{c(1,17)}")
genes_DIVA <- ens.genes %>% plyranges::filter(gene_id %in% DE_DIVA$rowname)%>%
  filter_by_non_overlaps(gamma_region)
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
  dat1
},mc.cores=length(mes_bw)) %>% bind_rows()

res <- as_tibble(mes_genes) %>% left_join(DE_DIVA_INFO,by = c("gene_id"="rowname"))%>% right_join(res.plot,by = c("name"="rowname"))

qDRIPgenes <- res %>% filter(seqnames %in% chr.to.study) %>%  group_by(wig) %>% arrange(desc(value)) %>% mutate(Group = ntile(value,100))

##Take only the top10/bottom10 of theses genes

qDRIPgenes <- qDRIPgenes%>% mutate(Group = case_when(
  Group %in% 1:10 ~ "Bottom10",
  Group %in% 91:100 ~ "Top10",
  TRUE ~ "Other"
)) %>% filter(Group != "Other")

first_gr <- qDRIPgenes %>% pull(gene_id) %>% unique()
first_gr <- mes_genes %>% filter(gene_id %in% first_gr) %>% sortSeqlevels() %>% sort()

qDRIP_Group <- qDRIPgenes %>% dplyr::select(gene_id,Group,wig)

##The second set of regions is composed of qDRIP top genes AND gamma_domain (2Mb size)
bless80 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed") %>% anchor_center() %>% mutate(width = 1000000)
second_gr <- c(bless80,first_gr) %>% sortSeqlevels() %>% sort()


binSize <- 50000
OEtype <- "observed"
reskb <- glue::glue("{binSize/10e2}kb")
filetypename <- glue::glue("dump_{OEtype}_KR_")

Chr.V <- paste0("chr",c(1:22,"X"))


SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
# Build the bins on chromosomes of interest

my.ranges <- tileGenome(SeqInfo[Chr.V],tilewidth=binSize, cut.last.tile.in.chrom=TRUE) %>%
  sortSeqlevels() %>% sort()

indices_ranges <- my.ranges%>% mutate(idx = 1:length(my.ranges)) %>% as_tibble() %>% dplyr::select(seqnames,start,idx)

combi_chromosomes <- combn(Chr.V,2) %>% t() %>% as_tibble() %>%
  rbind(cbind(V1=Chr.V,V2=Chr.V),.) %>% 
  unite("chr_combi",V1:V2) %>% pull(chr_combi)


# HiC.files.manipD <- list.files("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/TRANS/KR",pattern=glue::glue("dump_{OEtype}_.+_{reskb}.txt.gz"),full.names=T)
HiC.files.manipAline <- list.files("/home/rochevin/Documents/PROJET_INGE/HiC_Aline/dump/TRANS/KR",pattern=glue::glue("dump_{OEtype}_.+_{reskb}.txt.gz"),full.names=T)
HiC.files.manipAB <- list.files(glue::glue("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/matrix_generation/matrix_generation/data_retry/KR"),pattern=glue::glue("dump_{OEtype}_.+_{reskb}.txt.gz"),full.names=T)

HiC.files <- c(HiC.files.manipAline,HiC.files.manipAB)

combi_files <- combi_chromosomes %>% map(function(x){str_subset(HiC.files,pattern=glue::glue("{x}_{reskb}"))})



conditionList <-c("DIvA","OHT")
ctrlCOND <- conditionList[1]

res_hic <- lapply(combi_files,process_HiC,indices_ranges=indices_ranges,filters=first_gr,filters2=second_gr,binSize=binSize) 

res_hic <- res_hic %>%  plyr::compact()

total_counts <- res_hic %>% map(function(x){x$totals}) %>% Reduce("+",.)

res_hic <- res_hic %>% map(function(x){
  x$totals <- total_counts
  return(x)
})

res_hic <- do.call("rbind",res_hic)
colnames(res_hic) <- str_replace(colnames(res_hic),"DIVA","DIvA")


DSB_pos_1 <- anchors(res_hic)[[1]] %>% mutate(POS = paste(anchors(res_hic)[[1]])) %>% as_tibble() %>% 
  left_join(as_tibble(second_gr %>% mutate(POS = paste(second_gr)))) %>% as_granges()

DSB_pos_2 <- anchors(res_hic)[[2]] %>% mutate(POS = paste(anchors(res_hic)[[2]])) %>% as_tibble() %>% 
  left_join(as_tibble(first_gr %>% mutate(POS = paste(first_gr)))) %>% as_granges()

if(OEtype=="observed"){
  data.ratio.AB <- tibble(pos1=as.character(DSB_pos_1$name),pos2=as.character(DSB_pos_2$name)) %>% 
    cbind((assay(res_hic)*((1000000/res_hic$totals)))) %>% 
    group_by(pos1,pos2) %>% summarise_all(sum) %>% ungroup() %>% 
    gather(key=Condition,value = score,-pos1,-pos2) %>% 
    dplyr::rename(Experiment = "Condition") %>% 
    mutate(Condition=str_extract(Experiment,paste(glue::glue("{conditionList}"),collapse="|"))) %>% 
    mutate(Experiment=str_remove(Experiment,paste(glue::glue("{conditionList}"),collapse="|"))) %>% 
    mutate(Condition = ifelse(Condition == ctrlCOND,"CTRL",Condition))  %>% 
    mutate(score = score+1) %>%  
    spread(key = Condition,value = score) %>%
    gather(key = Condition,value = score,-pos2,-pos1,-CTRL,-Experiment) %>% 
    mutate(ratio = log2(score/CTRL)) %>%
    mutate(Condition = glue::glue("{Condition}/DIvA"))
  data.ratio.AB <- data.ratio.AB %>% filter(pos1!=pos2)
  
  data.ratio.AB <- lapply(split(qDRIP_Group,qDRIP_Group$wig),function(x){
    data.ratio.AB %>% left_join(x,by=c("pos2"="gene_id")) %>% drop_na()
  }) %>% bind_rows()
  
  p1 <- data.ratio.AB %>% ggplot(aes(x=Experiment,y=ratio,fill=Group)) + geom_boxplot() + theme_classic() + facet_grid(~wig)
  
  
  
}else{
  data.ratio.AB <- tibble(pos1=as.character(DSB_pos_1$name),pos2=as.character(DSB_pos_2$name)) %>% 
    cbind((assay(res_hic)*((1000000/res_hic$totals)))) %>% 
    group_by(pos1,pos2) %>% summarise_all(sum) %>% ungroup() %>% 
    gather(key=Condition,value = score,-pos1,-pos2) %>% 
    dplyr::rename(Experiment = "Condition") %>% 
    mutate(Condition=str_extract(Experiment,paste(glue::glue("{conditionList}"),collapse="|"))) %>% 
    mutate(Experiment=str_remove(Experiment,paste(glue::glue("{conditionList}"),collapse="|")))
  data.ratio.AB <- data.ratio.AB %>% filter(pos1!=pos2)
  data.ratio.AB <- lapply(split(qDRIP_Group,qDRIP_Group$wig),function(x){
    data.ratio.AB %>% left_join(x,by=c("pos2"="gene_id")) %>% drop_na()
  }) %>% bind_rows()
  
  p1 <- data.ratio.AB %>% ggplot(aes(x=Experiment,y=score,fill=Group)) + geom_boxplot() + theme_classic() + facet_grid(~wig) +scale_y_log10()
}


