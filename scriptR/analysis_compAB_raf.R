
# LIBRARIES ------------------------------------------------------------

library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(HiTC)

# Genome 
seqInfohg19=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
seqlen=seqlengths(seqInfohg19)

# Parameters
# binSize=50000
# binSizeTxt="50kb"

# Exp
expe="manipA"

# Chromosomes
SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
Chr.V=paste0("chr",c(1:22,"X"))
seqlen=seqlengths(SeqInfo[Chr.V])

options(scipen=999)

# ATACseq
ATAC.GR=import.bed("/media/HDD_ROCHER/PROJET_INGE/4CSeq_PROCESS_VINCENT/AFTER_EVA_PROCESSING/PROCESSED_21_01_19/peak_calling/PEAKS/HC3HCBGX9_ATACseqA_DIvA_18s005247-1-1_Clouaire_lane1ATACseqADIvA.rmdups.bam_summits.bed")

for(binSize in c(1000000,500000,100000)){

  binSizeTxt <- glue::glue("{binSize/10e2}kb")
  
  
  # Load compAB
  res <- lapply(1:length(Chr.V),function(i){
    chrendi=seqlengths(SeqInfo[Chr.V[i]])
    binendi=ceiling(chrendi/binSize)
    bini.GR=GRanges(Chr.V[i],IRanges(breakInChunks(totalsize=seqlengths(SeqInfo)[i],chunksize=binSize)))
    ATACbini=countOverlaps(bini.GR,ATAC.GR)
    mes_files <- list.files("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/AB_comp/KR",pattern=glue::glue("eigen_{str_remove(Chr.V[i],'chr')}_.+_{binSizeTxt}.txt"),full.names=T)
    
    sres <- lapply(mes_files,function(one_file){
      compDIvAi=read.table(one_file,header=F)[,1]
      compDIvAi=scale(compDIvAi,scale=T)
      
      if(cor(compDIvAi,ATACbini,use="pairwise.complete.obs")<0){
        compDIvAi=-compDIvAi
      }
      
      
      binCompDIvA.GRi=bini.GR
      binCompDIvA.GRi$score=compDIvAi
      binCompDIvA.GRi$score[is.na(binCompDIvA.GRi$score)]=0
      binCompDIvA.GRi$score=round(binCompDIvA.GRi$score,2)
      
    }) %>% setNames(str_remove(str_remove(basename(mes_files),glue::glue("_{binSizeTxt}.txt")),glue::glue("eigen_{str_remove(Chr.V[i],'chr')}_KR_"))) %>% bind_cols()
    
    bini.GR %>% as_tibble() %>% cbind(sres)
    
  }) %>% bind_rows() 
  res <- res %>% dplyr::select(-width,-strand) %>% gather(key=Condition,value=score,-seqnames,-start,-end) %>% group_by(Condition) %>% nest()
  
  walk(1:nrow(res),function(i){
    Cond <- res$Condition[[i]]
    data <- res$data[[i]] %>% as_granges()
    export.bedGraph(data,con=glue::glue("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/AB_comp/ABcomp_{Cond}_{binSizeTxt}.bedGraph"))
  })
  
  
  
  
}

