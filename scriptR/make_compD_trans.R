require(HiTC)
require(tidyverse)
require(plyranges)
require(data.table)
require(Matrix)
require(BSgenome.Hsapiens.UCSC.hg19)
setwd("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/DSB_compD/scriptR")
DSB174 <- read_bed("/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/ASIsites_hg19_174_clived_IQ1.5.bed")
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

read_Dumped_matrix <- function(fileCounti,binendi = NULL){
  datai=as.matrix(fread(paste0("zcat ",fileCounti),sep='\t',header=F))
  
  datai=datai[!is.na(datai[,3]),]
  datai=rbind(datai,c((binendi-1)*binSize,(binendi-1)*binSize,0))
  if(obserOE=="observed"){
    datai[,3]=round(datai[,3])
  }else{
    datai[,3]=round(datai[,3],2)
  }
  return(datai)
}

processingMATRIX <- function(x1,x2,datai,my.ranges,whichorder=getorder_for_our_Run){
  #For each chromosome, get chrom size and compute the last bin
  chrend1=seqlengths(SeqInfo[x1]) 
  binend1=ceiling(chrend1/binSize)
  chrend2=seqlengths(SeqInfo[x2])
  binend2=ceiling(chrend2/binSize)
  #Remove NA value (no interaction)
  datai=datai[!is.na(datai[,3]),]
  if(whichorder(x1,x2)){
    colnames.datai <- colnames(datai)
    datai <- datai[,c(2,1,3)]
    colnames(datai) <- colnames.datai
  }
  #Add the last bin (used to compute the SparseMatrix)
  datai=rbind(datai,c((binend1-1)*binSize,(binend2-1)*binSize,0))
  if(obserOE=="observed"){
    datai[,3]=round(datai[,3])
  }else{
    datai[,3]=round(datai[,3],2)
  }
  #Add 1 : conversion from base 0 to base 1
  datai[,1] <- datai[,1] + 1
  datai[,2] <- datai[,2] + 1
  # Get the index of the first chromosome
  indices_ranges <- my.ranges %>% filter(seqnames == x1) %>% as_tibble() %>% dplyr::select(start,idx)
  #Left join index and position
  datai <- datai %>% as_tibble() %>% left_join(indices_ranges,by = c("V1"="start")) %>%
    dplyr::rename(idx1 = idx) 
  # Get the index of the second chromosome
  indices_ranges <- my.ranges %>% filter(seqnames == x2) %>% as_tibble() %>% dplyr::select(start,idx)
  #Left join index and position
  datai <- datai %>% left_join(indices_ranges,by = c("V2"="start")) %>%
    dplyr::rename(idx2 = idx) 
  #Return the two index and the score
  datai %>% dplyr::select(idx1,idx2,V3)
}

# Resolution
binSize=100000
binSizeTxt=paste0(binSize/1e3,"kb")

SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)

my_chromosomes <- seqlevels(DSB174)
chrs <- seqlevels(DSB174)
Chr.V.combn <- crossing(col1=chrs,col2=chrs)

# Build the bins of all genome
my.ranges <- tileGenome(SeqInfo[chrs],tilewidth=binSize, cut.last.tile.in.chrom=TRUE)
# my.ranges <- my.ranges %>% as_tibble() %>% group_by(seqnames) %>% filter(end != max(end)) %>% as_granges()
# Create an index
my.ranges <- my.ranges %>% mutate(idx = 1:length(my.ranges))

getorder_for_our_Run <- function(x1,x2){
  my_order <- paste0("chr",c(1:7,"X",8:18,20,"Y",19,22,21))
  grep(paste("^",x1,"$", sep=""),my_order) > grep(paste("^",x2,"$", sep=""),my_order)
}

# Resolution
binSize=100000
binSizeTxt=paste0(binSize/1e3,"kb")


#Make files for manipD
dir.out <- "/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/TRANS/KR/"
expe <- c("HiC_D_DIvA","HiC_D_OHT","HiC_D_OHTATMi","HiC_D_OHTDNAPKi","HiC_D_OHTPARPi")
DIvA <- expe[1]
OHT <- expe[4]

obserOE <- "OE"

#For each combination of chromosome, do
# res <- apply(Chr.V.combn,1,function(x){
#   x1 <- x[1] %>% as.character()
#   x2 <- x[2] %>% as.character()
#   message(paste(x1,x2,sep="vs"))
#   
#   
#   
#   #Load the file
#   file_DIvA <- glue::glue("{dir.out}dump_{obserOE}_KR_{DIvA}_{x1}_{x2}_{binSizeTxt}.txt.gz")
#   #using fread from data.table
#   datai=as.matrix(fread(cmd=paste0("zcat ",file_DIvA),sep='\t',header=F))
#   data_DIvA=processingMATRIX(x1,x2,datai,my.ranges)  %>% 
#     as_tibble %>% 
#     setNames(c("bin1","bin2","mOHT"))
#   #Load the file
#   file_OHT <- glue::glue("{dir.out}dump_{obserOE}_KR_{OHT}_{x1}_{x2}_{binSizeTxt}.txt.gz")
#   #using fread from data.table
#   datai=as.matrix(fread(cmd=paste0("zcat ",file_OHT),sep='\t',header=F))
#   data_OHT=processingMATRIX(x1,x2,datai,my.ranges) %>% 
#     as_tibble %>% 
#     setNames(c("bin1","bin2","pOHT"))
#   
#   datai <- data_DIvA %>% inner_join(data_OHT,by = c("bin1","bin2"))
#   
#   datai <- datai %>% mutate(ratio = log2(pOHT/mOHT))
#   datai <- datai %>% mutate(ratio = ifelse(pOHT==0 & mOHT== 0,0,ratio))
#   datai <- datai[is.finite(datai$ratio),]
#   
# }) %>% bind_rows() %>% as.matrix()
# 
# #Create the sparse matrix from the data frame pos1 pos2 score
# data.Mati=sparseMatrix(i=res[,1],j=res[,2],x=res[,3])
# #Create a symmetrical matrix with only the first half 
# # Fill the other half with the transpose matrix
# data.MatSymi=data.Mati+t(data.Mati)
# # And divide by 2 the diag
# diag(data.MatSymi)=diag(data.MatSymi)/2
# 
# rm(res)
# cm <- HiTC:::sparseCor(data.Mati)
# cm[is.na(cm)] <- 0
# 
# rm(data.Mati)
# rm(data.MatSymi)
# saveRDS(cm,glue::glue("../../results/Correlation_matrix_Trans_{OHT}_{binSizeTxt}.rds"))
# saveRDS(cm,glue::glue("../../results/Correlation_matrix_Trans_{OHT}_{binSizeTxt}.rds"))
cm <- readRDS(glue::glue("../../results/Correlation_matrix_Trans_{OHT}_{binSizeTxt}.rds"))
require(bigmemory)
require(bigpca)
cm <- as.big.matrix(as.matrix(cm))
# princp = princomp(cm)
princp <- big.PCA(cm,pcs.to.keep = 1)
saveRDS(princp,glue::glue("../../results/princp_Trans_{OHT}_{binSizeTxt}.rds"))