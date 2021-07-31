require(HiTC)
require(tidyverse)
require(plyranges)
require(data.table)
require(Matrix)
require(BSgenome.Hsapiens.UCSC.hg19)
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

obserOE <- "OE"
expe <- c("HiC_D_DIvA","HiC_D_OHT","HiC_D_OHTATMi","HiC_D_OHTDNAPKi","HiC_D_OHTPARPi")
# Resolution
binSize=100000
binSizeTxt=paste0(binSize/1e3,"kb")

SeqInfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19)

my_chromosomes <- seqlevels(DSB174)

DIvA <- expe[1]
# OHT <- expe[4]
OHT <- expe[5]
res_pca <- lapply(my_chromosomes,function(Chr.V){
  message(Chr.V)
  chrendi=seqlengths(SeqInfo[Chr.V])
  binendi=ceiling(chrendi/binSize)
  
  fileCounti_DIvA=paste0("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/dump_",obserOE,"_KR_",DIvA,"_",Chr.V,"_",binSizeTxt,".txt.gz")
  # fileCounti_OHT=paste0("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/KR/dump_",obserOE,"_KR_",expe[3],"_",Chr.V,"_",Chr.V,"_",binSizeTxt,".txt.gz")
  fileCounti_OHT=paste0("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/dump_",obserOE,"_KR_",OHT,"_",Chr.V,"_",binSizeTxt,".txt.gz")
  
  
  
  data_mOHT <- read_Dumped_matrix(fileCounti_DIvA,binendi=binendi) %>% 
    as_tibble %>% 
    setNames(c("bin1","bin2","mOHT"))
  data_pOHT <- read_Dumped_matrix(fileCounti_OHT,binendi=binendi) %>% 
    as_tibble() %>%
    setNames(c("bin1","bin2","pOHT")) %>% inner_join(data_mOHT,by = c("bin1","bin2"))
  
  datai <- data_pOHT %>% mutate(ratio = log2(pOHT/mOHT))
  datai <- datai %>% mutate(ratio = ifelse(pOHT==0 & mOHT== 0,0,ratio))
  datai <- datai[is.finite(datai$ratio),]
  data.Mati=sparseMatrix(i=(datai$bin1/binSize)+1,j=(datai$bin2/binSize)+1,x=datai$ratio)
  data.MatSymi=data.Mati+t(data.Mati)
  diag(data.MatSymi)=diag(data.MatSymi)/2
  rm(data.Mati,datai)
  
  starti=seq(1,chrendi,by=binSize)
  endi=c(seq(binSize,chrendi,by=binSize),chrendi)
  chri.GR=GRanges(Chr.V,IRanges(starti,endi))
  names(chri.GR)=paste0("bin",1:length(chri.GR))
  
  HTC=HTCexp(data.MatSymi,chri.GR,chri.GR)
  
  cm <- HiTC:::sparseCor(intdata(HTC))
  cm[is.na(cm)] <- 0
  princp = princomp(cm)
  
  
  PCA_1 <-bind_cols(as_tibble(x_intervals(HTC)),as_tibble(loadings(princp)[])) %>%
    as_granges()%>% join_overlap_left(DSB174) 
  PCA_cor_DSB <- PCA_1 %>% mutate(DSB = ifelse(!is.na(name),1,0)) %>% as_tibble() %>% dplyr::select(-seqnames:-strand,-name,-score)
  
  
  PCA_cor_DSB <- PCA_cor_DSB %>% dplyr::select(-DSB) %>% apply(2,function(one_col){
    cc <- cor.test(one_col,PCA_cor_DSB$DSB)
    tibble(my_cor= cc$estimate,p.val = cc$p.value)
  }) %>% bind_rows() %>% mutate(name =colnames(PCA_cor_DSB)[-length(PCA_cor_DSB)])
  
  
  
  PCA_1 <- PCA_1 %>% as_tibble() %>%
    gather(key = PCs,value = value,-seqnames:-strand,-name:-score) %>% 
    left_join(PCA_cor_DSB,by = c("PCs"="name"))
  
  PCA_1 %>% filter(PCs %in% str_c("Comp",1:10,sep = ".")) 
  
}) %>% bind_rows()




res_pca <- res_pca  %>% as_granges()

# cc1 <- res_pca %>% filter_by_overlaps(DSB174 %>% anchor_center() %>% mutate(width = 10000)) %>% as_tibble() %>% group_by(seqnames) %>%
#   summarise(meanval = mean(value)) %>% filter(meanval < 0) %>% pull(seqnames) %>% as.character()
# 
# 
# res_pca <- res_pca %>% as_tibble() %>%split(.,.$seqnames)
# res_pca <- res_pca[-1]
# res_pca <- res_pca %>% map(function(x){
#   if(unique(x$seqnames) %in% cc1){
#     x %>% mutate(value = value * -1)
#   }else{
#     x
#   }
# })  



res_p <- lapply(res_pca[my_chromosomes],function(one_data){
  one_data <- one_data %>% filter(PCs %in% str_c("Comp",1,sep = "."))%>% mutate(my_cor = round(my_cor,3))%>% mutate(p.val = format.pval(p.val,3))%>% unite("PCs",PCs,my_cor,sep = " cor: ") %>% unite("PCs",PCs,p.val,sep = " p.val: ") 
  one_data %>%  ggplot(aes(x=start,y=value)) + geom_line() +
    facet_wrap(~PCs,scales="free",ncol=1) +
    geom_vline(data = filter(one_data,!is.na(name)),aes(xintercept=start),linetype ="dashed") +ggtitle(unique(one_data$seqnames))
})

pdf(glue::glue("../../results/{OHT}_test_pca_cor_alamano_100kb_1_with_reverse_signal_{binSizeTxt}.pdf"),height=4,width=12)
print(res_p)
dev.off()

##Extract to bigwig 
cc <- res_pca %>% map(filter,PCs %in% str_c("Comp",1,sep = ".")) %>% map(mutate,value=ifelse(my_cor < 0,value*-1,value)) %>% 
  bind_rows() %>% dplyr::select(seqnames:strand,value) %>% as_granges()
# cc <- res_pca %>% bind_rows() %>% filter(PCs == "Comp.1") %>% dplyr::select(seqnames:strand,value) %>% as_granges()

cc <- regioneR::filterChromosomes(cc,keep.chr=my_chromosomes)
seqlengths(cc) <- seqlengths(SeqInfo)[names(seqlengths(cc))]
cc.cov <- coverage(cc,weight="value")
export.bw(cc.cov,glue::glue("../../results/PC1_all_chr_log2ratio_{binSizeTxt}_{OHT}.bw"))


#Manip A / B 

expe=c("DIvA_manipA","OHT_manipA","DIvA_manipB","OHT_manipB")



# DIvA <- expe[1]
# OHT <- expe[2]
DIvA <- expe[3]
OHT <- expe[4]
res_pca <- lapply(my_chromosomes,function(Chr.V){
  message(Chr.V)
  chrendi=seqlengths(SeqInfo[Chr.V])
  binendi=ceiling(chrendi/binSize)
  
  fileCounti_DIvA=paste0("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/KR/dump_",obserOE,"_KR_",DIvA,"_",Chr.V,"_",binSizeTxt,".txt.gz")
  fileCounti_OHT=paste0("/home/rochevin/Documents/PROJET_INGE/HiC_Coline/data/KR/dump_",obserOE,"_KR_",OHT,"_",Chr.V,"_",binSizeTxt,".txt.gz")
  # fileCounti_OHT=paste0("/media/HDD_ROCHER/PROJET_THESE/PAPIER_COLINE_AUDE/data/HiC_D/KR/dump_",obserOE,"_KR_",expe[1],"_",Chr.V,"_",binSizeTxt,".txt.gz")
  
  
  data_mOHT <- read_Dumped_matrix(fileCounti_DIvA,binendi=binendi) %>% 
    as_tibble %>% 
    setNames(c("bin1","bin2","mOHT"))
  data_pOHT <- read_Dumped_matrix(fileCounti_OHT,binendi=binendi) %>% 
    as_tibble() %>%
    setNames(c("bin1","bin2","pOHT")) %>% inner_join(data_mOHT,by = c("bin1","bin2"))
  
  datai <- data_pOHT %>% mutate(ratio = log2(pOHT/mOHT))
  datai <- datai %>% mutate(ratio = ifelse(pOHT==0 & mOHT== 0,0,ratio))
  datai <- datai[is.finite(datai$ratio),]
  data.Mati=sparseMatrix(i=(datai$bin1/binSize)+1,j=(datai$bin2/binSize)+1,x=datai$ratio)
  data.MatSymi=data.Mati+t(data.Mati)
  diag(data.MatSymi)=diag(data.MatSymi)/2
  rm(data.Mati,datai)
  
  starti=seq(1,chrendi,by=binSize)
  endi=c(seq(binSize,chrendi,by=binSize),chrendi)
  chri.GR=GRanges(Chr.V,IRanges(starti,endi))
  names(chri.GR)=paste0("bin",1:length(chri.GR))
  
  HTC=HTCexp(data.MatSymi,chri.GR,chri.GR)
  
  cm <- HiTC:::sparseCor(intdata(HTC))
  cm[is.na(cm)] <- 0
  princp = princomp(cm)
  
  
  PCA_1 <-bind_cols(as_tibble(x_intervals(HTC)),as_tibble(loadings(princp)[])) %>%
    as_granges()%>% join_overlap_left(DSB174) 
  PCA_cor_DSB <- PCA_1 %>% mutate(DSB = ifelse(!is.na(name),1,0)) %>% as_tibble() %>% dplyr::select(-seqnames:-strand,-name,-score)
  
  
  PCA_cor_DSB <- PCA_cor_DSB %>% dplyr::select(-DSB) %>% apply(2,function(one_col){
    cc <- cor.test(one_col,PCA_cor_DSB$DSB)
    tibble(my_cor= cc$estimate,p.val = cc$p.value)
  }) %>% bind_rows() %>% mutate(name =colnames(PCA_cor_DSB)[-length(PCA_cor_DSB)])
  
  
  
  PCA_1 <- PCA_1 %>% as_tibble() %>%
    gather(key = PCs,value = value,-seqnames:-strand,-name:-score) %>% 
    left_join(PCA_cor_DSB,by = c("PCs"="name"))
  
  PCA_1 %>% filter(PCs %in% str_c("Comp",1:10,sep = ".")) 
  
}) %>% bind_rows()




res_pca <- res_pca  %>% as_granges()
res_pca <- res_pca %>% as_tibble() %>%split(.,.$seqnames)
res_pca <- res_pca[-1]
# cc1 <- res_pca %>% filter_by_overlaps(DSB174 %>% anchor_center() %>% mutate(width = 10000)) %>% as_tibble() %>% group_by(seqnames) %>%
#   summarise(meanval = mean(value)) %>% filter(meanval < 0) %>% pull(seqnames) %>% as.character()
# 
# 
# res_pca <- res_pca %>% as_tibble() %>%split(.,.$seqnames)
# res_pca <- res_pca[-1]
# res_pca <- res_pca %>% map(function(x){
#   if(unique(x$seqnames) %in% cc1){
#     x %>% mutate(value = value * -1)
#   }else{
#     x
#   }
# })  



res_p <- lapply(res_pca[my_chromosomes],function(one_data){
  one_data <- one_data %>% filter(PCs %in% str_c("Comp",1,sep = "."))%>% mutate(my_cor = round(my_cor,3))%>% mutate(p.val = format.pval(p.val,3))%>% unite("PCs",PCs,my_cor,sep = " cor: ") %>% unite("PCs",PCs,p.val,sep = " p.val: ") 
  one_data %>%  ggplot(aes(x=start,y=value)) + geom_line() +
    facet_wrap(~PCs,scales="free",ncol=1) +
    geom_vline(data = filter(one_data,!is.na(name)),aes(xintercept=start),linetype ="dashed") +ggtitle(unique(one_data$seqnames))
})

pdf(glue::glue("../../results/{OHT}_test_pca_cor_alamano_100kb_1_with_reverse_signal_{binSizeTxt}.pdf"),height=4,width=12)
print(res_p)
dev.off()

##Extract to bigwig 
cc <- res_pca %>% map(filter,PCs %in% str_c("Comp",1,sep = ".")) %>% map(mutate,value=ifelse(my_cor < 0,value*-1,value)) %>% 
  bind_rows() %>% dplyr::select(seqnames:strand,value) %>% as_granges()
# cc <- res_pca %>% bind_rows() %>% filter(PCs == "Comp.1") %>% dplyr::select(seqnames:strand,value) %>% as_granges()

cc <- regioneR::filterChromosomes(cc,keep.chr=my_chromosomes)
seqlengths(cc) <- seqlengths(SeqInfo)[names(seqlengths(cc))]
cc.cov <- coverage(cc,weight="value")
export.bw(cc.cov,glue::glue("../../results/PC1_all_chr_log2ratio_{binSizeTxt}_{OHT}.bw"))


#Correct what gaelle told me
# soucis avec le 18 -> manip D à l'envers OHT et parpi aussi surement
# manipD à l'envers chromosome 20
# 
# chr 21 Dnapki à l'envers et manip D OHT aussi 

files.PC1 <- list.files("../../results",full.names = T,pattern = "PC1_all_chr_log2ratio_100kb.+.bw") %>% 
  setNames(basename(.)) %>% 
  map(import.bw)

files.PC1[["PC1_all_chr_log2ratio_100kb_HiC_D_OHT.bw"]] <- files.PC1[["PC1_all_chr_log2ratio_100kb_HiC_D_OHT.bw"]] %>% mutate(score = ifelse(seqnames%in%c("chr18","chr20"),score*-1,score))
files.PC1[["PC1_all_chr_log2ratio_100kb_HiC_D_OHTDNAPKi.bw"]]  <- files.PC1[["PC1_all_chr_log2ratio_100kb_HiC_D_OHTDNAPKi.bw"]] %>% mutate(score = ifelse(seqnames%in%c("chr21"),score*-1,score))
files.PC1[["PC1_all_chr_log2ratio_100kb_HiC_D_OHTPARPi.bw"]] <- files.PC1[["PC1_all_chr_log2ratio_100kb_HiC_D_OHTPARPi.bw"]] %>% mutate(score = ifelse(seqnames%in%c("chr18"),score*-1,score))

for(one.f in names(files.PC1)){
  final_file <- glue::glue("../../results/{one.f}")
  export.bw(files.PC1[[one.f]],final_file)
}

#process trans matrix
#from PCA_trans_matrix.R 

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
res <- apply(Chr.V.combn,1,function(x){
  x1 <- x[1] %>% as.character()
  x2 <- x[2] %>% as.character()
  message(paste(x1,x2,sep="vs"))
  
  
  
  #Load the file
  file_DIvA <- glue::glue("{dir.out}dump_{obserOE}_KR_{DIvA}_{x1}_{x2}_{binSizeTxt}.txt.gz")
  #using fread from data.table
  datai=as.matrix(fread(cmd=paste0("zcat ",file_DIvA),sep='\t',header=F))
  data_DIvA=processingMATRIX(x1,x2,datai,my.ranges)  %>% 
    as_tibble %>% 
    setNames(c("bin1","bin2","mOHT"))
  #Load the file
  file_OHT <- glue::glue("{dir.out}dump_{obserOE}_KR_{OHT}_{x1}_{x2}_{binSizeTxt}.txt.gz")
  #using fread from data.table
  datai=as.matrix(fread(cmd=paste0("zcat ",file_OHT),sep='\t',header=F))
  data_OHT=processingMATRIX(x1,x2,datai,my.ranges) %>% 
    as_tibble %>% 
    setNames(c("bin1","bin2","pOHT"))
  
  datai <- data_DIvA %>% inner_join(data_OHT,by = c("bin1","bin2"))
  
  datai <- datai %>% mutate(ratio = log2(pOHT/mOHT))
  datai <- datai %>% mutate(ratio = ifelse(pOHT==0 & mOHT== 0,0,ratio))
  datai <- datai[is.finite(datai$ratio),]
  
}) %>% bind_rows() %>% as.matrix()

#Create the sparse matrix from the data frame pos1 pos2 score
data.Mati=sparseMatrix(i=res[,1],j=res[,2],x=res[,3])
#Create a symmetrical matrix with only the first half 
# Fill the other half with the transpose matrix
data.MatSymi=data.Mati+t(data.Mati)
# And divide by 2 the diag
diag(data.MatSymi)=diag(data.MatSymi)/2

rm(res)
cm <- HiTC:::sparseCor(data.Mati)
cm[is.na(cm)] <- 0

rm(data.Mati)
rm(data.MatSymi)
saveRDS(cm,glue::glue("../../results/Correlation_matrix_Trans_{OHT}_{binSizeTxt}.rds"))

require(bigmemory)
require(bigpca)
cm <- as.big.matrix(as.matrix(cm))
# princp = princomp(cm)
princp <- big.PCA(cm,pcs.to.keep = 1)

PCA_1 <-bind_cols(as_tibble(my.ranges),as_tibble(loadings(princp)[])) %>%
  as_granges()%>% join_overlap_left(DSB174) 
PCA_cor_DSB <- PCA_1 %>% mutate(DSB = ifelse(!is.na(name),1,0)) %>% as_tibble() %>% dplyr::select(-seqnames:-strand,-name,-score)


PCA_cor_DSB <- PCA_cor_DSB %>% dplyr::select(-DSB) %>% apply(2,function(one_col){
  cc <- cor.test(one_col,PCA_cor_DSB$DSB)
  tibble(my_cor= cc$estimate,p.val = cc$p.value)
}) %>% bind_rows() %>% mutate(name =colnames(PCA_cor_DSB)[-length(PCA_cor_DSB)])



PCA_1 <- PCA_1 %>% as_tibble() %>%
  gather(key = PCs,value = value,-seqnames:-strand,-name:-score) %>% 
  left_join(PCA_cor_DSB,by = c("PCs"="name"))

PCA_1 %>% filter(PCs %in% str_c("Comp",1:10,sep = ".")) 
res_pca <- res_pca  %>% as_granges()

# cc1 <- res_pca %>% filter_by_overlaps(DSB174 %>% anchor_center() %>% mutate(width = 10000)) %>% as_tibble() %>% group_by(seqnames) %>%
#   summarise(meanval = mean(value)) %>% filter(meanval < 0) %>% pull(seqnames) %>% as.character()
# 
# 
# res_pca <- res_pca %>% as_tibble() %>%split(.,.$seqnames)
# res_pca <- res_pca[-1]
# res_pca <- res_pca %>% map(function(x){
#   if(unique(x$seqnames) %in% cc1){
#     x %>% mutate(value = value * -1)
#   }else{
#     x
#   }
# })  



res_p <- lapply(res_pca[my_chromosomes],function(one_data){
  one_data <- one_data %>% filter(PCs %in% str_c("Comp",1,sep = "."))%>% mutate(my_cor = round(my_cor,3))%>% mutate(p.val = format.pval(p.val,3))%>% unite("PCs",PCs,my_cor,sep = " cor: ") %>% unite("PCs",PCs,p.val,sep = " p.val: ") 
  one_data %>%  ggplot(aes(x=start,y=value)) + geom_line() +
    facet_wrap(~PCs,scales="free",ncol=1) +
    geom_vline(data = filter(one_data,!is.na(name)),aes(xintercept=start),linetype ="dashed") +ggtitle(unique(one_data$seqnames))
})

pdf(glue::glue("../../results/{OHT}_test_pca_cor_alamano_100kb_1_with_reverse_signal_{binSizeTxt}.pdf"),height=4,width=12)
print(res_p)
dev.off()

##Extract to bigwig 
cc <- res_pca %>% map(filter,PCs %in% str_c("Comp",1,sep = ".")) %>% map(mutate,value=ifelse(my_cor < 0,value*-1,value)) %>% 
  bind_rows() %>% dplyr::select(seqnames:strand,value) %>% as_granges()
# cc <- res_pca %>% bind_rows() %>% filter(PCs == "Comp.1") %>% dplyr::select(seqnames:strand,value) %>% as_granges()

cc <- regioneR::filterChromosomes(cc,keep.chr=my_chromosomes)
seqlengths(cc) <- seqlengths(SeqInfo)[names(seqlengths(cc))]
cc.cov <- coverage(cc,weight="value")
export.bw(cc.cov,glue::glue("../../results/PC1_all_chr_log2ratio_{binSizeTxt}_{OHT}.bw"))




#Make files for manipA/B
dir.out <- "/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/DumpedMatrices/TRANS/"
expe <- c("DIvA_manipA","OHT_manipA","DIvA_manipB","OHT_manipB")
DIvA <- expe[3]
OHT <- expe[4]

obserOE <- "OE"
DSB174 <- read_bed("/home/rochevin/Documents/PROJET_THESE/PAPIER_COLINE_AUDE/data/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")

for(exp in expe){
  for(chr1 in chrs){
    for(chr2 in chrs){
      my_cmd <- glue::glue("java -jar /home/rochevin/Documents/PROJET_INGE/HiC_Coline/programs/juicer/juicer_tools.jar dump {obserOE} KR /mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC/hicfiles/{exp}/inter.hic {chr1} {chr2} BP {binSize} {dir.out}dump_{obserOE}_KR_{exp}_chr{chr1}_chr{chr2}_{binSizeTxt}.txt; gzip {dir.out}dump_{obserOE}_KR_{exp}_chr{chr1}_chr{chr2}_{binSizeTxt}.txt")
      
      system(my_cmd)
    }
  }
}

#For each combination of chromosome, do
res <- apply(Chr.V.combn,1,function(x){
  x1 <- x[1] %>% as.character()
  x2 <- x[2] %>% as.character()
  message(paste(x1,x2,sep="vs"))
  
  
  
  #Load the file
  file_DIvA <- glue::glue("{dir.out}dump_{obserOE}_KR_{DIvA}_{x1}_{x2}_{binSizeTxt}.txt.gz")
  #using fread from data.table
  datai=as.matrix(fread(cmd=paste0("zcat ",file_DIvA),sep='\t',header=F))
  data_DIvA=processingMATRIX(x1,x2,datai,my.ranges)  %>% 
    as_tibble %>% 
    setNames(c("bin1","bin2","mOHT"))
  #Load the file
  file_OHT <- glue::glue("{dir.out}dump_{obserOE}_KR_{OHT}_{x1}_{x2}_{binSizeTxt}.txt.gz")
  #using fread from data.table
  datai=as.matrix(fread(cmd=paste0("zcat ",file_OHT),sep='\t',header=F))
  data_OHT=processingMATRIX(x1,x2,datai,my.ranges) %>% 
    as_tibble %>% 
    setNames(c("bin1","bin2","pOHT"))
  
  datai <- data_DIvA %>% inner_join(data_OHT,by = c("bin1","bin2"))
  
  datai <- datai %>% mutate(ratio = log2(pOHT/mOHT))
  datai <- datai %>% mutate(ratio = ifelse(pOHT==0 & mOHT== 0,0,ratio))
  datai <- datai[is.finite(datai$ratio),]
  
}) %>% bind_rows() %>% as.matrix()

#Create the sparse matrix from the data frame pos1 pos2 score
data.Mati=sparseMatrix(i=res[,1],j=res[,2],x=res[,3])
#Create a symmetrical matrix with only the first half 
# Fill the other half with the transpose matrix
data.MatSymi=data.Mati+t(data.Mati)
# And divide by 2 the diag
diag(data.MatSymi)=diag(data.MatSymi)/2
#Return the range and the matrix sparse
res <- list(
  "range"=my.ranges,
  "matrix"=data.MatSymi
)

cm <- HiTC:::sparseCor(res$matrix)


