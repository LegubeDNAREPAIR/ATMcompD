require(tidyverse)
require(cowplot)
require(plyranges)
require(rtracklayer)
require(BSgenome.Hsapiens.UCSC.hg19.masked)
seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)


require(tidyverse)
require(cowplot)
require(plyranges)
require(rtracklayer)
require(BSgenome.Hsapiens.UCSC.hg19.masked)
seqlens <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
#TRANS ANALYSIS

#BOXPLOT SUR LES PROFILES 4C-seq
replaceName <- list(
  "4C-Seq_C_Legube_1_fastq_after_trimming"="C_mOHT"           
  ,"4C-Seq_C_Legube_2_fastq_after_trimming"="C_pOHT"           
  ,"4C-Seq_C_Legube_3_fastq_after_trimming"="C_pOHTpDNAPKi"           
  ,"4C-Seq_C_Legube_4_fastq_after_trimming"="C_pOHTpATMi"           
  ,"4C-seq_E-Legube_1_S1_all_R1_001_cutadapt"="E_mOHT"         
  ,"4C-seq_E-Legube_2_S2_all_R1_001_cutadapt"="E_pOHT24h"         
  ,"4C-Seq-D-Legube-Banque1_S1_all_R1_001_cutadapt"="D_mOHT"   
  ,"4C-Seq-D-Legube-Banque2_S2_all_R1_001_cutadapt"="D_pOHT"   
  ,"4C-Seq-D-Legube-Banque3_S3_all_R1_001_cutadapt"="D_pOHTpDNAPKi"  
  ,"4C-Seq-D-Legube-Banque4_S4_all_R1_001_cutadapt"="D_pOHTpATMi"   
  ,"4C-Seq-Dbis-Legube-Banque1_S1_all_R1_001_cutadapt"="Dbis_mOHT"
  ,"4C-Seq-Dbis-Legube-Banque2_S2_all_R1_001_cutadapt"="Dbis_pOHT"
  ,"4C-Seq-Dbis-Legube-Banque3_S3_all_R1_001_cutadapt"="Dbis_pOHTpDNAPKi"
  ,"4C-Seq-Dbis-Legube-Banque4_S4_all_R1_001_cutadapt"="Dbis_pOHTpATMi"
  ,"4C-seq-F-Legube_1_fastq_after_trimming"  ="F_G1mOHT"         
  ,"4C-seq-F-Legube_2_fastq_after_trimming" ="F_G1pOHT"          
  ,"4C-seq-F-Legube_3_fastq_after_trimming" ="F_G2mOHT"          
  ,"4C-seq-F-Legube_4_fastq_after_trimming" ="F_G2pOHT"          
  ,"4C-seq-F-Legube_5_fastq_after_trimming" ="F_SmOHT"          
  ,"4C-seq-F-Legube_6_fastq_after_trimming"  ="F_SpOHT"         
  ,"4Cseq-LBCMCP-index1_S1_all_R1_001_cutadapt"="A_mOHT"       
  ,"4Cseq-LBCMCP-index2_S2_all_R1_001_cutadapt"="A_pOHT"       
  ,"LEGU-7_1_fastq_after_trimming"="G_siCTRLmOHT"                    
  ,"LEGU-7_2_fastq_after_trimming" ="G_siCTRLpOHT"                   
  ,"LEGU-7_3_fastq_after_trimming"="G_siSCC1mOHT"                    
  ,"LEGU-7_4_fastq_after_trimming"="G_siSCC1pOHT"                     
  ,"LEGU-7_5_fastq_after_trimming"="G_siCTCFmOHT"                    
  ,"LEGU-7_6_fastq_after_trimming"="G_siCTCFpOHT"                    
  ,"LEGU-8_1_fastq_after_trimming"="H_siCTRLmOHT"                    
  ,"LEGU-8_2_fastq_after_trimming"="H_siCTRLpOHT"                    
  ,"LEGU-8_3_fastq_after_trimming"="H_siSCC1mOHT"                    
  ,"LEGU-8_4_fastq_after_trimming"="H_siSCC1pOHT"                     
  ,"LEGU-8_5_fastq_after_trimming"="H_siCTCFmOHT"                    
  ,"LEGU-8_6_fastq_after_trimming"="H_siCTCFpOHT" 
)

bless80 <- "/home/rochevin/Documents/PROJET_THESE/CLASSIF_HR_NHEJ/data/BED/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed" %>% import.bed()
path4C <- "/home/rochevin/Documents/PROJET_INGE/4CSeq_PROCESS_VINCENT/AFTER_EVA_PROCESSING/4C_ALN_TRANS/BIGWIG_NORMALIZED"
wigs.4C <- list.files(path4C)


Get1val <- function(my.wigs,one.w,x){
  lapply(split(x,droplevels(seqnames(x))),function(zz){
    message(unique(as.character(seqnames(zz))))
    cov <- one.w[[unique(as.character(seqnames(zz)))]]
    score <- Views( cov, start = start(zz), end = end(zz) ) %>% sum()
    tibble(wig = my.wigs,value = score,rowname = zz$name)
  }) %>% bind_rows()
}
to.plot.4C <- sres <- lapply(c(100000,500000,1000000),function(bin){
  values <-  bless80 %>% anchor_center() %>% mutate(width=bin)
  lapply(wigs.4C,function(wig){
    message(wig)
    one.w <- import.bw(str_c(path4C,wig,sep="/"),as="RleList")
    File <- replaceName[[str_remove(wig,"_chr.*.bw")]]
    if(is.null(File))
      return(NULL)
    vp <- str_extract(wig,"chr[0-9A-Z]+_[0-9]+-[0-9]+") %>% str_replace("_",":")
    x <- values %>% filter_by_non_overlaps(GRanges(vp))
    Get1val(File,one.w,x) %>% mutate(binsize = bin) %>% mutate(viewpoint = vp)
  }) %>% bind_rows()
  
}) %>% bind_rows()

to.plot.4C <- to.plot.4C %>% separate(wig,into = c("Manip","Condition"),sep= "_") 


#Barplot dataset all DSBs
pdf("BOXPLOT_TRANS/TRANS_bless80_on4C_seq_data_barplot_mean_all_dsb.pdf",height=12,width=16)
for(bin in unique(to.plot.4C$binsize)){
  for(onemanip in unique(to.plot.4C$Manip)){
    subdf <- to.plot.4C %>% filter(binsize == bin,viewpoint %in% c("chr20:30946314-30947710","chr1:89455867-89456712","chr17:45759770-45760603","chr17:57168614-57169531","chr21:33251469-33252587"),Manip == onemanip)
    
    p <- subdf %>%
      ggplot(aes(x=Condition,y=value,fill = Condition)) +  
      stat_summary(fun.y=mean,position=position_dodge(width=0.95),geom="bar")+
      stat_summary(fun.data=mean_se,position=position_dodge(0.95),geom="errorbar") + 
      coord_flip()+
      facet_grid(viewpoint~Manip,scales="free_x") +
      theme_minimal() + theme(legend.position = "none",axis.text.x = element_text(angle = 90)) +ggtitle(bin)
    print(p)
  }
  
  
}
dev.off()

#Update au 28/09/2020 pvals de la familia pour le papier de Coline
ccpval <- to.plot.4C %>% filter(viewpoint %in% c("chr20:30946314-30947710","chr1:89455867-89456712","chr17:45759770-45760603","chr17:57168614-57169531","chr21:33251469-33252587"),Manip == "D",binsize == 1e+06)
ccpval %>% filter(Condition %in% c("mOHT","pOHT")) %>% spread(key = Condition,value=value)  %>%
  group_by(viewpoint) %>% nest() %>%
  mutate(pvsmOHT = map_dbl(data,function(x){
    wilcox.test(x$pOHT,x$mOHT,paired=T)$p.value
  })) %>% dplyr::select(-data)



replaceName.compare <- list(
  "4C-Seq_A_pvsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
  ,"4C-Seq_A_pvsmOHT_Legube_chr17_57168614-57169531.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
  ,"4C-Seq_A_pvsmOHT_Legube_chr17_57578597-57579677.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
  ,"4C-Seq_A_pvsmOHT_Legube_chr1_80484021-80485207.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
  ,"4C-Seq_A_pvsmOHT_Legube_chr1_88985251-88986220.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
  ,"4C-Seq_A_pvsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
  ,"4C-Seq_A_pvsmOHT_Legube_chr1_89659070-89659811.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
  ,"4C-Seq_A_pvsmOHT_Legube_chr1_90042784-90045873.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
  ,"4C-Seq_A_pvsmOHT_Legube_chr21_32863772-32864618.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
  ,"4C-Seq_A_pvsmOHT_Legube_chr21_33251469-33252587.bw_pvsmOHT"="4C-Seq_A_pvsmOHT"
  ,"4C-Seq_C_pATMivsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_C_pATMivsmOHT"
  ,"4C-Seq_C_pATMivsmOHT_Legube_chr17_57168614-57169531_pvsmOHT"="4C-Seq_C_pATMivsmOHT"
  ,"4C-Seq_C_pATMivsmOHT_Legube_chr1_89455867-89456712_pvsmOHT"="4C-Seq_C_pATMivsmOHT"
  ,"4C-Seq_C_pATMivsmOHT_Legube_chr21_33251469-33252587_pvsmOHT"="4C-Seq_C_pATMivsmOHT"
  ,"4C-Seq_C_pDNAPKivsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_C_pDNAPKivsmOHT"
  ,"4C-Seq_C_pDNAPKivsmOHT_Legube_chr17_57168614-57169531_pvsmOHT"="4C-Seq_C_pDNAPKivsmOHT"
  ,"4C-Seq_C_pDNAPKivsmOHT_Legube_chr1_89455867-89456712_pvsmOHT"="4C-Seq_C_pDNAPKivsmOHT"
  ,"4C-Seq_C_pDNAPKivsmOHT_Legube_chr21_33251469-33252587_pvsmOHT"="4C-Seq_C_pDNAPKivsmOHT"
  ,"4C-Seq_C_pvsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_C_pvsmOHT"
  ,"4C-Seq_C_pvsmOHT_Legube_chr17_57168614-57169531_pvsmOHT"="4C-Seq_C_pvsmOHT"
  ,"4C-Seq_C_pvsmOHT_Legube_chr1_89455867-89456712_pvsmOHT"="4C-Seq_C_pvsmOHT"
  ,"4C-Seq_C_pvsmOHT_Legube_chr21_33251469-33252587_pvsmOHT"="4C-Seq_C_pvsmOHT"
  ,"4C-Seq_Dbis_pATMivsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_Dbis_pATMivsmOHT"
  ,"4C-Seq_Dbis_pATMivsmOHT_Legube_chr20_30946314-30947710.bw_pvsmOHT"="4C-Seq_Dbis_pATMivsmOHT"
  ,"4C-Seq_Dbis_pDNAPKivsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_Dbis_pDNAPKivsmOHT"
  ,"4C-Seq_Dbis_pDNAPKivsmOHT_Legube_chr20_30946314-30947710.bw_pvsmOHT"="4C-Seq_Dbis_pDNAPKivsmOHT"
  ,"4C-Seq_Dbis_pvsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_Dbis_pvsmOHT"
  ,"4C-Seq_Dbis_pvsmOHT_Legube_chr20_30946314-30947710.bw_pvsmOHT"="4C-Seq_Dbis_pvsmOHT"
  ,"4C-Seq_D_pATMivsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_D_pATMivsmOHT"
  ,"4C-Seq_D_pATMivsmOHT_Legube_chr17_57168614-57169531.bw_pvsmOHT"="4C-Seq_D_pATMivsmOHT"
  ,"4C-Seq_D_pATMivsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_D_pATMivsmOHT"
  ,"4C-Seq_D_pATMivsmOHT_Legube_chr21_33251469-33252587.bw_pvsmOHT"="4C-Seq_D_pATMivsmOHT"
  ,"4C-Seq_D_pDNAPKivsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_D_pDNAPKivsmOHT"
  ,"4C-Seq_D_pDNAPKivsmOHT_Legube_chr17_57168614-57169531.bw_pvsmOHT"="4C-Seq_D_pDNAPKivsmOHT"
  ,"4C-Seq_D_pDNAPKivsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_D_pDNAPKivsmOHT"
  ,"4C-Seq_D_pDNAPKivsmOHT_Legube_chr21_33251469-33252587.bw_pvsmOHT"="4C-Seq_D_pDNAPKivsmOHT"
  ,"4C-Seq_D_pvsmOHT_Legube_chr17_45759770-45760603.bw_pvsmOHT"="4C-Seq_D_pvsmOHT"
  ,"4C-Seq_D_pvsmOHT_Legube_chr17_57168614-57169531.bw_pvsmOHT"="4C-Seq_D_pvsmOHT"
  ,"4C-Seq_D_pvsmOHT_Legube_chr1_89455867-89456712.bw_pvsmOHT"="4C-Seq_D_pvsmOHT"
  ,"4C-Seq_D_pvsmOHT_Legube_chr21_33251469-33252587.bw_pvsmOHT"="4C-Seq_D_pvsmOHT"
  ,"4C-seq_E-Legube_p24vsmOHT_chr17_45759770-45760603_pvsmOHT"="4C-seq_E_p24vsmOHT"
  ,"4C-seq_E-Legube_p24vsmOHT_chr17_57168614-57169531_pvsmOHT"="4C-seq_E_p24vsmOHT"
  ,"4C-seq_E-Legube_p24vsmOHT_chr1_89455867-89456712_pvsmOHT"="4C-seq_E_p24vsmOHT"
  ,"4C-seq_E-Legube_p24vsmOHT_chr1_90292860-90295480_pvsmOHT"="4C-seq_E_p24vsmOHT"
  ,"4C-seq_E-Legube_p24vsmOHT_chr20_30946314-30947710_pvsmOHT"="4C-seq_E_p24vsmOHT"
  ,"4C-seq_E-Legube_p24vsmOHT_chr21_33251469-33252587_pvsmOHT"="4C-seq_E_p24vsmOHT"
  ,"4C-seq-F-Legube_pvsmOHT_G1_chr17_45759770-45760603_pvsmOHT"="4C-seq_F_pvsmOHT_G1"
  ,"4C-seq-F-Legube_pvsmOHT_G1_chr17_57168614-57169531_pvsmOHT"="4C-seq_F_pvsmOHT_G1"
  ,"4C-seq-F-Legube_pvsmOHT_G1_chr1_89455867-89456712_pvsmOHT"="4C-seq_F_pvsmOHT_G1"
  ,"4C-seq-F-Legube_pvsmOHT_G1_chr20_30946314-30947710_pvsmOHT"="4C-seq_F_pvsmOHT_G1"
  ,"4C-seq-F-Legube_pvsmOHT_G1_chr21_33251469-33252587_pvsmOHT"="4C-seq_F_pvsmOHT_G1"
  ,"4C-seq-F-Legube_pvsmOHT_G2_chr17_45759770-45760603_pvsmOHT"="4C-seq_F_pvsmOHT_G2"
  ,"4C-seq-F-Legube_pvsmOHT_G2_chr17_57168614-57169531_pvsmOHT"="4C-seq_F_pvsmOHT_G2"
  ,"4C-seq-F-Legube_pvsmOHT_G2_chr1_89455867-89456712_pvsmOHT"="4C-seq_F_pvsmOHT_G2"
  ,"4C-seq-F-Legube_pvsmOHT_G2_chr20_30946314-30947710_pvsmOHT"="4C-seq_F_pvsmOHT_G2"
  ,"4C-seq-F-Legube_pvsmOHT_G2_chr21_33251469-33252587_pvsmOHT"="4C-seq_F_pvsmOHT_G2"
  ,"4C-seq-F-Legube_pvsmOHT_S_chr17_45759770-45760603_pvsmOHT"="4C-seq_F_pvsmOHT_S"
  ,"4C-seq-F-Legube_pvsmOHT_S_chr17_57168614-57169531_pvsmOHT"="4C-seq_F_pvsmOHT_S"
  ,"4C-seq-F-Legube_pvsmOHT_S_chr1_89455867-89456712_pvsmOHT"="4C-seq_F_pvsmOHT_S"
  ,"4C-seq-F-Legube_pvsmOHT_S_chr20_30946314-30947710_pvsmOHT"="4C-seq_F_pvsmOHT_S"
  ,"4C-seq-F-Legube_pvsmOHT_S_chr21_33251469-33252587_pvsmOHT"="4C-seq_F_pvsmOHT_S"
  ,"LEGU-G_pvsmOHT_siCTCF_chr17_45759770-45760603_pvsmOHT"="4C-seq_G_pvsmOHT_siCTCF"
  ,"LEGU-G_pvsmOHT_siCTCF_chr17_57168614-57169531_pvsmOHT"="4C-seq_G_pvsmOHT_siCTCF"
  ,"LEGU-G_pvsmOHT_siCTCF_chr1_89455867-89456712_pvsmOHT"="4C-seq_G_pvsmOHT_siCTCF"
  ,"LEGU-G_pvsmOHT_siCTCF_chr20_30946314-30947710_pvsmOHT"="4C-seq_G_pvsmOHT_siCTCF"
  ,"LEGU-G_pvsmOHT_siCTCF_chr21_33251469-33252587_pvsmOHT"="4C-seq_G_pvsmOHT_siCTCF"
  ,"LEGU-G_pvsmOHT_siCTRL_chr17_45759770-45760603_pvsmOHT"="4C-seq_G_pvsmOHT_siCTRL"
  ,"LEGU-G_pvsmOHT_siCTRL_chr17_57168614-57169531_pvsmOHT"="4C-seq_G_pvsmOHT_siCTRL"
  ,"LEGU-G_pvsmOHT_siCTRL_chr1_89455867-89456712_pvsmOHT"="4C-seq_G_pvsmOHT_siCTRL"
  ,"LEGU-G_pvsmOHT_siCTRL_chr20_30946314-30947710_pvsmOHT"="4C-seq_G_pvsmOHT_siCTRL"
  ,"LEGU-G_pvsmOHT_siCTRL_chr21_33251469-33252587_pvsmOHT"="4C-seq_G_pvsmOHT_siCTRL"
  ,"LEGU-G_pvsmOHT_siSCC1_chr17_45759770-45760603_pvsmOHT"="4C-seq_G_pvsmOHT_siSCC1"
  ,"LEGU-G_pvsmOHT_siSCC1_chr17_57168614-57169531_pvsmOHT"="4C-seq_G_pvsmOHT_siSCC1"
  ,"LEGU-G_pvsmOHT_siSCC1_chr1_89455867-89456712_pvsmOHT"="4C-seq_G_pvsmOHT_siSCC1"
  ,"LEGU-G_pvsmOHT_siSCC1_chr20_30946314-30947710_pvsmOHT"="4C-seq_G_pvsmOHT_siSCC1"
  ,"LEGU-G_pvsmOHT_siSCC1_chr21_33251469-33252587_pvsmOHT"="4C-seq_G_pvsmOHT_siSCC1"
  ,"LEGU-H_pvsmOHT_siCTCF_chr17_45759770-45760603_pvsmOHT"="4C-seq_H_pvsmOHT_siCTCF"
  ,"LEGU-H_pvsmOHT_siCTCF_chr17_57168614-57169531_pvsmOHT"="4C-seq_H_pvsmOHT_siCTCF"
  ,"LEGU-H_pvsmOHT_siCTCF_chr1_89455867-89456712_pvsmOHT"="4C-seq_H_pvsmOHT_siCTCF"
  ,"LEGU-H_pvsmOHT_siCTCF_chr20_30946314-30947710_pvsmOHT"="4C-seq_H_pvsmOHT_siCTCF"
  ,"LEGU-H_pvsmOHT_siCTCF_chr21_33251469-33252587_pvsmOHT"="4C-seq_H_pvsmOHT_siCTCF"
  ,"LEGU-H_pvsmOHT_siCTRL_chr17_45759770-45760603_pvsmOHT"="4C-seq_H_pvsmOHT_siCTRL"
  ,"LEGU-H_pvsmOHT_siCTRL_chr17_57168614-57169531_pvsmOHT"="4C-seq_H_pvsmOHT_siCTRL"
  ,"LEGU-H_pvsmOHT_siCTRL_chr1_89455867-89456712_pvsmOHT"="4C-seq_H_pvsmOHT_siCTRL"
  ,"LEGU-H_pvsmOHT_siCTRL_chr20_30946314-30947710_pvsmOHT"="4C-seq_H_pvsmOHT_siCTRL"
  ,"LEGU-H_pvsmOHT_siCTRL_chr21_33251469-33252587_pvsmOHT"="4C-seq_H_pvsmOHT_siCTRL"
  ,"LEGU-H_pvsmOHT_siSCC1_chr17_45759770-45760603_pvsmOHT"="4C-seq_H_pvsmOHT_siSCC1"
  ,"LEGU-H_pvsmOHT_siSCC1_chr17_57168614-57169531_pvsmOHT"="4C-seq_H_pvsmOHT_siSCC1"
  ,"LEGU-H_pvsmOHT_siSCC1_chr1_89455867-89456712_pvsmOHT"="4C-seq_H_pvsmOHT_siSCC1"
  ,"LEGU-H_pvsmOHT_siSCC1_chr20_30946314-30947710_pvsmOHT"="4C-seq_H_pvsmOHT_siSCC1"
  ,"LEGU-H_pvsmOHT_siSCC1_chr21_33251469-33252587_pvsmOHT"="4C-seq_H_pvsmOHT_siSCC1"
  ,"4C-Seq_G_siCTCFvsCTRL_mOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_mOHT"
  ,"4C-Seq_G_siCTCFvsCTRL_mOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_mOHT"
  ,"4C-Seq_G_siCTCFvsCTRL_mOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_mOHT"
  ,"4C-Seq_G_siCTCFvsCTRL_mOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_mOHT"
  ,"4C-Seq_G_siCTCFvsCTRL_mOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_mOHT"
  ,"4C-Seq_G_siCTCFvsCTRL_pOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_pOHT"
  ,"4C-Seq_G_siCTCFvsCTRL_pOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_pOHT"
  ,"4C-Seq_G_siCTCFvsCTRL_pOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_pOHT"
  ,"4C-Seq_G_siCTCFvsCTRL_pOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_pOHT"
  ,"4C-Seq_G_siCTCFvsCTRL_pOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_G_siCTCFvsCTRL_pOHT"
  ,"4C-Seq_G_siRad21vsCTRL_mOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_mOHT"
  ,"4C-Seq_G_siRad21vsCTRL_mOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_mOHT"
  ,"4C-Seq_G_siRad21vsCTRL_mOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_mOHT"
  ,"4C-Seq_G_siRad21vsCTRL_mOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_mOHT"
  ,"4C-Seq_G_siRad21vsCTRL_mOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_mOHT"
  ,"4C-Seq_G_siRad21vsCTRL_pOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_pOHT"
  ,"4C-Seq_G_siRad21vsCTRL_pOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_pOHT"
  ,"4C-Seq_G_siRad21vsCTRL_pOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_pOHT"
  ,"4C-Seq_G_siRad21vsCTRL_pOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_pOHT"
  ,"4C-Seq_G_siRad21vsCTRL_pOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_G_siRad21vsCTRL_pOHT"
  ,"4C-Seq_H_siCTCFvsCTRL_mOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_mOHT"
  ,"4C-Seq_H_siCTCFvsCTRL_mOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_mOHT"
  ,"4C-Seq_H_siCTCFvsCTRL_mOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_mOHT"
  ,"4C-Seq_H_siCTCFvsCTRL_mOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_mOHT"
  ,"4C-Seq_H_siCTCFvsCTRL_mOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_mOHT"
  ,"4C-Seq_H_siCTCFvsCTRL_pOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_pOHT"
  ,"4C-Seq_H_siCTCFvsCTRL_pOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_pOHT"
  ,"4C-Seq_H_siCTCFvsCTRL_pOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_pOHT"
  ,"4C-Seq_H_siCTCFvsCTRL_pOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_pOHT"
  ,"4C-Seq_H_siCTCFvsCTRL_pOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_H_siCTCFvsCTRL_pOHT"
  ,"4C-Seq_H_siRad21vsCTRL_mOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_mOHT"
  ,"4C-Seq_H_siRad21vsCTRL_mOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_mOHT"
  ,"4C-Seq_H_siRad21vsCTRL_mOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_mOHT"
  ,"4C-Seq_H_siRad21vsCTRL_mOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_mOHT"
  ,"4C-Seq_H_siRad21vsCTRL_mOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_mOHT"
  ,"4C-Seq_H_siRad21vsCTRL_pOHT_chr17_45759770-45760603_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_pOHT"
  ,"4C-Seq_H_siRad21vsCTRL_pOHT_chr17_57168614-57169531_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_pOHT"
  ,"4C-Seq_H_siRad21vsCTRL_pOHT_chr1_89455867-89456712_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_pOHT"
  ,"4C-Seq_H_siRad21vsCTRL_pOHT_chr20_30946314-30947710_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_pOHT"
  ,"4C-Seq_H_siRad21vsCTRL_pOHT_chr21_33251469-33252587_pvsmOHT"="4C-Seq_H_siRad21vsCTRL_pOHT"
  #Derniers plots pour les dernières données (31 oct 19)
  
  ,"4C-Seq-I_pvsmOHT_G1_chr17_45759770-45760603_pvsmOHT"="4C-seq_I_pvsmOHT_G1"
  ,"4C-Seq-I_pvsmOHT_G1_chr17_57168614-57169531_pvsmOHT"="4C-seq_I_pvsmOHT_G1"
  ,"4C-Seq-I_pvsmOHT_G1_chr1_89455867-89456712_pvsmOHT"="4C-seq_I_pvsmOHT_G1"
  ,"4C-Seq-I_pvsmOHT_G1_chr20_30946314-30947710_pvsmOHT"="4C-seq_I_pvsmOHT_G1"
  ,"4C-Seq-I_pvsmOHT_G1_chr21_33251469-33252587_pvsmOHT"="4C-seq_I_pvsmOHT_G1"
  ,"4C-Seq-I_pvsmOHT_G2_chr17_45759770-45760603_pvsmOHT"="4C-seq_I_pvsmOHT_G2"
  ,"4C-Seq-I_pvsmOHT_G2_chr17_57168614-57169531_pvsmOHT"="4C-seq_I_pvsmOHT_G2"
  ,"4C-Seq-I_pvsmOHT_G2_chr1_89455867-89456712_pvsmOHT"="4C-seq_I_pvsmOHT_G2"
  ,"4C-Seq-I_pvsmOHT_G2_chr20_30946314-30947710_pvsmOHT"="4C-seq_I_pvsmOHT_G2"
  ,"4C-Seq-I_pvsmOHT_G2_chr21_33251469-33252587_pvsmOHT"="4C-seq_I_pvsmOHT_G2"
  ,"4C-Seq-I_pvsmOHT_S_chr17_45759770-45760603_pvsmOHT"="4C-seq_I_pvsmOHT_S"
  ,"4C-Seq-I_pvsmOHT_S_chr17_57168614-57169531_pvsmOHT"="4C-seq_I_pvsmOHT_S"
  ,"4C-Seq-I_pvsmOHT_S_chr1_89455867-89456712_pvsmOHT"="4C-seq_I_pvsmOHT_S"
  ,"4C-Seq-I_pvsmOHT_S_chr20_30946314-30947710_pvsmOHT"="4C-seq_I_pvsmOHT_S"
  ,"4C-Seq-I_pvsmOHT_S_chr21_33251469-33252587_pvsmOHT"="4C-seq_I_pvsmOHT_S"
  ,"4C-Seq-K_pvsmOHT_siCTCF_chr17_45759770-45760603_pvsmOHT"="4C-seq_K_pvsmOHT_siCTCF"
  ,"4C-Seq-K_pvsmOHT_siCTCF_chr17_57168614-57169531_pvsmOHT"="4C-seq_K_pvsmOHT_siCTCF"
  ,"4C-Seq-K_pvsmOHT_siCTCF_chr1_89455867-89456712_pvsmOHT"="4C-seq_K_pvsmOHT_siCTCF"
  ,"4C-Seq-K_pvsmOHT_siCTCF_chr20_30946314-30947710_pvsmOHT"="4C-seq_K_pvsmOHT_siCTCF"
  ,"4C-Seq-K_pvsmOHT_siCTCF_chr21_33251469-33252587_pvsmOHT"="4C-seq_K_pvsmOHT_siCTCF"
  ,"4C-Seq-K_pvsmOHT_siCTRL_chr17_45759770-45760603_pvsmOHT"="4C-seq_K_pvsmOHT_siCTRL"
  ,"4C-Seq-K_pvsmOHT_siCTRL_chr17_57168614-57169531_pvsmOHT"="4C-seq_K_pvsmOHT_siCTRL"
  ,"4C-Seq-K_pvsmOHT_siCTRL_chr1_89455867-89456712_pvsmOHT"="4C-seq_K_pvsmOHT_siCTRL"
  ,"4C-Seq-K_pvsmOHT_siCTRL_chr20_30946314-30947710_pvsmOHT"="4C-seq_K_pvsmOHT_siCTRL"
  ,"4C-Seq-K_pvsmOHT_siCTRL_chr21_33251469-33252587_pvsmOHT"="4C-seq_K_pvsmOHT_siCTRL"
  ,"4C-Seq-K_pvsmOHT_siSCC1_chr17_45759770-45760603_pvsmOHT"="4C-seq_K_pvsmOHT_siSCC1"
  ,"4C-Seq-K_pvsmOHT_siSCC1_chr17_57168614-57169531_pvsmOHT"="4C-seq_K_pvsmOHT_siSCC1"
  ,"4C-Seq-K_pvsmOHT_siSCC1_chr1_89455867-89456712_pvsmOHT"="4C-seq_K_pvsmOHT_siSCC1"
  ,"4C-Seq-K_pvsmOHT_siSCC1_chr20_30946314-30947710_pvsmOHT"="4C-seq_K_pvsmOHT_siSCC1"
  ,"4C-Seq-K_pvsmOHT_siSCC1_chr21_33251469-33252587_pvsmOHT"="4C-seq_K_pvsmOHT_siSCC1"
)

path4C <- "/home/rochevin/Documents/PROJET_INGE/4CSeq_PROCESS_VINCENT/AFTER_EVA_PROCESSING/4C_ALN_TRANS/BAMCOMPARE"
wigs.4C <- list.files(path4C,pattern="_normalized.bw")


Get1val <- function(my.wigs,one.w,x){
  lapply(split(x,droplevels(seqnames(x))),function(zz){
    message(unique(as.character(seqnames(zz))))
    cov <- one.w[[unique(as.character(seqnames(zz)))]]
    score <- Views( cov, start = start(zz), end = end(zz) ) %>% sum()
    tibble(wig = my.wigs,value = score,rowname = zz$name)
  }) %>% bind_rows()
}
to.plot.4C.compare <- sres <- lapply(c(100000,500000,1000000),function(bin){
  values <-  bless80 %>% anchor_center() %>% mutate(width=bin)
  lapply(wigs.4C,function(wig){
    one.w <- import.bw(str_c(path4C,wig,sep="/"),as="RleList")
    File <- replaceName.compare[[str_remove(wig,"_normalized.bw")]]
    vp <- str_extract(wig,"chr[0-9A-Z]+_[0-9]+-[0-9]+") %>% str_replace("_",":")
    x <- values %>% filter_by_non_overlaps(GRanges(vp))
    Get1val(File,one.w,x) %>% mutate(binsize = bin) %>% mutate(viewpoint = vp)
  }) %>% bind_rows()
  
}) %>% bind_rows()

to.plot.4C.compare <- to.plot.4C.compare %>%
  mutate(Manip = str_extract(wig,"4C-[sS]eq_[A-Za-z0-9]+_")) %>%
  mutate(Condition = str_remove(wig,"4C-[sS]eq_[A-Za-z0-9]+_"))


filterme <- to.plot.4C.compare %>%
  filter(Manip %in% c("4C-Seq_D_","4C-Seq_Dbis_"),Condition == "pvsmOHT") %>% group_by(binsize,viewpoint,Manip)  %>%
  nest() %>% group_by(binsize) %>% nest() %>%
  mutate(data = purrr::map(data,. %>% distinct(viewpoint,.keep_all = T))) %>%
  unnest() %>%
  mutate(data = purrr::map(data,. %>% filter(value > quantile(value,0.75)))) %>% unnest()

#Barplot dataset bamcompare for all DSBs last manip I et K

dataIK <- to.plot.4C.compare %>% filter(Manip %in% c("4C-seq_I_","4C-seq_K_"))

pdf("BOXPLOT_TRANS/Bamcompare_TRANS_bless80_on4C_seq_data_barplot_bamcompare_value_IK.pdf",height=12,width=16)
for(bin in unique(dataIK$binsize)){
  for(onemanip in unique(dataIK$Manip)){
    subdf <- dataIK %>% filter(binsize == bin,viewpoint %in% c("chr20:30946314-30947710","chr1:89455867-89456712","chr17:45759770-45760603","chr17:57168614-57169531","chr21:33251469-33252587"),Manip == onemanip)
    p <- subdf %>%
      ggplot(aes(x=Condition,y=value,fill = Condition)) +  
      stat_summary(fun.y=mean,position=position_dodge(width=0.95),geom="bar")+
      stat_summary(fun.data=mean_se,position=position_dodge(0.95),geom="errorbar") + 
      coord_flip()+
      facet_grid(viewpoint~Manip,scales="free_x") +
      theme_minimal() + theme(legend.position = "none",axis.text.x = element_text(angle = 90)) 
    print(p)
  }
  
  
}
dev.off()

#Barplot dataset bamcompare for all DSBs 
pdf("BOXPLOT_TRANS/Bamcompare_TRANS_bless80_on4C_seq_data_barplot_bamcompare_value.pdf",height=12,width=16)
for(bin in unique(to.plot.4C.compare$binsize)){
  for(onemanip in unique(to.plot.4C.compare$Manip)){
    subdf <- to.plot.4C.compare %>% filter(binsize == bin,viewpoint %in% c("chr20:30946314-30947710","chr1:89455867-89456712","chr17:45759770-45760603","chr17:57168614-57169531","chr21:33251469-33252587"),Manip == onemanip)
    p <- subdf %>%
      ggplot(aes(x=Condition,y=value,fill = Condition)) +  
      stat_summary(fun.y=mean,position=position_dodge(width=0.95),geom="bar")+
      stat_summary(fun.data=mean_se,position=position_dodge(0.95),geom="errorbar") + 
      coord_flip()+
      facet_grid(viewpoint~Manip,scales="free_x") +
      theme_minimal() + theme(legend.position = "none",axis.text.x = element_text(angle = 90)) +ggtitle(bin)
    print(p)
  }
  
  
}
dev.off()
#Update au 28/09/2020 pvals de la familia pour le papier de Coline
ccpval <- to.plot.4C.compare %>% filter(viewpoint %in% c("chr20:30946314-30947710","chr1:89455867-89456712","chr17:45759770-45760603","chr17:57168614-57169531","chr21:33251469-33252587"),Manip == "4C-seq_F_",binsize == 1e+06)
ccpval %>% dplyr::select(-wig) %>% spread(key = Condition,value=value)  %>%
  group_by(viewpoint) %>% nest() %>%
  mutate(pvalG1vsG2 = map_dbl(data,function(x){
    wilcox.test(x$pvsmOHT_G1,x$pvsmOHT_G2,paired=T)$p.value
  })) %>% 
  mutate(pvalG1vsS = map_dbl(data,function(x){
    wilcox.test(x$pvsmOHT_G1,x$pvsmOHT_S,paired=T)$p.value
  })) %>% 
  mutate(pvalG2vsS = map_dbl(data,function(x){
    wilcox.test(x$pvsmOHT_G2,x$pvsmOHT_S,paired=T)$p.value
  })) %>% dplyr::select(-data)

ccpval <- to.plot.4C.compare %>% filter(viewpoint %in% c("chr20:30946314-30947710","chr1:89455867-89456712","chr17:45759770-45760603","chr17:57168614-57169531","chr21:33251469-33252587"),Manip == "4C-Seq_D_",binsize == 1e+06)
ccpval %>% dplyr::select(-wig) %>% spread(key = Condition,value=value)  %>%
  group_by(viewpoint) %>% nest() %>%
  mutate(pDNAPKivspOHT = map_dbl(data,function(x){
    wilcox.test(x$pDNAPKivsmOHT,x$pvsmOHT,paired=T)$p.value
  })) %>% dplyr::select(-data)

ccpval %>% dplyr::select(-wig) %>% spread(key = Condition,value=value)  %>%
  group_by(viewpoint) %>% nest() %>%
  mutate(pDNAPKivspOHT = map_dbl(data,function(x){
    wilcox.test(x$pDNAPKivsmOHT,x$pvsmOHT)$p.value
  }))%>% dplyr::select(-data)

#Barplot dataset bamcompare for UP manip D/Dbis DSBs 
pdf("BOXPLOT_TRANS/Bamcompare_TRANS_bless80_on4C_seq_data_barplot_bamcompare_value_75quant_dsb_from_manipD.pdf",height=12,width=16)
for(bin in unique(to.plot.4C.compare$binsize)){
  for(onemanip in unique(to.plot.4C.compare$Manip)){
    subdf <- to.plot.4C.compare %>% filter(binsize == bin,viewpoint %in% c("chr20:30946314-30947710","chr1:89455867-89456712","chr17:45759770-45760603","chr17:57168614-57169531","chr21:33251469-33252587"),Manip == onemanip)
    subdf <- subdf %>% semi_join(filter(filterme,binsize == bin),by = c("rowname","viewpoint"))
    p <- subdf %>%
      ggplot(aes(x=Condition,y=value,fill = Condition)) +  
      stat_summary(fun.y=mean,position=position_dodge(width=0.95),geom="bar")+
      stat_summary(fun.data=mean_se,position=position_dodge(0.95),geom="errorbar") + 
      coord_flip()+
      facet_grid(viewpoint~Manip,scales="free_x") +
      theme_minimal() + theme(legend.position = "none",axis.text.x = element_text(angle = 90)) 
    print(p)
  }
  
  
}
dev.off()
