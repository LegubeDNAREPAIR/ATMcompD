

mes_bw <- c("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_HLNKYBGXC_SCC1//PROCESSED/mapping/EXPERIMENT/BIGWIG/HLNKYBGXC_Pool_ChIP-seq_legube_19s004478-1-1_Clouaire_lane1Rad21DIvA_sequence.exp_spikeinfactor.bw",
            "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_HLNKYBGXC_SCC1//PROCESSED/mapping/EXPERIMENT/BIGWIG/HLNKYBGXC_Pool_ChIP-seq_legube_19s004478-1-1_Clouaire_lane1Rad21OHT_sequence.exp_spikeinfactor.bw")
names(mes_bw) <- basename(mes_bw) %>% 
  str_extract("SMC.+_sequence|Rad21DIvA|Rad21OHT|[A-za-z0-9]+SCC1.+_sequence") %>% str_remove("_sequence") %>% str_remove("1_Clouaire_lane1")
mes_bw <- mes_bw %>% map(import.bw,as="RleList")


bless80 = read_bed("/mnt/NAS/DATA/AsiSI/BLESS_80best_JunFragPE_Rmdups_pm500bp.bed")
# bless80 = read_bed("/mnt/NAS/DATA/AsiSI/ASIsites_hg19_174_clived_IQ1.5.bed")
my_window <- 1000000
mes_seqlens <- seqlens %>% enframe()
mon_bed <- bless80 %>%
  anchor_center() %>% mutate(width = my_window*2) 

#load rad21 peaks
rad21_peaks <- read_narrowpeaks("/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/ChIP-Seq/Clouaire_HLNKYBGXC_SCC1/PROCESSED/mapping/EXPERIMENT/MACS/Rad21DIvA.exp_peaks.narrowPeak")
rad21_peaks_inside <- rad21_peaks %>% join_overlap_inner(mon_bed) 
rad21_peaks_inside$name <- rad21_peaks_inside$name.x
rad21_peaks_outside <- rad21_peaks %>% filter(!name %in% rad21_peaks_inside$name.x)


mes_peaks_list <- list(
  "rad21.peaks_damaged" = rad21_peaks_inside,
  "rad21.peaks_undamaged" = rad21_peaks_outside
)


res.boxplot <- lapply(names(mes_bw),function(one.w){
  message(one.w)
  lapply(mes_peaks_list,function(one_bed){
    PhDfunc::Get1val(one.w,mes_bw[[one.w]],one_bed) %>% left_join(as_tibble(one_bed),by=c("rowname"="name"))
  }) %>% bind_rows(.id="Type")
}) %>% bind_rows() %>% 
  mutate(Condition = str_extract(wig,"DIvA|OHT")) %>% 
  mutate(wig = str_remove(wig,"DIvA|OHT")) %>% 
  separate("Type",into = c("peak","Damaged"),sep="_")

#boxplot sur loop anchor et sur rad21
DSB_n_chr <- str_c("chr",c(1:22,"X")) %>% enframe(value = "seqnames") %>% dplyr::select(-name) %>% left_join(bless80 %>% as_tibble() %>% dplyr::count(seqnames)) %>% mutate(n = replace_na(n,0)) %>% arrange(n) 
peak_rad21_by_chr <- res.boxplot %>% filter(peak == "rad21.peaks",Damaged == "undamaged") %>% 
  left_join(DSB_n_chr,by = "seqnames") %>% 
  ggplot(aes(x=as.factor(n),y=value,fill=Condition)) +
  scale_y_log10() +
  geom_boxplot() + facet_wrap(~wig,scales="free_y",ncol=1)
